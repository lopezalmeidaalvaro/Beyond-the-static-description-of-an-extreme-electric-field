"""Sistema ACOPLADO completo: cascaron cargado dinamico + campo electrico radial + plasma (dos flujos nulos) + descarga local, en calibre polar-areal.
El cascaron es MATERIA DINAMICA del mismo solver (sin masa de prueba): una particula-cascaron de masa en reposo mu~=2G mu (=C), momento radial ortonormal w, carga q,
que entra en los vinculos de Einstein con su energia E_d = mu~ sqrt(1+w^2)/a, su flujo de momento P_d = E_d w^2/(1+w^2) y su flujo de energia E_d dR/dt:
    m(r)  : salto E_d (reparto CIC en dos celdas),            phi_r += (1/2) a^2 P_d/(r dr)  (salto de lapso por la tension radial del cascaron en movimiento),
    lam_t += -(1/2) a^2 E_d Rdot/(r dr)                       (T^r_t del cascaron);
su ecuacion de movimiento es la geodesica radial en la metrica externa (sin su propio delta) mas la aceleracion de Coulomb propia a_p = q^2 g^2/(C R^2):
    dR/dt = alpha w/(a sqrt(1+w^2)),   dw/dt = -alpha phi_r sqrt(1+w^2)/a - lam_t w + alpha a_p.
El campo exterior (carga sobre el cascaron, interior sin campo) es rho~_EM = q^2 g(s_c)^2/r^4 para r>R(t); el plasma lo atraviesa libremente y
se crea 1/2+1/2 por conversion local (mismo cierre que dyn1d).  Unidades r0=1, c=1, compacidad m=2GM."""
import numpy as np
from dyn1d import _g_of, _minmod


def run_shell(q, C, s_par, curve, kappa, Nr=2000, Rmax=16.0, Tend=14.0, r_core=0.0, instant=False, plasma=True, no_coulomb=False, cfl=0.4, ah_stop=0.999,
              R_escape=7.0, record=True, eom_scale=1.0, plasma_gravity=True, freeze=False):
    dr = Rmax/Nr
    r = (np.arange(Nr) + 0.5)*dr
    jc = int(round(r_core/dr))
    s0 = curve[0][-1]; Lam = kappa*s_par
    sc = np.zeros(Nr)
    Wp = np.zeros(Nr); Wm = np.zeros(Nr)
    Rs, ws, tau_s, Rdot = 1.0, 0.0, 0.0, 0.0
    t = 0.0
    E_out_b = 0.0; E_in_b = 0.0
    out = dict(status='ok', t_AH=None, r_AH=None)
    rec = dict(t=[], R=[], w=[], tau=[], Ed=[], Cplus=[], mtot=[], Cmax=[], E_out=[], E_in=[])

    def build(Rs_, ws_, Rdot_):
        frac = np.clip((r + 0.5*dr - Rs_)/dr, 0.0, 1.0)
        if instant: g = np.zeros(Nr)
        else: g = _g_of(sc, curve)
        rho_geom = frac*q**2/r**4
        rho_em = rho_geom*g**2
        rho_pl = (Wp + Wm)/r**2 if plasma else np.zeros(Nr)
        w0 = r**2*(rho_pl + rho_em)*dr
        mf0 = np.concatenate(([0.0], np.cumsum(w0)))
        m0c = 0.5*(mf0[:-1] + mf0[1:])
        m_minus = float(np.interp(Rs_, r, m0c))
        gam = np.sqrt(1.0 + ws_*ws_)
        E = C*gam
        for _ in range(8):
            E = C*gam*np.sqrt(max(1.0 - (m_minus + 0.5*E)/Rs_, 1e-6))
        x = (Rs_ - r[0])/dr; j0 = int(np.floor(x)); fr_ = x - j0
        j0 = min(max(j0, 0), Nr - 2)
        wsh = np.zeros(Nr); Psh = np.zeros(Nr); fsh = np.zeros(Nr)
        V2 = ws_*ws_/(1.0 + ws_*ws_)
        for jj, wt in ((j0, 1.0 - fr_), (j0 + 1, fr_)):
            wsh[jj] = E*wt; Psh[jj] = E*V2*wt; fsh[jj] = E*Rdot_*wt
        w = w0 + wsh
        mface = np.concatenate(([0.0], np.cumsum(w)))
        m_c = 0.5*(mface[:-1] + mface[1:])
        a2 = 1.0/np.maximum(1.0 - m_c/r, 1e-6)
        m_r = r**2*(rho_pl + rho_em) + wsh/dr
        p_r = rho_pl - rho_em
        phi_r_ext = 0.5*a2*r*p_r + (a2 - 1.0)/(2.0*r)
        phi_r = phi_r_ext + 0.5*a2*Psh/(r*dr)
        lam_r = a2*(m_r/(2.0*r) - m_c/(2.0*r**2))
        phi_N = -0.5*np.log(a2[-1])
        seg = 0.5*(phi_r[:-1] + phi_r[1:])*dr
        phi = np.empty(Nr); phi[-1] = phi_N; phi[:-1] = phi_N - np.cumsum(seg[::-1])[::-1]
        alpha = np.exp(phi); v = alpha/np.sqrt(a2)
        lam_t_pl = -0.5*a2*v*(Wp - Wm)/r if plasma else np.zeros(Nr)
        lam_t = lam_t_pl - 0.5*a2*fsh/(r*dr)
        return dict(m_c=m_c, a2=a2, phi_r=phi_r, phi_r_ext=phi_r_ext, lam_r=lam_r, v=v, lam_t=lam_t, lam_t_ext=lam_t_pl, alpha=alpha, g=g, frac=frac, rho_geom=rho_geom, E=E)

    def rhs(Wp_, Wm_, M):
        v = M['v']; phi_r = M['phi_r']; lam_r = M['lam_r']; lam_t = M['lam_t']
        n = Nr - jc
        Wa_p = Wp_[jc:]; Wa_m = Wm_[jc:]; va = v[jc:]
        gl_p = Wa_m[1::-1]; gl_m = Wa_p[1::-1]
        def faces(W, direction, ghost, right):
            Wg = np.concatenate((ghost, W, [right, right]))
            vg = np.concatenate((va[1::-1], va, [va[-1], va[-1]]))
            d = np.diff(Wg)
            sl = np.zeros_like(Wg); sl[1:-1] = _minmod(d[:-1], d[1:])
            Wf = Wg[:-1] + 0.5*sl[:-1] if direction > 0 else Wg[1:] - 0.5*sl[1:]
            vf = 0.5*(vg[:-1] + vg[1:])
            return (direction*vf*Wf)[1:n+2]
        Fpa = faces(Wa_p, +1, gl_p, Wa_p[-1]); Fma = faces(Wa_m, -1, gl_m, 0.0)
        Fpa[0] = -Fma[0]
        Fp = np.zeros(Nr + 1); Fm = np.zeros(Nr + 1); Fp[jc:] = Fpa; Fm[jc:] = Fma
        coef = v*(phi_r + lam_r)
        dWp = -(Fp[1:] - Fp[:-1])/dr - Wp_*(coef + 2.0*lam_t)
        dWm = -(Fm[1:] - Fm[:-1])/dr + Wm_*(coef - 2.0*lam_t)
        if jc: dWp[:jc] = 0.0; dWm[:jc] = 0.0
        return dWp, dWm, Fp, Fm

    def convert(dt_c, M):
        nonlocal sc, Wp, Wm
        if instant or not plasma: return
        mask = M['frac'] > 0
        mask[:jc] = False
        sc_new = np.where(mask, sc + dt_c*M['alpha']*Lam/r, sc)           # E/E_surf = 1/r^2 en el exterior del cascaron
        g_old = _g_of(sc, curve); g_new = _g_of(sc_new, curve)
        d_em = M['rho_geom']*(g_old**2 - g_new**2)
        sc = sc_new
        add = 0.5*r**2*d_em
        Wp = Wp + add; Wm = Wm + add

    def sdot(R_, w_, M, Rdot_unused=None):
        aa = np.sqrt(np.interp(R_, r, M['a2'])); al = np.interp(R_, r, M['alpha'])
        fr = np.interp(R_, r, M['phi_r_ext']); lt = np.interp(R_, r, M['lam_t_ext'])
        gR = float(np.interp(R_ + 0.5*dr, r, M['g']))
        ap = 0.0 if no_coulomb else q**2*gR**2/(C*R_**2)
        gm = np.sqrt(1.0 + w_*w_)
        return al*w_/(aa*gm), eom_scale*(-al*fr*gm/aa) - lt*w_ + al*ap, al/gm

    exposed = (np.clip((r + 0.5*dr - Rs)/dr, 0.0, 1.0) > 0)
    nstep = 0
    while t < Tend:
        M = build(Rs, ws, Rdot)
        vmax = float(np.max(M['v']))
        dt = cfl*dr/vmax
        if t + dt > Tend: dt = Tend - t
        lim = max(Rmax - t, 2.0)
        mask = r < lim
        C_ = M['m_c']/r
        if np.nanmax(C_[mask]) >= ah_stop:
            j = int(np.argmax(np.where(mask, C_, -1)))
            out.update(status='AH', t_AH=t, r_AH=float(r[j]), R_shell=float(Rs), tau_AH=float(tau_s)); break
        if Rs >= R_escape and ws > 0:
            out.update(status='escaped', t_esc=t); break
        if Rs < max(r_core, 3*dr) + dr:
            out.update(status='inner', t_in=t); break
        convert(0.5*dt, M)
        M = build(Rs, ws, Rdot); Ma = M
        if plasma:
            dWp1, dWm1, Fp1, Fm1 = rhs(Wp, Wm, M)
            Wp_s = Wp + dt*dWp1; Wm_s = Wm + dt*dWm1
            sv = (Wp, Wm); Wp, Wm = Wp_s, Wm_s
            M2 = build(Rs, ws, Rdot)
            dWp2, dWm2, Fp2, Fm2 = rhs(Wp, Wm, M2)
            Wp = 0.5*(sv[0] + Wp_s + dt*dWp2); Wm = 0.5*(sv[1] + Wm_s + dt*dWm2)
            Fpf = 0.5*(Fp1 + Fp2); Fmf = 0.5*(Fm1 + Fm2)
            E_out_b += dt*float(Fpf[Nr]); E_in_b += dt*float(-Fmf[Nr])
        Mb = build(Rs, ws, Rdot)
        k1 = sdot(Rs, ws, Ma); k2 = sdot(Rs + dt*k1[0], ws + dt*k1[1], Mb)
        Rdot = 0.5*(k1[0] + k2[0])
        if freeze: k1 = (0.0, 0.0, k1[2]); k2 = (0.0, 0.0, k2[2]); Rdot = 0.0
        Rs += dt*Rdot; ws += 0.5*dt*(k1[1] + k2[1]); tau_s += 0.5*dt*(k1[2] + k2[2])
        fr_new = np.clip((r + 0.5*dr - Rs)/dr, 0.0, 1.0) > 0
        newly = fr_new & ~exposed
        if newly.any():                                     # celdas que el cascaron deja atras al contraerse: el campo exterior hereda el apantallamiento de la celda vecina (continuidad de Gauss)
            jref = int(np.max(np.where(newly)[0])) + 1
            sc[newly] = sc[jref]
        exposed = fr_new
        Mc = build(Rs, ws, Rdot)
        convert(0.5*dt, Mc)
        t += dt; nstep += 1
        if record and nstep % 4 == 0:
            Mr = build(Rs, ws, Rdot)
            rec['t'].append(t); rec['R'].append(Rs); rec['w'].append(ws); rec['tau'].append(tau_s); rec['Ed'].append(Mr['E'])
            jp = min(int(Rs/dr) + 2, Nr - 1)
            rec['Cplus'].append(float(Mr['m_c'][jp]/r[jp])); rec['mtot'].append(float(Mr['m_c'][-1]))
            rec['Cmax'].append(float(np.max((Mr['m_c']/r)[r < max(Rmax - t, 1.0)]))); rec['E_out'].append(E_out_b); rec['E_in'].append(E_in_b)
    out.update(t_end=t, nstep=nstep, R_end=float(Rs), w_end=float(ws), tau_end=float(tau_s), E_out_b=E_out_b, E_in_b=E_in_b, m_outer=float(build(Rs, ws, Rdot)['m_c'][-1]))
    if record:
        for k in rec: rec[k] = np.array(rec[k])
        out['rec'] = rec
    return out
