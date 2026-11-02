"""Motor 1+1D (calibre polar-areal) para la descarga de Schwinger no homogenea con fuga/inflow radiativo en Relatividad General.
Unidades: r0=1 (radio de la bola cargada), c=1, y todas las energias en compacidad  m=2G M  (longitud);  rho~ = 8 pi G rho.
Metrica  ds^2 = -alpha^2 dt^2 + a^2 dr^2 + r^2 dOmega^2,  a^-2 = 1 - m/r,  m_r = r^2 rho~_tot  (vinculo hamiltoniano, derivado en 06b: D_MS1).
Contenido:  - campo electrico radial de una bola uniforme: rho~_geom = q^2 r^2 (r<1),  q^2/r^4 (r>1),  q^2 = r_Q^2/r0^2 = 2 G U_out / r0
            - dos flujos nulos radiales (plasma colineal, w_par ~ 1: resultado de 05c): W_+ saliente, W_- entrante, W = r^2 rho~
              d_t W_+ + d_r(v W_+) = -W_+ [v(phi_r+lam_r) + 2 lam_t],   d_t W_- - d_r(v W_-) = W_- [v(phi_r+lam_r) - 2 lam_t]    (D_MS4)
              v = alpha/a;  lam_t = -(1/2) a^2 v (W_+ - W_-)/r;  phi_r = (a^2/2) r p~_r + (a^2-1)/(2r),  p~_r = rho~_+ + rho~_- - rho~_EM   (D_MS1/D_MS4)
            - descarga local: reloj s_c' = alpha*Lambda*sqrt(E/E_surf), Lambda = kappa*s_par, s_par=sqrt(N/N*);  rho~_EM = rho~_geom g(s_c)^2
              con g(s) = curva universal EXACTA de la primera oscilacion (05c); la energia convertida se reparte 1/2 saliente + 1/2 entrante (momento neto nulo)
            - masa 'mantenida' M_n: offset estatico de prueba m_mat = C_n min(1,r^3)  (aproximacion de masa de prueba declarada; ver gate D_MS5b)
Opciones: G_on=False (limite plano), absorb (materia opaca), instant (descarga a t=0), interior/exterior (canales de conversion)."""
import numpy as np


def _g_of(sc, curve):
    s, g = curve
    return np.where(sc >= s[-1], 0.0, np.interp(sc, s, g))


def metric(r, dr, W_p, W_m, rho_em, m_mat_c, mr_mat, G_on=True, rho_heat=None, alpha_outer_from_a=True):
    """m, a^2, phi_r, lam_r, v, lam_t en centros de celda."""
    rho_pl = (W_p + W_m)/r**2
    rho_tot = rho_pl + rho_em + (0.0 if rho_heat is None else rho_heat)
    w = r**2*rho_tot*dr
    mface = np.concatenate(([0.0], np.cumsum(w)))
    m_c = 0.5*(mface[:-1] + mface[1:]) + m_mat_c
    if not G_on:
        z = np.zeros_like(r)
        return m_c, np.ones_like(r), z, z, np.ones_like(r), z, np.ones_like(r)
    a2 = 1.0/np.maximum(1.0 - m_c/r, 1e-6)
    m_r = r**2*rho_tot + mr_mat
    p_r = rho_pl - rho_em
    phi_r = 0.5*a2*r*p_r + (a2 - 1.0)/(2.0*r)
    lam_r = a2*(m_r/(2.0*r) - m_c/(2.0*r**2))
    phi = np.empty_like(r)
    phi[-1] = -0.5*np.log(a2[-1])
    for j in range(len(r) - 2, -1, -1):
        phi[j] = phi[j+1] - 0.5*(phi_r[j] + phi_r[j+1])*dr
    alpha = np.exp(phi)
    v = alpha/np.sqrt(a2)
    lam_t = -0.5*a2*v*(W_p - W_m)/r
    return m_c, a2, phi_r, lam_r, v, lam_t, alpha


def _minmod(a, b):
    return np.where(a*b > 0, np.where(np.abs(a) < np.abs(b), a, b), 0.0)


def run_dyn(q, s_par, curve, kappa, Cn=0.0, Nr=900, Rmax=6.0, Tend=2.5, G_on=True, absorb=False, instant=False,
            interior=True, exterior=True, cfl=0.4, ah_stop=0.999, record=True, no_balance=0.0, center_reflect=True, r_core=0.0, particle=False, mass_profile='ball'):
    """Evoluciona y devuelve diagnosticos. no_balance=delta>0: control negativo (el plasma recibe (1+delta) de lo que pierde el campo)."""
    dr = Rmax/Nr
    r = (np.arange(Nr) + 0.5)*dr
    iface = int(round(1.0/dr))                         # indice de la cara r=1: entre celdas iface-1 e iface
    assert abs(iface*dr - 1.0) < 1e-12, "r=1 debe ser una cara de la malla"
    rho_geom = np.where(r < 1.0, q**2*r**2, q**2/r**4)
    Er = np.where(r < 1.0, r, 1.0/r**2)
    Er_sqrt = np.sqrt(Er)
    s0 = curve[0][-1]
    Lam = kappa*s_par
    if mass_profile == 'point':
        m_mat_c = Cn*np.ones_like(r); mr_mat = np.zeros_like(r)
    else:
        m_mat_c = Cn*np.minimum(1.0, r**3)
        mr_mat = np.where(r < 1.0, 3.0*Cn*r**2, 0.0)
    jc = int(round(r_core/dr))
    inside = r < 1.0
    sc = np.zeros(Nr)
    conv_mask = np.ones(Nr, bool)
    if not interior: conv_mask &= ~inside
    if not exterior: conv_mask &= inside
    conv_mask[:jc] = False
    Wp = np.zeros(Nr); Wm = np.zeros(Nr); Wh = np.zeros(Nr)      # Wh: calor depositado (materia opaca)
    if instant:
        sc = np.where(conv_mask, 2*s0, 0.0)
        e0 = r**2*rho_geom*np.where(conv_mask, 1.0, 0.0)
        Wp = 0.5*e0*(1.0 + no_balance); Wm = 0.5*e0*(1.0 + no_balance)
        if absorb: Wh = np.where(inside, Wp + Wm, 0.0); Wp = np.where(inside, 0.0, Wp); Wm = np.where(inside, 0.0, Wm)
    t = 0.0
    rec = dict(t=[], dm_face=[], Cb=[], Cmax=[], rmax=[], Fin=[], Fout=[], Ein=[], Eout=[], mface1=[], residual=[])
    E_in = 0.0; E_out = 0.0; out = dict(status='ok', t_AH=None, r_AH=None, t_b=None)

    def state_metric():
        rho_em = rho_geom*_g_of(sc, curve)**2
        return metric(r, dr, Wp, Wm, rho_em, m_mat_c, mr_mat, G_on, rho_heat=Wh/r**2), rho_em

    def rhs(Wp_, Wm_, mt):
        (m_c, a2, phi_r, lam_r, v, lam_t, alpha), rho_em = mt
        n = Nr - jc
        Wa_p = Wp_[jc:]; Wa_m = Wm_[jc:]; va = v[jc:]
        gl_p = Wa_m[1::-1]; gl_m = Wa_p[1::-1]            # fantasma izquierdo: reflexion en el centro (jc=0) o en el nucleo r=r_core (jc>0)
        def faces(W, direction, ghost, right):
            Wg = np.concatenate((ghost, W, [right, right]))
            vg = np.concatenate((va[1::-1], va, [va[-1], va[-1]]))
            d = np.diff(Wg)
            sl = np.zeros_like(Wg); sl[1:-1] = _minmod(d[:-1], d[1:])
            Wf = Wg[:-1] + 0.5*sl[:-1] if direction > 0 else Wg[1:] - 0.5*sl[1:]
            vf = 0.5*(vg[:-1] + vg[1:])
            return (direction*vf*Wf)[1:n+2]
        Fpa = faces(Wa_p, +1, gl_p, Wa_p[-1]); Fma = faces(Wa_m, -1, gl_m, 0.0)
        Fpa[0] = -Fma[0] if center_reflect else 0.0       # el flujo entrante que sale por el centro/nucleo reaparece como saliente
        Fp = np.zeros(Nr + 1); Fm = np.zeros(Nr + 1); Fp[jc:] = Fpa; Fm[jc:] = Fma
        coef = v*(phi_r + lam_r)
        dWp = -(Fp[1:] - Fp[:-1])/dr - Wp_*(coef + 2.0*lam_t)
        dWm = -(Fm[1:] - Fm[:-1])/dr + Wm_*(coef - 2.0*lam_t)
        if jc: dWp[:jc] = 0.0; dWm[:jc] = 0.0
        return dWp, dWm, Fp, Fm, v

    def convert(dt_c, mt):
        nonlocal sc, Wp, Wm, Wh
        (m_c, a2, phi_r, lam_r, v, lam_t, alpha), rho_em = mt
        if instant: return
        sc_new = np.where(conv_mask, sc + dt_c*alpha*Lam*Er_sqrt, sc)
        d_em = rho_geom*(_g_of(sc, curve)**2 - _g_of(sc_new, curve)**2)
        sc = sc_new
        add = 0.5*r**2*d_em*(1.0 + no_balance)
        if absorb:
            Wh = Wh + np.where(inside, 2*add, 0.0); add = np.where(inside, 0.0, add)
        Wp = Wp + add; Wm = Wm + add

    rp, wp, tau_p = 1.0, 0.0, 0.0; p_active = bool(particle and G_on); t_ptrap = None; tau_ptrap = None
    rec.update(rp=[], Cp=[], tau=[], wp=[])
    nstep = 0
    while t < Tend:
        mt = state_metric()
        (m_c, a2, phi_r, lam_r, v, lam_t, alpha), rho_em = mt
        vmax = float(np.max(v)) if G_on else 1.0
        dt = cfl*dr/vmax
        if t + dt > Tend: dt = Tend - t
        # diagnosticos (antes del paso)
        lim = Rmax - t
        mask = r < lim
        C = m_c/r
        if G_on and np.nanmax(C[mask]) >= ah_stop:
            j = int(np.argmax(np.where(mask, C, -1)))
            out.update(status='AH', t_AH=t, r_AH=float(r[j]))
            break
        if G_on and out['t_b'] is None and float(np.interp(1.0, r, C)) >= 1.0 - 1e-9:
            out['t_b'] = t
        # --- Strang: conv(dt/2) transporte(dt) conv(dt/2)
        convert(0.5*dt, mt)
        mt = state_metric(); mt_a = mt
        dWp1, dWm1, Fp1, Fm1, v1 = rhs(Wp, Wm, mt)
        Wp_s = Wp + dt*dWp1; Wm_s = Wm + dt*dWm1
        # etapa 2 (Heun)
        sv = (Wp, Wm); Wp, Wm = Wp_s, Wm_s
        mt2 = state_metric()
        dWp2, dWm2, Fp2, Fm2, v2 = rhs(Wp, Wm, mt2)
        Wp = 0.5*(sv[0] + Wp_s + dt*dWp2); Wm = 0.5*(sv[1] + Wm_s + dt*dWm2)
        Fp_f = 0.5*(Fp1 + Fp2); Fm_f = 0.5*(Fm1 + Fm2)
        # flujo a traves de r=1 (cara iface): saliente (Fp) y entrante (|Fm|)
        fo = float(Fp_f[iface]); fi = float(-Fm_f[iface])
        E_out += dt*fo; E_in += dt*fi
        if absorb:                                                    # la materia opaca absorbe el plasma que entra en r<1
            Wh = Wh + np.where(inside, Wp + Wm, 0.0); Wp = np.where(inside, 0.0, Wp); Wm = np.where(inside, 0.0, Wm)
        mt = state_metric()
        if p_active:                                               # elemento de materia neutro (geodesica radial) en la metrica dinamica; Heun con metrica inicial/final
            def pdot(rp_, w_, M):
                (m_c_, a2_, phi_r_, lam_r_, v_, lam_t_, al_) = M[0]
                aa = np.sqrt(np.interp(rp_, r, a2_)); al = np.interp(rp_, r, al_); fr = np.interp(rp_, r, phi_r_); lt = np.interp(rp_, r, lam_t_)
                g_ = np.sqrt(1.0 + w_*w_)
                return al*w_/(aa*g_), -al*fr*g_/aa - lt*w_, al/g_
            k1 = pdot(rp, wp, mt_a); k2 = pdot(rp + dt*k1[0], wp + dt*k1[1], mt)
            rp += 0.5*dt*(k1[0] + k2[0]); wp += 0.5*dt*(k1[1] + k2[1]); tau_p += 0.5*dt*(k1[2] + k2[2])
            Cp_ = float(np.interp(rp, r, mt[0][0]/r))
            if Cp_ >= 0.999 and t_ptrap is None: t_ptrap = t + dt; tau_ptrap = tau_p
            if rp < max(r[jc], 2*dr): p_active = False
        convert(0.5*dt, mt)
        t += dt; nstep += 1
        if record:
            mt = state_metric(); m_c2 = mt[0][0]; Cc = m_c2/r
            m1 = float(np.interp(1.0, r, m_c2))
            rec['t'].append(t); rec['mface1'].append(m1); rec['Cb'].append(m1); rec['Ein'].append(E_in); rec['Eout'].append(E_out)
            rec['Cmax'].append(float(np.max(Cc[r < Rmax - t]))); rec['rmax'].append(float(r[np.argmax(np.where(r < Rmax - t, Cc, -1))]))
            rec['Fin'].append(fi); rec['Fout'].append(fo)
            rec['rp'].append(rp); rec['wp'].append(wp); rec['tau'].append(tau_p); rec['Cp'].append(float(np.interp(rp, r, mt[0][0]/r)))
    out.update(t_end=t, nstep=nstep, E_in=E_in, E_out=E_out, Nr=Nr, particle_t_trap=t_ptrap, particle_tau_trap=tau_ptrap, particle_r=rp, particle_tau=tau_p)
    if record:
        for k in rec: rec[k] = np.array(rec[k])
        out['rec'] = rec
    mt = state_metric(); out['m1_final'] = float(np.interp(1.0, r, mt[0][0]))
    out['m1_initial'] = float(Cn + q**2/5.0)                       # m(r=1) inicial: masa mantenida + energia del campo interior
    out['r'] = r
    return out
