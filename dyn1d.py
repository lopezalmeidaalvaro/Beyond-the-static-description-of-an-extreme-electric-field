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
            interior=True, exterior=True, cfl=0.4, ah_stop=0.999, record=True, no_balance=0.0, center_reflect=True):
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
    m_mat_c = Cn*np.minimum(1.0, r**3)
    mr_mat = np.where(r < 1.0, 3.0*Cn*r**2, 0.0)
    inside = r < 1.0
    sc = np.zeros(Nr)
    conv_mask = np.ones(Nr, bool)
    if not interior: conv_mask &= ~inside
    if not exterior: conv_mask &= inside
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
        # fantasma izquierdo (centro): W_+(-r)=W_-(r), W_-(-r)=W_+(r)
        gl_p = Wm_[1::-1]; gl_m = Wp_[1::-1]
        # flujos conservativos con MUSCL
        def faces(W, direction, ghost, right):
            Wg = np.concatenate((ghost, W, [right, right]))
            vg = np.concatenate((v[1::-1], v, [v[-1], v[-1]]))
            d = np.diff(Wg)
            sl = np.zeros_like(Wg); sl[1:-1] = _minmod(d[:-1], d[1:])
            if direction > 0:
                Wf = Wg[:-1] + 0.5*sl[:-1]
            else:
                Wf = Wg[1:] - 0.5*sl[1:]
            vf = 0.5*(vg[:-1] + vg[1:])
            F = direction*vf*Wf                      # caras entre Wg[k] y Wg[k+1]; Wg[2] es la celda 0
            return F[1:Nr+2]                         # caras: izquierda de celda 0 .. derecha de celda Nr-1
        Fp = faces(Wp_, +1, gl_p, Wp_[-1]); Fm = faces(Wm_, -1, gl_m, 0.0)
        Fp[0] = -Fm[0] if center_reflect else 0.0       # el flujo entrante que sale por r=0 reaparece como saliente (paso por el centro)
        coef = v*(phi_r + lam_r)
        dWp = -(Fp[1:] - Fp[:-1])/dr - Wp_*(coef + 2.0*lam_t)
        dWm = -(Fm[1:] - Fm[:-1])/dr + Wm_*(coef - 2.0*lam_t)
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
        mt = state_metric()
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
        convert(0.5*dt, mt)
        t += dt; nstep += 1
        if record:
            mt = state_metric(); m_c2 = mt[0][0]; Cc = m_c2/r
            m1 = float(np.interp(1.0, r, m_c2))
            rec['t'].append(t); rec['mface1'].append(m1); rec['Cb'].append(m1); rec['Ein'].append(E_in); rec['Eout'].append(E_out)
            rec['Cmax'].append(float(np.max(Cc[r < Rmax - t]))); rec['rmax'].append(float(r[np.argmax(np.where(r < Rmax - t, Cc, -1))]))
            rec['Fin'].append(fi); rec['Fout'].append(fo)
    out.update(t_end=t, nstep=nstep, E_in=E_in, E_out=E_out, Nr=Nr)
    if record:
        for k in rec: rec[k] = np.array(rec[k])
        out['rec'] = rec
    mt = state_metric(); out['m1_final'] = float(np.interp(1.0, r, mt[0][0]))
    out['m1_initial'] = float(Cn + q**2/5.0)                       # m(r=1) inicial: masa mantenida + energia del campo interior
    out['r'] = r
    return out
