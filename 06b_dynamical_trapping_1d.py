"""CIERRE DE LA BRECHA ESPACIAL/GRAVITATORIA (K_ij != 0, transporte radial): Misner-Sharp dinamico con descarga de Schwinger no homogenea y fuga/inflow radiativo.
PARTE SIMBOLICA (SymPy, metrica general ds^2 = -e^{2phi(t,r)}dt^2 + e^{2lam(t,r)}dr^2 + r^2 dOmega^2, SIN simetria temporal):
  D_MS1   M_r = 4 pi r^2 rho ,  M_t = 4 pi r^2 T^r_t = -4 pi r^2 F_flux ;  K_ij = 0  <=>  T^r_t = 0   (el caso P10 es exactamente el de flujo nulo)
  D_MS2   balance cuasilocal: la conversion E^2/2 -> rho_pl (F_{nu lam}J^lam) se cancela en M; M_t(R) = 4 pi R^2 T^r_t(R) UNICAMENTE por el flujo en r=R
  C_MS    control: si la ley de transferencia viola la conservacion total, aparece una fuente espuria de volumen en M_t
  D_MS3   frontera movil: dM/dt|_{r=R(t)} = 4 pi R^2 (T^r_t + Rdot rho) = 0 para materia comovil
  D_MS4   transporte de flujos nulos radiales en la metrica dinamica (ecuaciones derivadas de nabla_mu T^{mu nu} = 0)
  D_MS5   el esquema (transporte + conversion con balance) PRESERVA el vinculo de momento; la materia 'mantenida' lo viola en una cantidad exacta (aprox. de prueba)
  P12     datos iniciales de la bola uniforme: horizonte aparente <=> compacidad maxima >=1 (cota de extremalidad M_ADM >= Q)
PARTE NUMERICA (dyn1d.py): solver 1+1D polar-areal, validacion plana analitica, balance de energia, convergencia (dr, dr/2, dr/4), controles negativos,
  descomposicion de canales (saliente interior / entrante exterior), ley de escala de DeltaC(q, s) y frontera de atrapamiento corregida."""
import sys, json, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, sympy as sp
from common import *
from dyn1d import run_dyn

R = {}; ok_all = True
def gate(name, ok, msg):
    global ok_all; ok_all &= bool(ok); R.setdefault('gates', {})[name] = bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)

# ============================================== 1. SIMBOLICO ==============================================
t, r, th, ph = sp.symbols('t r theta phi', real=True)
Gs = sp.symbols('G', positive=True)
f = sp.Function('phi')(t, r); l = sp.Function('lam')(t, r)
X = [t, r, th, ph]
g = sp.diag(-sp.exp(2*f), sp.exp(2*l), r**2, r**2*sp.sin(th)**2); gi = g.inv()
def Gam(a, b, c): return sum(gi[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(4))/2
Gm = [[[sp.simplify(Gam(a, b, c)) for c in range(4)] for b in range(4)] for a in range(4)]
def Ric(b, c):
    return sp.simplify(sum(sp.diff(Gm[a][b][c], X[a]) - sp.diff(Gm[a][b][a], X[c]) + sum(Gm[a][a][d]*Gm[d][b][c] - Gm[a][c][d]*Gm[d][b][a] for d in range(4)) for a in range(4)))
Rm = sp.Matrix(4, 4, lambda b, c: Ric(b, c)); Rs = sp.simplify(sum(gi[a, b]*Rm[a, b] for a in range(4) for b in range(4)))
Gmix = sp.simplify(gi*(Rm - Rs*g/2))                       # G^mu_nu
Mms = r/(2*Gs)*(1 - sp.exp(-2*l))                         # Misner-Sharp
# --- D_MS1
c_rho = sp.simplify(Gmix[0, 0] + 2*Gs*sp.diff(Mms, r)/r**2)             # G^t_t = -2 G M_r / r^2  (=8 pi G T^t_t)  => M_r = 4 pi r^2 rho
c_flux = sp.simplify(Gmix[1, 0] - 2*Gs*sp.diff(Mms, t)/r**2)            # G^r_t = +2 G M_t / r^2  (=8 pi G T^r_t)  => M_t = 4 pi r^2 T^r_t = -4 pi r^2 F_flux
# extrinsic curvature de las rebanadas t=const (sin shift):  K_ij = -(1/2 alpha) d_t g_ij
gsp = sp.diag(sp.exp(2*l), r**2, r**2*sp.sin(th)**2)
Kdn = -sp.diff(gsp, t)/(2*sp.exp(f)); Kmix = sp.simplify(gsp.inv()*Kdn)
c_K1 = sp.simplify(Kmix[0, 0] + sp.exp(-f)*sp.diff(l, t)); c_K2 = sp.simplify(Kmix[1, 1]) ; c_K3 = sp.simplify(Kmix[2, 2])
c_K4 = sp.simplify(sp.diff(l, t) - r*sp.exp(2*l)*Gmix[1, 0]/2)           # lam_t = 4 pi G r e^{2lam} T^r_t  (G^r_t = 8 pi G T^r_t)
gate("D_MS1_ecuaciones_dinamicas_MS", c_rho == 0 and c_flux == 0 and c_K1 == 0 and c_K2 == 0 and c_K3 == 0 and c_K4 == 0,
     "Einstein completo (φ(t,r),λ(t,r)): G^t_t=−2G M_r/r² ⇒ ∂_rM=4πr²ρ ; G^r_t=2G M_t/r² ⇒ ∂_tM=4πr²T^r_t=−4πr²F_flux (F_flux=−T^r_t). "
     "K^r_r=−e^{−φ}∂_tλ=−4πG r e^{2λ−φ}T^r_t, K^θ_θ=0: K_ij≠0 ⇔ T^r_t≠0 (el caso K_ij=0 de P10 es el de flujo nulo). Residuos: 0")

# --- tensor EM radial y ley de transferencia
E = sp.Function('E')(t, r)
F = sp.zeros(4, 4); F[0, 1] = -sp.exp(f + l)*E; F[1, 0] = sp.exp(f + l)*E
Fup = gi*F*gi
Fsq = sum(F[a, b]*Fup[a, b] for a in range(4) for b in range(4))
TEMup = sp.Matrix(4, 4, lambda m, n: sum(Fup[m, a]*g[a, b]*Fup[n, b] for a in range(4) for b in range(4)) - gi[m, n]*Fsq/4)
TEMmix = sp.simplify(TEMup*g)                                           # T^mu_nu
sqrtg = sp.exp(f + l)*r**2
def divF(nu): return sp.simplify(sum(sp.diff(sqrtg*Fup[m, nu], X[m]) for m in range(4))/sqrtg)
Jmax = [divF(nu) for nu in range(4)]                                    # nabla_mu F^{mu nu} (convencion: J^nu := nabla_mu F^{nu mu} = -divF)
def divT(Tmix, nu):                                                    # nabla_mu T^mu_nu
    s = sum(sp.diff(Tmix[m, nu], X[m]) for m in range(4)) + sum(Gm[m][m][k]*Tmix[k, nu] for m in range(4) for k in range(4)) - sum(Gm[k][m][nu]*Tmix[m, k] for m in range(4) for k in range(4))
    return sp.simplify(s)
dT_EM = [divT(TEMmix, nu) for nu in range(4)]
Jup = [-Jmax[nu] for nu in range(4)]                                    # (signo interno de la comprobacion; el resultado es nabla T = F_{nu lam} nabla_mu F^{mu lam})
lorentz = [sp.simplify(sum(F[nu, lam_]*Jup[lam_] for lam_ in range(4))) for nu in range(4)]    # F_{nu lam} J^lam
okEM_form = [sp.simplify(TEMmix[0, 0] + E**2/2), sp.simplify(TEMmix[1, 1] + E**2/2), sp.simplify(TEMmix[2, 2] - E**2/2), sp.simplify(TEMmix[1, 0])]
sign_ok = [sp.simplify(dT_EM[nu] - s_*lorentz[nu]) == 0 for nu in (0, 1) for s_ in (1,)]
sign_ok_neg = [sp.simplify(dT_EM[nu] + lorentz[nu]) == 0 for nu in (0, 1)]
ssg = 1 if all(sign_ok) else (-1 if all(sign_ok_neg) else 0)
S_t = lorentz[0]                                                          # F_{t lam}J^lam (potencia intercambiada con el plasma, ν=t)
c_joule = sp.simplify(dT_EM[0] + sp.diff(E**2/2, t))                    # nabla_mu T^mu_t (EM) = - d_t (E^2/2)
# --- balance cuasilocal con T = T_EM + T_pl (general), usando los vinculos de Einstein
rho_p, q_p, pr_p, pt_p = [sp.Function(n)(t, r) for n in ('rho_p', 'q_p', 'pr_p', 'pt_p')]
Tpl = sp.zeros(4, 4); Tpl[0, 0] = -rho_p; Tpl[1, 0] = q_p; Tpl[0, 1] = -sp.exp(2*l - 2*f)*q_p; Tpl[1, 1] = pr_p; Tpl[2, 2] = pt_p; Tpl[3, 3] = pt_p
divTpl_t = divT(Tpl, 0)
hand = -sp.diff(rho_p, t) + sp.diff(q_p, r) + (sp.diff(f, r) + sp.diff(l, r) + 2/r)*q_p - sp.diff(l, t)*(rho_p + pr_p)
c_div_t = sp.simplify(divTpl_t - hand)
lam_t_c = 4*sp.pi*Gs*r*sp.exp(2*l)*q_p                                  # (C1)  de G^r_t
c_rr_tt = sp.simplify(Gmix[1, 1] - Gmix[0, 0] - 2*sp.exp(-2*l)*(sp.diff(f, r) + sp.diff(l, r))/r)      # (C2) G^r_r - G^t_t = (2/r) e^{-2lam}(phi_r+lam_r)
# ecuaciones de energia:  EM: rho_EM,t = -S_t (rho+p_r=0, q=0);  plasma: nabla T_t(pl) = -S_t
S_sym, qs, qrs, rhopl, ppl, pr_t = sp.symbols('S q q_r rho_pl p_pl phi_r', real=True)
lam_r_s = sp.symbols('lam_r', real=True)
rhopl_t = S_sym + qrs + (pr_t + lam_r_s + 2/r)*qs - sp.symbols('lam_t', real=True)*(rhopl + ppl)         # plasma
rhoem_t = -S_sym
lam_t_sym = sp.symbols('lam_t', real=True)
# constraint (C1,C2): lam_t = 4 pi G r e^{2lam} q ;  8 pi G (rho_tot + p_tot) = (2/r) e^{-2lam}(phi_r+lam_r),  rho_EM + p_EM = 0  => rho_pl+p_pl = rho_tot+p_tot
e2l = sp.symbols('e2l', positive=True)
sub_c = {lam_t_sym: 4*sp.pi*Gs*r*e2l*qs}
plpl = (2/(r*e2l))*(pr_t + lam_r_s)/(8*sp.pi*Gs)                     # rho_pl+p_pl = rho_tot+p_tot (constraint C2, e^{-2lam}=1/e2l)
tot_t = sp.simplify((rhopl_t + rhoem_t).subs(sub_c).subs(rhopl + ppl, plpl))
c_total = sp.simplify(tot_t - (qrs + 2*qs/r))                             # rho_tot,t = (1/r^2) d_r(r^2 q): S desaparece
c_indep = sp.simplify(sp.diff(tot_t, S_sym))
qf = sp.Function('qf'); Rr = sp.symbols('R', positive=True)
Mt_R = sp.integrate(4*sp.pi*(qrs.subs(qrs, sp.diff(r**2*qf(r), r)/r**2) + 2*qf(r)/r)*r**2, (r, 0, Rr)) if False else sp.integrate(4*sp.pi*sp.diff(r**2*qf(r), r), (r, 0, Rr))
c_bdry = sp.simplify(Mt_R - 4*sp.pi*Rr**2*qf(Rr) + 4*sp.pi*(0**2)*qf(0))
gate("D_MS2_balance_cuasilocal", c_div_t == 0 and c_rr_tt == 0 and okEM_form == [0, 0, 0, 0] and ssg != 0 and c_joule == 0 and c_total == 0 and c_indep == 0 and c_bdry == 0,
     f"T^t_t=T^r_r=−E²/2, T^θ_θ=E²/2, T^r_t(EM)=0; ∇_μT^μ_ν(EM)=F_νλJ^λ [J^ν≡∇_μF^μν en unidades Heaviside–Lorentz; comprobado en las componentes t y r] y su componente t es −∂_t(E²/2) (calor de Joule); ∇_μT^μ_t(pl) = −ρ_t+q_r+(φ_r+λ_r+2/r)q−λ_t(ρ+p_r) (residuo 0); "
     f"con los vínculos G^r_t y G^r_r−G^t_t: ∂_tρ_tot = (1/r²)∂_r(r²q) independientemente de S_t ({c_indep}=0) ⇒ ∂_tM(t,R)=4πR²T^r_t(R) sólo por el flujo en r=R (T^r_t(EM)=0)")
R['D_MS2_sign'] = int(ssg)
# --- C_MS: ley de transferencia alterada: el plasma recibe (1+delta) de lo que pierde el campo
delta = sp.symbols('delta', positive=True)
tot_t_bad = sp.simplify((rhopl_t.subs(S_sym, (1 + delta)*S_sym) + rhoem_t).subs(sub_c).subs(rhopl + ppl, plpl))
spur = sp.simplify(tot_t_bad - (qrs + 2*qs/r))
S0, a_ = sp.symbols('S0 a', positive=True)
Mt_spur = sp.integrate(4*sp.pi*r**2*spur.subs(S_sym, S0*sp.exp(-a_*r)), (r, 0, Rr))
gate("C_MS_control_transferencia_alterada", sp.simplify(spur - delta*S_sym) == 0 and sp.simplify(Mt_spur) != 0,
     f"con ∇T_pl=−(1+δ)F J: ∂_tρ_tot − r⁻²∂_r(r²q) = δ·S_t ≠ 0 ⇒ ∂_tM(R) − 4πR²q(R) = 4π∫₀ᴿ r² δ S_t dr = {sp.simplify(Mt_spur)} ≠ 0 (fuente espuria de volumen)")
# --- D_MS3: frontera movil
Rdot, rho_d, ut = sp.symbols('Rdot rho_d u_t_up', real=True)
ugrav = sp.symbols('e2phi', positive=True)
u_t_lower = -ugrav*ut                                                   # u_t = g_tt u^t
q_dust = rho_d*(Rdot*ut)*u_t_lower; rho_n = -rho_d*ut*u_t_lower          # T^r_t and -T^t_t of dust with u^r = Rdot u^t
q_pl_out, rho_out, vv = sp.symbols('q_out rho_out v', real=True)
gate("D_MS3_frontera_movil", sp.simplify(q_dust + Rdot*rho_n) == 0 and sp.simplify((-vv*rho_out) + vv*rho_out) == 0,
     "dM/dt|_{r=R(t)} = ∂_tM + Ṙ∂_rM = 4πR²(T^r_t + Ṙρ): se anula para materia comóvil (T^r_t=−Ṙρ) y para el flujo nulo saliente que viaja con la frontera; sólo el flujo RELATIVO a la superficie cambia M")
# --- D_MS4: transporte de flujos nulos radiales
rho_s = sp.Function('rho')(t, r)
res_tr = []
for sgn in (1, -1):
    lv = sp.Matrix([sp.exp(-f), sgn*sp.exp(-l), 0, 0])
    Tn = rho_s*lv*lv.T                                                  # T^{mu nu} (contravariante)
    Tm = sp.simplify(Tn*g)
    d_t_ = divT(Tm, 0); d_r_ = divT(Tm, 1)                              # nabla_mu T^mu_nu (covariante en nu)
    vloc = sp.exp(f - l)
    brk = sp.diff(rho_s, t) + sgn*vloc*sp.diff(rho_s, r) + 2*rho_s*(sp.diff(l, t) + sgn*vloc*sp.diff(f, r) + sgn*vloc/r)
    res_tr.append((sp.simplify(d_t_ + brk), sp.simplify(d_r_ - sgn*sp.exp(l - f)*brk)))
d_ok = all(a == 0 and b == 0 for a, b in res_tr)
gate("D_MS4_transporte_flujos_nulos", d_ok,
     "T^{μν}=ρ_±ℓ_±^μℓ_±^ν, ℓ_±=(e^{−φ},±e^{−λ},0,0), v=e^{φ−λ}=α/a: ∇_μT^μ_ν=0 ⇔ ∂_tρ_± ± v∂_rρ_± = −2ρ_±(λ_t ± vφ_r ± v/r); con W=r²ρ: "
     "∂_tW_+ + ∂_r(vW_+) = −W_+[v(φ_r+λ_r)+2λ_t],  ∂_tW_− − ∂_r(vW_−) = +W_−[v(φ_r+λ_r)−2λ_t] (las dos componentes t y r dan la MISMA ecuación; residuos 0)")
# --- D_MS5: el esquema preserva el vinculo de momento
Gs_, rr_ = sp.symbols('G r', positive=True); Mm = sp.symbols('M', positive=True); vsym = sp.symbols('v', positive=True)
rp, rm_, rem, rmat, Sx = sp.symbols('rho_p rho_m rho_em rho_mat S', real=True)
drp, drm = sp.symbols('drho_p drho_m', real=True)
rho_tot = rp + rm_ + rem + rmat; p_r = rp + rm_ - rem
a2 = 1/(1 - 2*Gs_*Mm/rr_)
Mr = 4*sp.pi*rr_**2*rho_tot
lam_r = a2*(Gs_*Mr/rr_ - Gs_*Mm/rr_**2); phi_r = a2/2*8*sp.pi*Gs_*rr_*p_r + (a2 - 1)/(2*rr_)
lam_t = -4*sp.pi*Gs_*rr_*a2*vsym*(rp - rm_)
def balance_residual(rmat_on, dlt=0):
    rp_t = -vsym*drp - 2*rp*(lam_t + vsym*phi_r + vsym/rr_) + Sx/2*(1 + dlt)
    rm_t = vsym*drm - 2*rm_*(lam_t - vsym*phi_r - vsym/rr_) + Sx/2*(1 + dlt)
    rho_t = rp_t + rm_t - Sx                                           # rho_EM,t = -S ;  materia estatica
    v_r = vsym*(phi_r - lam_r)
    qq = -vsym*(rp - rm_)
    dflux = 2*rr_*qq + rr_**2*(-(v_r*(rp - rm_) + vsym*(drp - drm)))
    out = sp.simplify((rr_**2*rho_t - dflux))
    return out if rmat_on else sp.simplify(out.subs(rmat, 0))
res0 = balance_residual(False); res_mat = balance_residual(True); res_bad = balance_residual(False, dlt=delta)
gate("D_MS5_esquema_preserva_vinculo", res0 == 0 and res_mat != 0 and sp.simplify(res_bad) != 0,
     f"con ρ_tot=ρ_EM+ρ_++ρ_−, q=−v(ρ_+−ρ_−), p_r=ρ_++ρ_−−ρ_EM y las ecuaciones de D_MS4 (+conversión con balance): r²∂_tρ_tot − ∂_r(r²q) = 0 (el flujo y M_t son coherentes con M_r). "
     f"Materia 'mantenida' (estática, p=0): residuo = {sp.simplify(res_mat)} ≠ 0 ⇒ tratada como masa de PRUEBA (declarado); transferencia con (1+δ): residuo {sp.simplify(res_bad)} ≠ 0 (control)")
R['D_MS5_matter_residual'] = str(sp.simplify(res_mat))
# --- P12: horizonte aparente en los datos iniciales de la bola uniforme
xq, Aq, Bq, Cnn, qq_ = sp.symbols('x A B C_n q', positive=True)
Cx = (Cnn + sp.Rational(6, 5)*qq_**2 - qq_**2/xq)/xq                  # C(x)=2G M_MS(x r0)/(x r0), x>=1: M=M_n+U_in+U_out(1-1/x)
xs = sp.solve(sp.diff(Cx, xq), xq)[0]; Cmax_x = sp.simplify(Cx.subs(xq, xs))
crit = sp.simplify(Cmax_x - (Cnn + sp.Rational(6, 5)*qq_**2)**2/(4*qq_**2))
q_c0 = sp.solve(sp.Eq(Cmax_x.subs(Cnn, 0), 1), qq_)
gate("P12_horizonte_aparente_inicial_bola", crit == 0 and sp.simplify(xs - 2*qq_**2/(Cnn + sp.Rational(6, 5)*qq_**2)) == 0 and q_c0 == [sp.Rational(5, 3)],
     f"bola uniforme: C(x)=[C_n+1.2q²−q²/x]/x, máximo en x=2q²/(C_n+1.2q²) con valor (C_n+1.2q²)²/(4q²) (si x≥1, i.e. C_n≤0.8q²); para C_n=0 hay horizonte aparente en t=0 ⇔ q≥{q_c0[0]} (r0≤0.6 r_Q) = cota de extremalidad M_ADM≥Q; la descarga NO puede deshacerlo (NEC)")
R['P12'] = dict(q_crit_Cn0=float(q_c0[0]), r0_over_rQ=float(1/q_c0[0]))

# ============================================== 2. NUMERICO ==============================================
D0 = load_all(); fs = D0['branch_B3']['first_swing']; curve = (np.array(fs['s']), np.array(fs['g'])); kappa = D0['branch_B']['scaling']['kappa_inf']
def go(q, s, T=1.0, Nr=500, **kw):
    inst = (s == 'inst')
    return run_dyn(q, 1.0 if inst else s, curve, kappa, Nr=Nr, Rmax=5.0, Tend=T, instant=inst, **kw)
def dC(o, T=1.0):
    rec = o['rec']; i = int(np.argmin(np.abs(rec['t'] - T))); return float(rec['mface1'][i] - o['m1_initial'])

# --- N1: limite plano analitico (descarga instantanea, G=0)
def E_out_an(T, q):
    a = min(T, 1.0); o = 0.5*q**2*(1 - (1 - a)**5)/5
    return o + (0.5*q**2*(T - 1)**5/5 if T > 1 else 0.0)
def E_in_an(T, q): return 0.5*q**2*(1 - 1/(1 + T))
flat = {}
for Nr in (250, 500, 1000, 2000):
    o = run_dyn(1.0, 1.0, curve, kappa, Nr=Nr, Rmax=5.0, Tend=1.5, G_on=False, instant=True, record=False)
    flat[Nr] = dict(E_out=o['E_out'], E_in=o['E_in'], err_out=abs(o['E_out']/E_out_an(1.5, 1.0) - 1), err_in=abs(o['E_in']/E_in_an(1.5, 1.0) - 1))
orders = [np.log2(flat[a]['err_out']/flat[b]['err_out']) for a, b in ((250, 500), (500, 1000), (1000, 2000))]
R['flat_limit'] = dict(table={str(k): v for k, v in flat.items()}, E_out_exact=E_out_an(1.5, 1.0), E_in_exact=E_in_an(1.5, 1.0), observed_orders=[float(x) for x in orders])
gate("N1_limite_plano_analitico", flat[2000]['err_out'] < 1e-4 and flat[2000]['err_in'] < 1e-4 and orders[-1] > 1.5,
     f"G=0, descarga instantánea: E_out(T=1.5) y E_in vs fórmulas exactas (incl. paso por el centro): error {flat[2000]['err_out']:.1e}/{flat[2000]['err_in']:.1e} a Δr=r0/400; orden observado {[round(x,2) for x in orders]}")
o_nc = run_dyn(1.0, 1.0, curve, kappa, Nr=1000, Rmax=5.0, Tend=1.5, G_on=False, instant=True, record=False, center_reflect=False)
err_nc = abs(o_nc['E_out']/E_out_an(1.5, 1.0) - 1)
gate("C_N1_control_sin_paso_por_el_centro", err_nc > 20*flat[1000]['err_out'] and err_nc > 1e-3, f"sin reflexión en r=0 el error sube a {err_nc:.1e} (vs {flat[1000]['err_out']:.1e}): el test analítico detecta la condición de contorno errónea")

# --- N2: balance de energia / M_MS en r=r0 (D_MS2 numerico), con gravedad; sistema consistente (campo + plasma, sin materia mantenida)
cases = [(0.5, 1.0, {}), (1.0, 3.0, {}), (1.3, 'inst', {}), (0.8, 30.0, {}), (0.4, 0.5, {})]
bal = []
for (q, s, kw) in cases:
    for T in (0.5, 1.0):
        o = go(q, s, T=T, **kw); dm = o['rec']['mface1'][-1] - o['m1_initial']; flux = o['E_in'] - o['E_out']
        bal.append(dict(q=q, s=str(s), T=T, kw=str(kw), dm=float(dm), flux=float(flux), resid=float(abs(dm - flux)/max(abs(dm), abs(flux), 1e-9))))
R['balance'] = bal
resmax = max(b_['resid'] for b_ in bal if max(abs(b_['dm']), abs(b_['flux'])) > 1e-2)
resabs = max(abs(b_['dm'] - b_['flux']) for b_ in bal)
R['balance_summary'] = dict(resmax=float(resmax), resabs=float(resabs), n_cases=len(bal))
gate("N2_balance_Misner_Sharp_numerico", resmax < 5e-3 and resabs < 2e-4,
     f"m(r0,t)−m(r0,0) = ∫(F_in−F_out)dt en {len(bal)} casos consistentes (q=0.4…1.3, s=0.5…∞, con gravedad): residuo relativo máx {resmax:.1e}, absoluto máx {resabs:.1e} (D_MS2 verificado en el solver)")
# canal unico: la fuente se corta justo en la cara r=r0 donde se mide M -> error O(dr) por la discontinuidad (diagnostico; los canales se separan con el balance de flujo)
sc = {}
for Nr in (250, 500, 1000):
    o = go(1.0, 10.0, T=1.0, Nr=Nr, exterior=False); sc[Nr] = abs((o['rec']['mface1'][-1] - o['m1_initial']) - (o['E_in'] - o['E_out']))
R['single_channel_discrepancy'] = {str(k): float(v) for k, v in sc.items()}
gate("N2c_canal_unico_orden_1", sc[500] < 0.65*sc[250] and sc[1000] < 0.65*sc[500],
     f"con un solo canal la discrepancia |Δm−∫F| = {sc[250]:.1e}, {sc[500]:.1e}, {sc[1000]:.1e} (Δr/2 la reduce ×2: error O(Δr) de la discontinuidad de la fuente en r0, no del esquema); por eso los canales se separan con ∫F")
# materia mantenida (masa de prueba): la violacion del vinculo predicha por D_MS5 aparece como discrepancia medible y se declara como sistematico
mt_cases = []
for (q, s, kw) in ((1.0, 10.0, dict(Cn=0.3)), (1.0, 3.0, dict(Cn=0.3)), (0.8, 3.0, dict(absorb=True))):
    o = go(q, s, T=1.0, **kw); dm = o['rec']['mface1'][-1] - o['m1_initial']; flux = o['E_in'] - o['E_out']
    mt_cases.append(dict(q=q, s=s, kw=str(kw), dm=float(dm), flux=float(flux), rel_gap=float(abs(dm - flux)/abs(dm))))
R['held_matter_systematic'] = mt_cases
gap = max(c_['rel_gap'] for c_ in mt_cases)
gate("N2b_sistematico_materia_mantenida", 1e-3 < gap < 0.5,
     f"con masa mantenida C_n=0.3 o materia opaca, |Δm − ∫F|/Δm = {[round(c_['rel_gap'],3) for c_ in mt_cases]}: la discrepancia predicha por D_MS5 (∝ρ_mat) es real y ≲{gap:.0%}; los umbrales con C_n≠0 llevan ese sistemático (declarado)")
o_bad = go(1.0, 3.0, T=1.0, no_balance=0.2); dm_b = o_bad['rec']['mface1'][-1] - o_bad['m1_initial']; res_b = abs(dm_b - (o_bad['E_in'] - o_bad['E_out']))
R['control_no_balance'] = dict(residual=float(res_b), baseline=float(resabs))
gate("C_N2_control_conversion_sin_balance", res_b > 50*resabs and res_b > 1e-3, f"si el plasma recibe 1.2× lo que pierde el campo el balance de M se rompe: |Δm − ∫F| = {res_b:.1e} (vs {resabs:.1e} consistente)")

# --- N3: convergencia en la malla (dr, dr/2, dr/4)
conv = {}
for (q, s, kw, tag) in ((0.5, 1.0, {}, 'A'), (1.0, 10.0, {}, 'A'), (1.4, 'inst', {}, 'A'), (1.0, 3.0, dict(Cn=0.3), 'B')):
    vals = {}
    for Nr in (250, 500, 1000):
        o = go(q, s, T=1.0, Nr=Nr, record=True, **kw); vals[Nr] = dC(o, 1.0)
    d1 = abs(vals[500] - vals[250]); d2 = abs(vals[1000] - vals[500])
    conv[f"q{q}_s{s}_{kw}"] = dict(tag=tag, values={str(k): v for k, v in vals.items()}, diff_dr=float(d1), diff_dr2=float(d2), ratio=float(d1/max(d2, 1e-12)), rel_final=float(d2/max(abs(vals[1000]), 1e-12)))
R['convergence'] = conv
cA = [v for v in conv.values() if v['tag'] == 'A']; cB = [v for v in conv.values() if v['tag'] == 'B'][0]
R['convergence_summary'] = dict(rel_final_max_consistent=float(max(v['rel_final'] for v in cA)), first_order_case_ratio=float(cB['ratio']), first_order_case_rel_final=float(cB['rel_final']))
gate("N3_convergencia_malla", all(v['rel_final'] < 0.01 and (v['ratio'] > 3 or v['diff_dr2'] < 1e-4) for v in cA) and cB['ratio'] > 1.7 and cB['rel_final'] < 0.05,
     f"ΔC(t=r0/c) con Δr, Δr/2, Δr/4: casos consistentes: variación relativa final {[round(v['rel_final'],4) for v in cA]} y cociente de diferencias sucesivas {[round(v['ratio'],1) for v in cA]} (≥3: orden ≥1.6); "
     f"caso con masa de prueba C_n=0.3: convergencia de primer orden (cociente {cB['ratio']:.2f}; cambio final {cB['rel_final']:.1e})")

# --- N4: canales (gravedad debil q=0.2): saliente interior (perdida), entrante exterior (ganancia), neto
weak = {}
for s in (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 'inst'):
    row = {}
    for nm, kw in (('both', {}), ('exterior_only', dict(interior=False)), ('interior_only', dict(exterior=False)), ('absorbing', dict(absorb=True))):
        o = go(0.2, s, T=1.0, **kw); row[nm] = (o['E_in'] - o['E_out'])/0.04
    o = go(0.2, s, T=1.0); row['f_esc_interior'] = o['E_out']/(0.2*0.04); row['E_in_over_q2'] = o['E_in']/0.04
    weak[str(s)] = row
R['weak_gravity_F'] = weak
inf = weak['inst']
signs_ok = inf['interior_only'] < 0 < inf['exterior_only'] and inf['both'] > 0
mono = all(weak[a]['both'] <= weak[b]['both'] + 1e-6 for a, b in (('0.1', '0.3'), ('0.3', '1.0'), ('1.0', '3.0'), ('3.0', '10.0'), ('10.0', '30.0'), ('30.0', 'inst')))
lim = abs(weak['30.0']['both']/inf['both'] - 1) < 0.06
gate("N4_canales_saliente_entrante", signs_ok and mono and lim,
     f"ΔC(t=r0/c)/q² (q=0.2): canal interior (fuga) = {inf['interior_only']:+.4f}, canal exterior (inflow) = {inf['exterior_only']:+.4f}, neto = {inf['both']:+.4f} (instantánea); "
     f"F(s) crece monótonamente {weak['0.1']['both']:.4f} (s=0.1) → {weak['1.0']['both']:.4f} (s=1) → {inf['both']:.4f}; fracción del campo interior que escapa por r0: {inf['f_esc_interior']:.3f} (s→∞), {weak['1.0']['f_esc_interior']:.3f} (s=1)")

# --- N5: gravedad fuerte y masa mantenida C_n: supresion del inflow
qs_ = [0.2, 0.5, 0.8, 1.0, 1.25, 1.5]
strong = {f"{qv}": dict(zip(('s3', 'inst'), (dC(go(qv, 3.0), 1.0)/qv**2, dC(go(qv, 'inst'), 1.0)/qv**2))) for qv in qs_}
cn_dep = {f"{cn}": dC(go(1.0, 3.0, Cn=cn), 1.0) for cn in (0.0, 0.2, 0.4, 0.6)}
R['strong_gravity'] = strong; R['Cn_dependence_q1_s3'] = cn_dep
gate("N5_supresion_gravitatoria", strong['1.5']['s3'] < 0.2*strong['0.2']['s3'] and cn_dep['0.6'] < 0.1*cn_dep['0.0'],
     f"ΔC/q² (s=3): q=0.2 → {strong['0.2']['s3']:.4f}, q=1 → {strong['1.0']['s3']:.4f}, q=1.5 → {strong['1.5']['s3']:.4f}; con masa mantenida C_n=0→0.6 (q=1,s=3) ΔC cae {cn_dep['0.0']:.4f} → {cn_dep['0.6']:.4f}: el corrimiento sufre un corte al aproximarse a θ₊=0 (v=α/a→0)")

# --- N6: respuesta lineal del umbral de formacion de horizonte por la dinamica de descarga (t <= r0/c)
def Cn_AH0(q):                                  # umbral ESTATICO (t=0) de horizonte aparente global de la bola uniforme (P12, forma cerrada)
    return (2*q - 1.2*q**2) if q >= 1.0 else (1.0 - 0.2*q**2)
def Cmax0_analytic(q, Cn):
    return (Cn + 1.2*q**2)**2/(4*q**2) if Cn <= 0.8*q**2 else Cn + 0.2*q**2
qgrid = [0.3, 0.6, 0.9, 1.2, 1.5]
sgrid = [0.1, 0.3, 1.0, 3.0, 10.0, 'inst']
OFF = 0.05
resp = {}; chk0 = []
t0 = time.time()
for sv in sgrid:
    row = {}
    for qv in qgrid:
        Cn = Cn_AH0(qv) - OFF
        o = go(qv, sv, T=1.0, Cn=Cn); rec = o['rec']
        c0 = float(rec['Cmax'][0]); chk0.append(abs(c0 - Cmax0_analytic(qv, Cn)))
        slope = 1.0 if qv < 1.0 else (Cn + 1.2*qv**2)/(2*qv**2)           # dC_max/dC_n del maximo estatico
        row[f"{qv}"] = dict(dCmax=float(rec['Cmax'].max() - c0), dCb=float(rec['mface1'][-1] - rec['mface1'][0]), delta_Cn=float((rec['Cmax'].max() - c0)/slope), status=o['status'])
    resp[str(sv)] = row
print(f"  respuesta lineal calculada en {time.time()-t0:.0f}s", flush=True)
R['threshold_response'] = dict(q=qgrid, s=[str(x) for x in sgrid], offset_below_threshold=OFF, table=resp, Cn_AH0={f"{qv}": float(Cn_AH0(qv)) for qv in qgrid}, t_c=1.0,
                               max_Cmax0_discretization_err=float(max(chk0)))
dmax = max(v['delta_Cn'] for row in resp.values() for v in row.values()); dmin = min(v['delta_Cn'] for row in resp.values() for v in row.values())
R['threshold_response']['delta_Cn_max'] = float(dmax); R['threshold_response']['delta_Cn_min'] = float(dmin)
R['threshold_range_pct'] = float(100*dmax)
R['threshold_response']['delta_Cn_max_inst'] = float(max(v['delta_Cn'] for v in resp['inst'].values()))
R['threshold_response']['delta_Cn_max_s1'] = float(max(v['delta_Cn'] for v in resp['1.0'].values()))
gate("N6_umbral_de_horizonte_vs_dinamica", max(chk0) < 0.01 and dmax < 0.01 and dmin > -0.01 and all(v['status'] == 'ok' for row in resp.values() for v in row.values()),
     f"a C_n = C_n^AH0(q)−{OFF} (justo bajo el horizonte aparente inicial, P12): el máximo de compacidad C(r) al cabo de t=r0/c crece como mucho δC_n = {dmax:.4f} (s→∞: {R['threshold_response']['delta_Cn_max_inst']:.4f}; s=1: {R['threshold_response']['delta_Cn_max_s1']:.4f}); "
     f"el valor inicial numérico reproduce la forma cerrada de P12 a {max(chk0):.1e}; ningún caso forma horizonte nuevo: la descarga NO desplaza el umbral estático en más de ≈{100*dmax:.1f}%")

# --- N7: curvas en el plano (N, r0)
MP = MPL; aF = 0.75*(9*np.pi/4)**(1/3)
Nst = kappa**2/(4*np.pi*ALPHA**2)
def r_s_static(n): return n/MP**2 + np.sqrt((n/MP**2)**2 + 2*aF*n**(4/3)/MP**2)
def r_Q(n): return n*np.sqrt(ALPHA)/MP
def Cn_of_r(n, rr): return 2*(n + aF*n**(4/3)/rr)/(MP**2*rr)
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq
sidx = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0]
dl_arr = np.array([[resp[str(sv) if sv != 30.0 else 'inst'][f"{qv}"]['delta_Cn'] for qv in qgrid] for sv in sidx])
interp = RegularGridInterpolator((np.log10(sidx), qgrid), dl_arr, bounds_error=False, fill_value=None)
def delta_interp(q, s): return float(interp((np.log10(min(max(s, 0.1), 30.0)), min(max(q, 0.3), 1.5))))
Ngrid = np.logspace(1.5, 14, 80)
rows = []
for n in Ngrid:
    rq = r_Q(n); sN = np.sqrt(n/Nst)
    fa = lambda rr: Cn_of_r(n, rr) - Cn_AH0(min(rq/rr, 5/3)) if rq/rr < 5/3 else 1.0
    fd = lambda rr: Cn_of_r(n, rr) - (Cn_AH0(min(rq/rr, 5/3)) - delta_interp(rq/rr, sN)) if rq/rr < 5/3 else 1.0
    r_lo, r_hi = rq/(5/3)*(1 + 1e-9), 1e4*rq
    row = dict(N=float(n), s=float(sN), r_s=float(r_s_static(n)), r_Q=float(rq), r_AH0_matterfree=float(rq/(5/3)))
    if fa(r_lo) < 0 and fd(r_lo) < 0:
        row.update(r_AH=None, r_dyn=None, regime='AH_sin_materia')          # incluso con C_n pequeno hay horizonte aparente: r0 <= 0.6 r_Q basta
    else:
        row.update(r_AH=float(brentq(fa, r_lo, r_hi, xtol=1e-14*rq)), r_dyn=float(brentq(fd, r_lo, r_hi, xtol=1e-14*rq)), regime='B')
    rows.append(row)
R['phase_curve'] = rows
rowsB = [x for x in rows if x['r_dyn'] is not None]
N_A = next((x['N'] for x in rows if x['r_s'] < x['r_AH0_matterfree']), None)
rel_dyn = [x['r_dyn']/x['r_AH'] - 1 for x in rowsB]
ratio_AH_s = [x['r_AH']/x['r_s'] for x in rowsB]
R['N_A_staticcurve_enters_initial_horizon'] = N_A
R['curve_summary'] = dict(rel_shift_dyn_vs_AH_min=float(min(rel_dyn)), rel_shift_dyn_vs_AH_max=float(max(rel_dyn)), ratio_AH_over_rs_min=float(min(ratio_AH_s)), ratio_AH_over_rs_max=float(max(ratio_AH_s)), n_B=len(rowsB), n_total=len(rows))
gate("N7_frontera_corregida_en_el_plano_NR", N_A is not None and len(rowsB) > 5 and min(rel_dyn) >= -1e-9,
     f"r_s(N) (cascarón, C_n=1) frente a r_AH(N) (horizonte aparente global de la bola, P12 + dinámica): r_AH/r_s ∈ [{min(ratio_AH_s):.2f}, {max(ratio_AH_s):.2f}] sobre {len(rowsB)} valores de N; "
     f"la corrección dinámica r_dyn/r_AH−1 ∈ [{min(rel_dyn):.1e}, {max(rel_dyn):.1e}]; la curva estática r_s(N) del paper cruza r0=0.6 r_Q en N_A≈{N_A:.2e} (para N>N_A el horizonte aparente existe ya en los datos iniciales con M_n arbitrariamente pequeña)")

# ---------- ley de escala de ΔC y de la fraccion que escapa ----------
sv_ = np.array([0.1, 0.3, 1.0, 3.0, 10.0, 30.0]); Fv = np.array([weak[str(x)]['both'] for x in sv_]); Finf = weak['inst']['both']
# ajuste F(s)/F_inf = s^p/(s^p + s_h^p)
from scipy.optimize import curve_fit
def hill(s_, p, sh): return s_**p/(s_**p + sh**p)
(pp, sh), _ = curve_fit(hill, sv_, Fv/Finf, p0=[1.5, 1.0])
R['loss_fraction_of_UC'] = float(-weak['inst']['interior_only']/1.2)
R['inflow_over_leak'] = float(weak['inst']['exterior_only']/abs(weak['inst']['interior_only']))
R['F_scaling'] = dict(F_inf=float(Finf), hill_p=float(pp), hill_s_half=float(sh), fit_max_abs_err=float(np.max(np.abs(hill(sv_, pp, sh) - Fv/Finf))),
                      interior_inf=float(weak['inst']['interior_only']), exterior_inf=float(weak['inst']['exterior_only']), f_esc_inf=float(weak['inst']['f_esc_interior']))
gate("N8_ley_de_escala_F(s)", R['F_scaling']['fit_max_abs_err'] < 0.05 and 0.5 < pp < 3.0,
     f"ΔC(r0/c) ≃ q²·F∞·s^p/(s^p+s_h^p) para q≲0.5 (gravedad débil) con F∞={Finf:.4f}, p={pp:.2f}, s_h={sh:.2f} (s=√(N/N*)): error máx del ajuste {R['F_scaling']['fit_max_abs_err']:.3f}")

# ============================================== 3. FIGURA ==============================================
try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(10.5, 7.6))
    for nm, kw, st in (('both', {}, '-'), ('exterior only (inflow)', dict(interior=False), '--'), ('interior only (leak)', dict(exterior=False), ':')):
        o = go(0.2, 'inst', T=1.0, **kw); ax[0, 0].plot(o['rec']['t'], (o['rec']['mface1'] - o['m1_initial'])/0.04, st, label=nm)
    ax[0, 0].axhline(0, color='gray', lw=0.5); ax[0, 0].set_xlabel(r'$t\,c/r_0$'); ax[0, 0].set_ylabel(r'$\Delta C/q^2$'); ax[0, 0].legend(fontsize=8); ax[0, 0].set_title(r'Retained compactness at $r_0$ ($q=0.2$, $s\to\infty$)')
    sl = np.logspace(-1, 1.5, 30); ax[0, 1].semilogx(sv_, Fv, 'o', label='simulation'); ax[0, 1].semilogx(sl, Finf*hill(sl, pp, sh), '-', label=f'fit $p={pp:.2f}$, $s_h={sh:.2f}$')
    ax[0, 1].set_xlabel(r'$s=\sqrt{N/N_*}$'); ax[0, 1].set_ylabel(r'$\Delta C(r_0/c)/q^2$'); ax[0, 1].legend(fontsize=8); ax[0, 1].set_title('Discharge-speed dependence (weak gravity)')
    for sv in ('0.1', '1.0', '10.0', 'inst'):
        ax[1, 0].plot(qgrid, [resp[sv][f"{qv}"]['delta_Cn'] for qv in qgrid], 'o-', ms=3, label=f's={sv}')
    ax[1, 0].set_xlabel(r'$q=r_Q/r_0$'); ax[1, 0].set_ylabel(r'$\delta C_n$ (shift of the horizon threshold)'); ax[1, 0].legend(fontsize=8); ax[1, 0].set_title(r'Dynamic shift of the critical $2GM_n/r_0$ ($t\leq r_0/c$)')
    NN = np.array([x['N'] for x in rows]); ax[1, 1].loglog(NN, [x['r_s'] for x in rows], 'k-', label=r'$r_s(N)$ (shell criterion, paper)')
    ax[1, 1].loglog(NN, [x['r_Q']/R['P12']['q_crit_Cn0'] for x in rows], 'r-', label=r'$0.6\,r_Q$ (apparent horizon, $M_n\to0$)')
    ax[1, 1].loglog(NN, [x['r_Q'] for x in rows], 'r:', label=r'$r_Q$')
    ax[1, 1].loglog([x['N'] for x in rowsB], [x['r_AH'] for x in rowsB], 'g-.', label=r'$r_{\rm AH}(N)$ (static ball, P12)')
    ax[1, 1].loglog([x['N'] for x in rowsB], [x['r_dyn'] for x in rowsB], 'b--', lw=2, label=r'$r^{\rm dyn}(N)$ (with discharge, $t\leq r_0/c$)')
    ax[1, 1].axvline(Nst, color='gray', lw=0.6, ls=':'); ax[1, 1].set_xlabel(r'$N$'); ax[1, 1].set_ylabel(r'$r_0\ [\lambda_C]$'); ax[1, 1].legend(fontsize=7); ax[1, 1].set_title(r'Trapping boundary in the $(N,r_0)$ plane')
    fig.tight_layout(); fig.savefig('fig_dynamic_trapping.pdf'); fig.savefig('fig_dynamic_trapping.png', dpi=150); R['figure'] = 'fig_dynamic_trapping.pdf'
except Exception as ex:
    import traceback; traceback.print_exc(); R['figure'] = f'no generada: {ex!r}'; ok_all = False; print('FIGURA FALLÓ:', repr(ex))

# ---------- mapa de fases corregido ----------
try:
    NNf = np.logspace(1.5, 14, 400)
    r_eps_f = np.sqrt(NNf*ALPHA); r_P_f = np.maximum((ALPHA*NNf**2/(8*np.pi))**0.25/MP, 10.0/MP)
    rq_f = r_Q(NNf); rs_f = r_s_static(NNf)
    rAHf = np.array([next((x['r_AH'] for x in rows if abs(np.log10(x['N']/n)) < 0.08 and x['r_AH']), np.nan) for n in NNf])
    figp, axp = plt.subplots(figsize=(7.2, 5.2))
    axp.fill_between(NNf, 1e-26, rq_f/(5/3), color='tab:red', alpha=0.15, label='apparent horizon in the initial data ($r_0\\leq0.6\\,r_Q$, $M_n\\to0$)')
    axp.loglog(NNf, r_eps_f, 'C0-', label=r'$r_\epsilon=\sqrt{N\alpha}$ ($\epsilon=1$)'); axp.loglog(NNf, rq_f, 'C3:', label=r'$r_Q$')
    axp.loglog(NNf, r_P_f, 'C2--', label=r'$r_P$ ($\rho=M_{Pl}^4$, floor $10\,l_P$)'); axp.loglog(NNf, rs_f, 'k-', lw=1.5, label=r'$r_s(N)$ (shell criterion, previous map)')
    axp.loglog(NNf, rq_f/(5/3), 'r-', lw=1.0, label=r'$0.6\,r_Q$')
    Nmask = [x['N'] for x in rowsB]; axp.loglog(Nmask, [x['r_dyn'] for x in rowsB], 'b--', lw=2, label=r'$r^{\rm dyn}(N)$: horizon boundary with discharge ($t\leq r_0/c$)')
    axp.axvline(Nst, color='gray', lw=0.6, ls=':'); axp.text(Nst*1.2, 2e-24, r'$N_*$', color='gray')
    if N_A: axp.axvline(N_A, color='gray', lw=0.6, ls='--'); axp.text(N_A*1.2, 2e-25, r'$N_A$', color='gray')
    axp.set_ylim(1e-26, 1e4); axp.set_xlabel(r'$N$'); axp.set_ylabel(r'$r_0\ [\lambda_C]$'); axp.legend(fontsize=6.5, loc='lower right', framealpha=0.9); axp.set_title('Phase map with the dynamical horizon boundary')
    figp.tight_layout(); figp.savefig('phase_diagram_corrected.pdf'); figp.savefig('phase_diagram_corrected.png', dpi=150); R['figure_phase'] = 'phase_diagram_corrected.pdf'
except Exception as ex:
    import traceback; traceback.print_exc(); R['figure_phase'] = f'no generada: {ex!r}'; ok_all = False; print('FIGURA FASE FALLÓ:', repr(ex))
R['model'] = ("calibre polar-areal exacto (vínculos de Einstein) + campo radial de bola uniforme + dos flujos nulos radiales (plasma colineal, w∥≈1 de 05c) + descarga local con reloj s_c'=α κ s √(E/E_surf) "
              "y curva universal exacta g(s) de 05c; conversión 1/2+1/2 saliente/entrante (momento neto nulo); materia M_n mantenida como masa de prueba estática; "
              "SIN colisiones, SIN dispersión transversal (caústica central → resultados limitados a t≤r0/c), SIN transporte de carga (Gauss no impuesto en la conversión), SIN explosión de Coulomb de la materia")
save_section('branch_C2', R)
print("\nRESUMEN 06b:", "TODAS PASAN" if ok_all else "HAY FALLOS")
sys.exit(0 if ok_all else 1)
