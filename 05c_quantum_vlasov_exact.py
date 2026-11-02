"""CINETICA CUANTICA EXACTA (no markoviana): sistema de Bloch/Bogoliubov (f,u,v) acoplado a Maxwell. Sustituye el ansatz de Markov de 05b.
Para cada modo (P, s=p_perp^2): H_modo = p sigma_z + eps_perp sigma_x (reduccion exacta de Dirac), p=P+Phi, Phi'=x=eE, eps_perp^2=m^2+s, om=sqrt(p^2+eps_perp^2).
  f' = (lam/2) u ;  u' = lam (1-2f) - 2 om v ;  v' = 2 om u ;  lam = x eps_perp/om^2        (DERIVADO en SymPy desde n' = 2 h x n: gate Q0)
  x' = -2 e^2 sum_modos W [ (p/om) f + (eps_perp/(2 om)) u ]                                 (corriente de conduccion + de polarizacion coherente)
  rho_pl = 2 sum W om f ;  P_par = sum W p V,  V = 2 (p/om) f + (eps_perp/om) u   (velocidad <d_p H> sustraida del vacio)
SIN terminos fuente ad hoc, SIN LCFA, SIN ancho sigma: la coherencia de fase entre cruces y el bloqueo de Pauli (0<=f<=1) son exactos.
NOTA DE ALCANCE (declarada, ver tambien el paper): el enunciado del encargo escribe el acoplamiento como f'=(eE m/om^2)u, u'=(2eE m/om^2)(1-2f)-2 om v
(factor 2 de mas); eso conserva pureza y energia pero NO reproduce Schwinger-Nikishov: el gate C11 lo demuestra y Q0 deriva el factor correcto."""
import sys, os, time, ctypes, subprocess, warnings
warnings.filterwarnings("ignore")
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
from common import *
import markov_engine as ME

R = {}; ok_all = True
def gate(name, ok, msg):
    global ok_all; ok_all &= bool(ok); R.setdefault('gates', {})[name] = bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}")
e2 = E_CHARGE**2

# ================================ 1. SIMBOLICO ================================
t_ = sp.symbols('t', real=True)
x_, eps_, e_, c_, W_ = sp.symbols('x epsilon e c W', positive=True)
p_t = sp.Function('p', real=True)(t_)
om_t = sp.sqrt(p_t**2 + eps_**2)
Rx, Ry, Rz = [sp.Function(n, real=True)(t_) for n in ('Rx', 'Ry', 'Rz')]
h_vec = sp.Matrix([eps_, 0, p_t]); R_vec = sp.Matrix([Rx, Ry, Rz])
Rdot = 2*h_vec.cross(R_vec)                                  # n' = 2 h x n  (H = h.sigma);  R=-n  (vacio: z=+1)
hh = h_vec/om_t; e1 = sp.Matrix([p_t, 0, -eps_])/om_t; ey = sp.Matrix([0, 1, 0])
z_t = (R_vec.T*hh)[0]; u_t = (R_vec.T*e1)[0]; v_t = (R_vec.T*ey)[0]
def Dt(expr):
    d = sp.diff(expr, t_)
    d = d.subs({sp.Derivative(Rx, t_): Rdot[0], sp.Derivative(Ry, t_): Rdot[1], sp.Derivative(Rz, t_): Rdot[2]})
    return d.subs(sp.Derivative(p_t, t_), x_)
lam_t = x_*eps_/om_t**2
r1 = sp.simplify(Dt(z_t) + lam_t*u_t); r2 = sp.simplify(Dt(u_t) - (lam_t*z_t - 2*om_t*v_t)); r3 = sp.simplify(Dt(v_t) - 2*om_t*u_t)
f_t = (1 - z_t)/2
Vexp = sp.simplify(-Rz - (-p_t/om_t))                         # <sigma_z> = n_z = -R_z ; sustraccion del vacio (z=1,u=0 -> n_z=-p/om)
Vclaim = sp.simplify(2*p_t/om_t*f_t + eps_/om_t*u_t)
Eexp = sp.simplify(-(h_vec.T*R_vec)[0] - (-om_t))             # <H> sustraido = om(1-z) = 2 om f
Eclaim = sp.simplify(2*om_t*f_t)
gate("Q0_Dirac_a_Bloch", r1 == 0 and r2 == 0 and r3 == 0 and sp.simplify(Vexp - Vclaim) == 0 and sp.simplify(Eexp - Eclaim) == 0,
     "H=pσ_z+ε⊥σ_x, n'=2h×n; f=(1−z)/2, u=R·e1, v=R_y ⇒ f'=(λ/2)u, u'=λ(1−2f)−2ωv, v'=2ωu con λ=xε⊥/ω² (¡sin factor 2!); "
     "⟨∂_pH⟩_sust = 2(p/ω)f+(ε⊥/ω)u y ⟨H⟩_sust = 2ωf (residuos 0)")

f_, u_, v_, p_, W2 = sp.symbols('f u v p W2', real=True)
om_ = sp.sqrt(p_**2 + eps_**2); lam_ = c_*x_*eps_/om_**2                    # c=1 exacto; c=2 es el acoplamiento literal del enunciado
fdot = lam_*u_/2; udot = lam_*(1 - 2*f_) - 2*om_*v_; vdot = 2*om_*u_
def Dm(expr, fd=fdot, ud=udot, vd=vdot):
    return sp.diff(expr, f_)*fd + sp.diff(expr, u_)*ud + sp.diff(expr, v_)*vd + sp.diff(expr, p_)*x_
pur = (1 - 2*f_)**2 + u_**2 + v_**2
res_pur = sp.simplify(Dm(pur))
ineq = sp.solveset((1 - 2*f_)**2 <= 1, f_, sp.S.Reals)
ineq_ok = (ineq == sp.Interval(0, 1))
gate("Q1_pureza_Pauli_exacto", res_pur == 0 and ineq_ok,
     f"d/dt[(1−2f)²+u²+v²] = {res_pur} (idénticamente, para todo c) ; vacío (f,u,v)=(0,0,0) ⇒ pureza=1 ⇒ (1−2f)²≤1 ⇔ {ineq}: 0≤f≤1 sin truncar f=½")
res_bad = sp.simplify(Dm(pur, fd=lam_*u_))                                    # control: f' sin el 1/2
gate("C9_control_pureza_sin_1/2", res_bad != 0, f"con f'=λu (sin ½) la pureza se rompe: residuo = {sp.factor(res_bad)} ≠ 0")

# Energia: x'=-2e^2 W [ (p/om) f + c (eps/(2 om)) u ] ; rho_mode = 2 W om f ; E_campo=x^2/(2e^2)
xdot_full = -2*e_**2*W_*((p_/om_)*f_ + c_*(eps_/(2*om_))*u_)
def Dtot(xd):
    Efield = x_**2/(2*e_**2); rho = 2*W_*om_*f_
    return sp.diff(Efield, x_)*xd + sp.diff(rho, f_)*fdot + sp.diff(rho, p_)*x_
res_en = sp.simplify(Dtot(xdot_full))
gate("Q2_energia_exacta", res_en == 0, f"d/dt[E²/2 + 2∫ωf] = {res_en} (idénticamente, c genérico; con c=2 la corriente es 2e[(p/ω)f+(ε⊥/ω)u]: el sistema del enunciado también la conserva)")
xdot_nopol = -2*e_**2*W_*((p_/om_)*f_)
res_np = sp.simplify(Dtot(xdot_nopol))
c2 = sp.simplify(res_np.subs(c_, 2) - 2*W_*x_*eps_*u_/om_)
gate("C10_sin_corriente_de_polarizacion", res_np != 0 and c2 == 0,
     f"omitiendo (ε⊥/ω)u: residuo = {sp.simplify(res_np)} ≠ 0; para c=2 es 2W·eE·ε⊥u/ω (magnitud del enunciado, signo: dE_tot/dt sin esa corriente)")

# Traza: rho - P_par - 2 P_perp = sum W m^2 (2f - p u/eps)/om   (P_perp via p_perp d_pperp H, <sigma_x> sustraido = 2 f eps/om - u p/om)
Ppar_ = W_*(2*p_**2*f_/om_ + p_*eps_*u_/om_)
s_ = eps_**2 - 1                                                              # m=1
sx = 2*f_*eps_/om_ - u_*p_/om_
Pperp2 = W_*(s_/eps_)*sx
trace = sp.simplify(2*W_*om_*f_ - Ppar_ - Pperp2 - W_*(1/om_)*(2*f_ - p_*u_/eps_))
gate("Q4_traza_del_plasma", trace == 0, "ρ_pl − P_∥ − 2P_⊥ = Σ W m²(2f − p u/ε⊥)/ω (residuo 0): P_∥ y P_⊥ vienen de p·⟨∂_pH⟩ (Hellmann–Feynman); → 0 si m→0 (plasma conforme)")

xc = sp.solve(sp.exp(-sp.pi/sp.Symbol('X', positive=True)) - sp.Rational(1, 2), sp.Symbol('X', positive=True))
x_cross = float(xc[0])
gate("Q3a_umbral_f>1/2", abs(x_cross - float(sp.pi/sp.log(2))) < 1e-12 and 4.5 < x_cross < 4.6,
     f"exp(−π/ε)=½ ⇔ ε=π/ln2 = {x_cross:.5f} (ε>4.5: ocupación exacta > ½ para ε⊥=m)")
R['x_cross'] = x_cross
fz = sp.symbols('f0', positive=True)
defM = lambda f0: 1 - (1 - sp.exp(-2*f0))/(2*f0)
R['markov_single_crossing_deficit'] = dict(at_f0_0p9=float(defM(sp.Rational(9, 10))), limit_f0_1=float(defM(1)), x0_for_f0_0p9=float(sp.pi/(-sp.log(sp.Rational(9, 10)))))

from qvlasov_engine import *

# ================================ 3. LANDAU-ZENER EXACTO (Q3) ================================
def lz_scipy(x, eps2, c=1.0, Lfac=200.0, dressed=False):
    eperp = np.sqrt(eps2); L = Lfac*np.sqrt(x); T = L/x
    def rhs(t, y):
        f, u, v = y; p = x*t; om = np.sqrt(p*p + eps2); lam = c*x*eperp/om**2
        return [0.5*lam*u, lam*(1 - 2*f) - 2*om*v, 2*om*u]
    s = solve_ivp(rhs, [-T, T], [0, 0, 0], method='DOP853', rtol=1e-12, atol=1e-14, max_step=0.05/np.sqrt(x))
    f, u, v = s.y[:, -1]; return float(f), float((1 - 2*f)**2 + u*u + v*v)
lz = {}
for xv in (5.0, 10.0, 30.0, 100.0):
    for e2v in (1.0, 9.0):
        fi, pu = lz_scipy(xv, e2v); fe = float(np.exp(-np.pi*e2v/xv)); lz[f"x{xv:g}_eps2_{e2v:g}"] = dict(f_num=fi, f_exact=fe, rel_err=abs(fi/fe - 1), purity_dev=abs(pu - 1))
R['landau_zener'] = lz
R['n_lz_cases'] = len(lz); R['lz_fmax'] = max(v['f_num'] for v in lz.values()); R['lz_max_rel_err'] = max(v['rel_err'] for v in lz.values())
gt = [v for k, v in lz.items() if v['f_exact'] > 0.5]
gate("Q3_Landau-Zener_exacto", max(v['rel_err'] for v in lz.values()) < 1e-5 and len(gt) >= 4 and all(v['f_num'] > 0.5 for v in gt) and max(v['purity_dev'] for v in lz.values()) < 1e-8,
     f"integrador de (f,u,v) con E=const: máx |f/e^(−πε⊥²/eE) − 1| = {max(v['rel_err'] for v in lz.values()):.1e} en {len(lz)} casos; {len(gt)} con f_∞>½ (hasta {max(v['f_num'] for v in lz.values()):.4f}), donde Markov satura en ½; pureza {max(v['purity_dev'] for v in lz.values()):.1e}")
lz_bad = {k: lz_scipy(float(k.split('_')[0][1:]), float(k.split('_')[2]), c=2.0)[0] for k in ('x10_eps2_1', 'x30_eps2_1')}
bad_err = {k: abs(v/lz[k]['f_exact'] - 1) for k, v in lz_bad.items()}
R['landau_zener_literal_factor2'] = dict(f_num=lz_bad, rel_err=bad_err)
gate("C11_control_acoplamiento_literal_x2", min(bad_err.values()) > 0.3, f"con f'=(eEm/ω²)u, u'=(2eEm/ω²)(1−2f)−2ωv (literal del enunciado) f_∞ = {lz_bad} frente a {[lz[k]['f_exact'] for k in lz_bad]}: error relativo ≥ {min(bad_err.values()):.2f} (detectado)")

# nucleo C vs exacto (e2=0 => x constante, Phi = x t)
def lz_kernel(x, eps2s, Lfac=120.0):
    L = Lfac*np.sqrt(x); T = 2*L/x; om_max = L; dt = 0.1/om_max; ns = int(T/dt)
    P = np.full(len(eps2s), -L); s = np.array(eps2s) - 1.0; W = np.ones(len(eps2s)); st = np.array([x, 0.0]); f = np.zeros(len(P)); u = np.zeros(len(P)); v = np.zeros(len(P))
    o = kernel(P, s, W, st, f, u, v, dt, ns, ns, e2_=0.0)
    return f.copy(), o
fk, ok_ = lz_kernel(10.0, [1.0, 2.0, 4.0, 9.0])
fx = np.exp(-np.pi*np.array([1.0, 2.0, 4.0, 9.0])/10.0)
R['kernel_lz'] = dict(f_kernel=[float(a) for a in fk], f_exact=[float(a) for a in fx], rel_err=float(np.max(np.abs(fk/fx - 1))))
gate("Q3b_nucleo_C_vs_Schwinger", R['kernel_lz']['rel_err'] < 1e-5, f"núcleo C (RK4, E=const) vs e^(−πε⊥²/eE) en ε⊥²∈{{1,2,4,9}}, eE=10: error relativo máx {R['kernel_lz']['rel_err']:.1e}")

# ================================ 5. REGRESION DEL MOTOR DE MARKOV (M0) ================================
D0 = load_all(); mk_ref = D0['branch_B2']['sweep']['100']
mk100 = ME.analyze(ME.run(100.0, ncross=10, block=True))
gate("M0_regresion_Markov_05b", abs(mk100['tau0']/mk_ref['tau0'] - 1) < 1e-9 and abs(mk100['omega_late']/mk_ref['omega_late'] - 1) < 1e-9,
     f"markov_engine reproduce 05b a ε0=100: τ0 {mk100['tau0']:.8f} vs {mk_ref['tau0']:.8f}")

# ================================ 6. BARRIDO EN eps0: EXACTO vs MARKOV ================================
X0S = [4.0, 10.0, 30.0, 100.0, 300.0, 1000.0]
NCR = {4.0: 4}                                   # subumbral (ε0<π/ln2): la dinamica es lenta y la malla ancha; 4 cruces = 2 oscilaciones
exact = {}; markov = {}; series = {}
for x0 in X0S:
    r = ex_run(x0, ncross=NCR.get(x0, 10)); a, ser = ex_analyze(r); exact[f"{x0:g}"] = a; series[x0] = ser
    try:
        m = ME.analyze(ME.run(x0, ncross=NCR.get(x0, 10), block=True)); m['tau0_e_sqrtx0'] = m['tau0']*E_CHARGE*np.sqrt(x0)
        m['omega_over_esqrtx0'] = m['omega_late']/(E_CHARGE*np.sqrt(x0)) if m.get('omega_late') else None
        m['gamma_over_omega'] = m['gamma']/m['omega_late'] if m.get('gamma') and m.get('omega_late') else None
    except Exception as ex:
        m = dict(error=str(ex))
    markov[f"{x0:g}"] = m
    print(f"  ε0={x0:>6g}: exacto n_cross={a['n_cross']} τ0·e√ε0={a['tau0_e_sqrtx0']:.3f} ω_late/(e√ε0)={a['omega_over_esqrtx0']} γ/ω={a.get('gamma_over_omega')} fmax={a['fmax']:.5f} |ΔE|/E0={a['Emax_err']:.1e} ⟨w_∥⟩_pl={a.get('w_pl_avg_last_period')} ({a['runtime']:.0f}s) | Markov τ0·e√ε0={m.get('tau0_e_sqrtx0')} ω_late/(e√ε0)={m.get('omega_over_esqrtx0')} fmax={m.get('fmax')}", flush=True)
R['exact_sweep'] = exact; R['markov_sweep'] = markov

ncr_min = min(v['n_cross'] for k, v in exact.items() if float(k) >= 10); ncr_sub = min(v['n_cross'] for k, v in exact.items() if float(k) < 10); eerr = max(v['Emax_err'] for v in exact.values())
gate("Q2b_energia_numerica_multioscilacion", ncr_min >= 10 and ncr_sub >= 4 and eerr < 1e-6, f"≥{ncr_min} cruces por cero (≥5 oscilaciones) en ε0≥10 y ≥{ncr_sub} (2 oscilaciones) en el punto subumbral; máx |ΔE|/E0 = {eerr:.1e} < 1e-6 (el sistema discreto conserva la energía idénticamente; RK4 la preserva a ese nivel)")
R['max_energy_err'] = eerr
fmx = max(v['fmax'] for v in exact.values()); fmn = min(v['fmin'] for v in exact.values())
f_sup = {k: v['fmax'] for k, v in exact.items() if float(k) >= 10}
mk_f = max(v['fmax'] for v in markov.values() if 'fmax' in v)
gate("Q1c_0<=f<=1_cruza_1/2", fmn >= -1e-6 and fmx <= 1 + 1e-6 and all(v > 0.9 for v in f_sup.values()) and mk_f <= 0.5 + 1e-9,
     f"exacto: {fmn:.1e} ≤ f ≤ {fmx:.6f}; en ε0≥10 f_max = {min(f_sup.values()):.4f}…{max(f_sup.values()):.6f} (supera ½); Markov: f_max = {mk_f:.9f} ≤ ½ (truncado)")
# (la pureza numerica se evalua tras la convergencia en dt: gate Q1d mas abajo)
cvs = [v['cov'] for v in exact.values()]
gate("Q5_cobertura_de_malla", max(c['Phimax_over_Pa'] for c in cvs) < 0.92 and max(c['minusPhimin_over_Pb'] for c in cvs) < 0.92 and max(max(c['edge_f_left'], c['edge_f_right']) for c in cvs) < 1e-6,
     f"|Φ|max/P_a ≤ {max(c['Phimax_over_Pa'] for c in cvs):.2f}, −Φmin/P_b ≤ {max(c['minusPhimin_over_Pb'] for c in cvs):.2f}; f en los bordes ≤ {max(max(c['edge_f_left'], c['edge_f_right']) for c in cvs):.1e}")

# ================================ 7. CONVERGENCIA ================================
base = exact['100']; conv = {}
def metrics(a): return dict(tau0=a['tau0_e_sqrtx0'], hp=a['half_periods_u'][:5], omega_late=a['omega_over_esqrtx0'], pur_w=a['pur_w'])
for name, kw in (('dP/2', dict(dPf=0.02)), ('J=60', dict(J=60)), ('dt/2', dict(dtf=0.1)), ('Pa+30%', dict(Paf=52.0))):
    a, _ = ex_analyze(ex_run(100.0, ncross=6, **kw)); conv[name] = metrics(a)
dev_tau = max(abs(v['tau0']/base['tau0_e_sqrtx0'] - 1) for v in conv.values())
dev_hp = max(abs(np.array(v['hp'][:4])/np.array(base['half_periods_u'][:4]) - 1).max() for v in conv.values())
conv['max_rel_change_tau0'] = float(dev_tau); conv['max_rel_change_halfperiods_1to4'] = float(dev_hp)
# convergencia a eps0=10 (regimen no relativista de masa) con P y s mas finos
a10, _ = ex_analyze(ex_run(10.0, ncross=6, J=60)); b10 = exact['10']
conv['eps10_J=60'] = dict(tau0=a10['tau0_e_sqrtx0'], hp=a10['half_periods_u'][:5])
dev10 = max(abs(a10['tau0_e_sqrtx0']/b10['tau0_e_sqrtx0'] - 1), abs(np.array(a10['half_periods_u'][:4])/np.array(b10['half_periods_u'][:4]) - 1).max())
conv['eps10_max_rel_change'] = float(dev10); R['convergence'] = conv
gate("K6_convergencia", dev_tau < 5e-3 and dev_hp < 1.5e-2 and dev10 < 1.5e-2,
     f"ε0=100: variación máx. de τ0 = {dev_tau:.1e}, semiperiodos 1–4 = {dev_hp:.1e} bajo dP/2, J×2, Δt/2 y P_a+30%; ε0=10 (J×2): {dev10:.1e}")
pw, pw2 = base['pur_w'], conv['dt/2']['pur_w']
R['purity_numeric'] = dict(weighted_mean_dt=pw, weighted_mean_dt_half=pw2, ratio=pw2/pw, max_end_dt=base['pur_max_end'])
gate("Q1d_pureza_numerica", pw < 1e-2 and pw2 < pw/8, f"desviación media ponderada de (1−2f)²+u²+v²−1 (ε0=100): {pw:.1e} con Δt y {pw2:.1e} con Δt/2 (cociente {pw2/pw:.3f} ≈ 2⁻⁵: error de fase/amplitud de RK4 en modos de ω alto, que desaparece al refinar; la pureza exacta es el gate simbólico Q1)")
a_small, _ = ex_analyze(ex_run(100.0, ncross=3, Paf=20.0, autogrow=False))                                  # control: malla demasiado estrecha
dev_small = abs(a_small['tau0_e_sqrtx0']/base['tau0_e_sqrtx0'] - 1)
R['coverage_control'] = dict(tau0_narrow=a_small['tau0_e_sqrtx0'], rel_change=float(dev_small), cov=a_small['cov'])
gate("C12_control_malla_estrecha", a_small['cov']['Phimax_over_Pa'] > 0.92 or dev_small > 0.01, f"con P_a=20√ε0 la cobertura falla (|Φ|/P_a={a_small['cov']['Phimax_over_Pa']:.2f}) y τ0 cambia {dev_small:.1%}: el gate Q5 tiene capacidad de detección")

# ================================ 8. SENSIBILIDAD AL ESTADO INICIAL ================================
ad, _ = ex_analyze(ex_run(100.0, ncross=6, init='dressed'))
R['init_sensitivity'] = dict(tau0_vac=base['tau0_e_sqrtx0'], tau0_dressed=ad['tau0_e_sqrtx0'], rel_change_tau0=float(ad['tau0_e_sqrtx0']/base['tau0_e_sqrtx0'] - 1),
                             hp_vac=base['half_periods_u'][:4], hp_dressed=ad['half_periods_u'][:4], fmax_dressed=ad['fmax'], Emax_err_dressed=ad['Emax_err'])
dev_init = max(abs(R['init_sensitivity']['rel_change_tau0']), abs(np.array(ad['half_periods_u'][:4])/np.array(base['half_periods_u'][:4]) - 1).max())
R['init_sensitivity']['max_rel_change'] = float(dev_init)
gate("K8_sensibilidad_estado_inicial", dev_init < 0.05, f"vacío de orden cero (f=u=v=0) frente a vacío superadiabático de 1.er orden: variación máx. de τ0 y semiperiodos = {dev_init:.1%} (ambigüedad de los modos con |p|≲(eE ε⊥)^(1/3) en t=0; parámetro del ESTADO INICIAL, no de la dinámica)")

# ================================ 9. CONTROLES NEGATIVOS NUMERICOS ================================
c_np, _ = ex_analyze(ex_run(100.0, ncross=3, cpol=0.0, TmaxU=80.0))
gate("C13_sin_corriente_de_polarizacion_numerico", c_np['Emax_err'] > 1e-3 and c_np['Emax_err'] > 1e4*eerr, f"sin la corriente coherente: |ΔE|/E0 = {c_np['Emax_err']:.1e} ≫ {eerr:.1e}")
c_cl, _ = ex_analyze(ex_run(100.0, ncross=3, clam=2.0, cpol=2.0, TmaxU=80.0))
R['literal_factor2_plasma'] = dict(tau0_e_sqrtx0=c_cl['tau0_e_sqrtx0'], fmax=c_cl['fmax'], Emax_err=c_cl['Emax_err'])
gate("C14_acoplamiento_x2_en_el_plasma", abs(c_cl['tau0_e_sqrtx0']/base['tau0_e_sqrtx0'] - 1) > 0.05 and c_cl['Emax_err'] < 1e-6,
     f"con λ→2λ y corriente consistente (conserva energía: {c_cl['Emax_err']:.1e}) τ0·e√ε0 = {c_cl['tau0_e_sqrtx0']:.3f} frente a {base['tau0_e_sqrtx0']:.3f}: la conservación NO distingue el factor; solo Landau–Zener (Q3) lo fija")

# ================================ 10. COMPARACION MARKOV vs EXACTO ================================
cmp_ = {}
for k in exact:
    a, m = exact[k], markov[k]
    if 'tau0_e_sqrtx0' not in m: continue
    cmp_[k] = dict(tau0_ratio_exact_over_markov=a['tau0_e_sqrtx0']/m['tau0_e_sqrtx0'], omega_ratio=a['omega_over_esqrtx0']/m['omega_over_esqrtx0'],
                   gamma_over_omega_exact=a.get('gamma_over_omega'), gamma_over_omega_markov=m.get('gamma_over_omega'),
                   fmax_exact=a['fmax'], fmax_markov=m['fmax'], npairs_ratio_end=None)
R['comparison'] = cmp_
gate("D1_Markov_subestima_en_supercritico", all(cmp_[k]['fmax_exact'] > 0.9 and cmp_[k]['fmax_markov'] <= 0.5 + 1e-9 for k in cmp_ if float(k) >= 10) and abs(R['markov_single_crossing_deficit']['at_f0_0p9'] - 0.54) < 0.01,
     f"ε0≥10: f_max exacto > 0.9 frente a ≤ ½ de Markov; déficit de un cruce con f0=0.9 (ε≈{R['markov_single_crossing_deficit']['x0_for_f0_0p9']:.1f}) = {R['markov_single_crossing_deficit']['at_f0_0p9']:.1%} (el «54 %» del enunciado), límite f0→1: {R['markov_single_crossing_deficit']['limit_f0_1']:.1%}")

# ley de escala del plasma exacto (ε0 ≥ 30: regimen sin masa)
xs = np.array([v for v in X0S if v >= 30]); wl = np.array([exact[f"{v:g}"]['omega_late'] for v in xs])
expo = float(np.polyfit(np.log(xs), np.log(wl), 1)[0])
cs = [exact[f"{v:g}"]['omega_over_esqrtx0'] for v in xs]; gs = [exact[f"{v:g}"].get('gamma_over_omega') for v in xs]
R['scaling'] = dict(omega_exponent=expo, x0s=[float(v) for v in xs], omega_over_esqrtx0=[float(c) for c in cs], gamma_over_omega=[float(g) for g in gs],
                    omega_over_esqrtx0_markov=[markov[f"{v:g}"]['omega_over_esqrtx0'] for v in xs], tau0_e_sqrtx0=[exact[f"{v:g}"]['tau0_e_sqrtx0'] for v in xs])
gate("K7_escala_omega_p_exacta", abs(expo - 0.5) < 0.03 and (max(cs) - min(cs))/np.mean(cs) < 0.05,
     f"ω_p ∝ ε0^{expo:.3f}; ω_p/(e√ε0) = {[round(c,4) for c in cs]} (ε0∈{[int(v) for v in xs]}); γ/ω_p = {[round(g,4) for g in gs]}")

# ecuacion de estado
wavg = {k: dict(wpl=v.get('w_pl_avg_last_period'), wtot=v.get('w_tot_avg_last_period'), t_cross0=v['t_wtot_cross_0_u'], t_third=v['t_wpl_cross_third_u'], wpl_max=v['wpl_max']) for k, v in exact.items()}
R['eos'] = wavg
gate("Q6_ecuacion_de_estado", all(abs(v['wtot_t0'] + 1) < 1e-12 for v in exact.values()) and all(v['w_pl_avg_last_period'] is not None for v in exact.values()),
     f"w_tot=(P_∥−E²/2)/(ρ+E²/2) = −1 exacto en t=0 (campo puro, vacío-like); w_∥^pl=P_∥/ρ_pl promediado en el último periodo: {[round(v['wpl'],3) for v in wavg.values()]} (ε0={list(wavg.keys())})")

# ================================ 11. FIGURA ================================
try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(10.5, 7.2))
    for x0, c in ((10.0, 'C0'), (100.0, 'C1')):
        ser = series[x0]; ax[0, 0].plot(ser['t'], ser['x'], color=c, lw=1.1, label=rf'exact, $\epsilon_0={x0:g}$')
        mk = ME.run(x0, ncross=10, block=True); ax[0, 0].plot(mk['T']*E_CHARGE*np.sqrt(x0), mk['X']/x0, color=c, ls='--', lw=1.0, label=rf'Markov, $\epsilon_0={x0:g}$')
    ax[0, 0].set_xlim(0, 120); ax[0, 0].set_xlabel(r'$t\,e\sqrt{\epsilon_0}$'); ax[0, 0].set_ylabel(r'$\epsilon(t)/\epsilon_0$'); ax[0, 0].axhline(0, color='gray', lw=0.5); ax[0, 0].legend(fontsize=7); ax[0, 0].set_title('Field oscillations: exact vs Markov')
    xx = np.array(X0S); fm = [exact[f"{v:g}"]['fmax'] for v in xx]; fk = [markov[f"{v:g}"]['fmax'] for v in xx]
    ax[0, 1].semilogx(xx, fm, 'o-', label='exact (Bloch)'); ax[0, 1].semilogx(xx, fk, 's--', label='Markov (Pauli-blocked)'); ax[0, 1].axhline(0.5, color='gray', lw=0.6, ls=':')
    ax[0, 1].axvline(x_cross, color='gray', lw=0.6, ls=':'); ax[0, 1].set_xlabel(r'$\epsilon_0$'); ax[0, 1].set_ylabel(r'$\max f$'); ax[0, 1].legend(fontsize=8); ax[0, 1].set_title(r'Maximum occupation ($\epsilon=\pi/\ln 2$ dotted)')
    ax[1, 0].semilogx(xx, [exact[f"{v:g}"]['omega_over_esqrtx0'] for v in xx], 'o-', label=r'exact $\omega_p/(e\sqrt{\epsilon_0})$')
    ax[1, 0].semilogx(xx, [markov[f"{v:g}"].get('omega_over_esqrtx0') for v in xx], 's--', label='Markov')
    ax[1, 0].set_xlabel(r'$\epsilon_0$'); ax[1, 0].set_ylabel(r'$\omega_p/(e\sqrt{\epsilon_0})$'); ax[1, 0].legend(fontsize=8); ax[1, 0].set_title('Plasma frequency (mean of last 3 half-periods)')
    for x0, c in ((10.0, 'C0'), (100.0, 'C1')):
        ser = series[x0]; ax[1, 1].plot(ser['t'], ser['wtot'], color=c, lw=1.0, label=rf'$w_{{tot}}$, $\epsilon_0={x0:g}$')
        ax[1, 1].plot(ser['t'], ser['wpl'], color=c, ls=':', lw=1.0, label=rf'$w_\parallel^{{pl}}$, $\epsilon_0={x0:g}$')
    ax[1, 1].set_xlim(0, 120); ax[1, 1].set_ylim(-1.1, 1.3); ax[1, 1].axhline(1/3, color='gray', lw=0.5); ax[1, 1].set_xlabel(r'$t\,e\sqrt{\epsilon_0}$'); ax[1, 1].legend(fontsize=7); ax[1, 1].set_title('Equation of state (parallel)')
    fig.tight_layout(); fig.savefig('fig_markov_vs_exact.pdf'); fig.savefig('fig_markov_vs_exact.png', dpi=150); R['figure'] = 'fig_markov_vs_exact.pdf'
except Exception as ex:
    R['figure'] = f'no generada: {ex}'; ok_all = False; print('FIGURA FALLÓ:', ex)

# curva universal de la primera oscilacion (campo normalizado vs t e√ε0 a ε0=10³; invariante de escala para ε0≳30): cierre local de 06b
_s = series[1000.0]['t']; _g = series[1000.0]['x']; _i0 = int(np.argmax(_g <= 0)); _s0 = float(_s[_i0-1] + (_s[_i0]-_s[_i0-1])*_g[_i0-1]/(_g[_i0-1]-_g[_i0]))
_sg = np.linspace(0.0, _s0, 241); _gg = np.interp(_sg, _s[:_i0+1], np.clip(_g[:_i0+1], 0, None)); _gg[-1] = 0.0
R['first_swing'] = dict(s=[float(v) for v in _sg], g=[float(v) for v in _gg], s_zero=_s0, eps0=1000.0)
gate("Q7_curva_universal_primera_oscilacion", abs(_gg[0] - 1) < 1e-12 and _gg[-1] == 0.0 and bool(np.all(np.diff(_gg) <= 1e-9)) and abs(_s0 - exact['1000']['tau0_e_sqrtx0']) < 1e-6,
     f"campo normalizado g(s)=ε(t)/ε0 con s=t·e√ε0 (ε0=10³): monótona 1→0 en s0={_s0:.4f}={exact['1000']['tau0_e_sqrtx0']:.4f} (primer cero); se usa como cierre local de la descarga en 06b")
R['pct_tau0_shift_1000'] = (1 - cmp_['1000']['tau0_ratio_exact_over_markov'])*100
R['pct_omega_shift_1000'] = (cmp_['1000']['omega_ratio'] - 1)*100
R['min_fmax_ge10'] = min(f_sup.values())
R['max_wpl'] = max(v['wpl_max'] for v in exact.values())
R['tau0_vs_kappa_pct'] = abs(exact['1000']['tau0_e_sqrtx0']/D0['branch_B']['scaling']['kappa_inf'] - 1)*100
R['controls_C13_energy'] = c_np['Emax_err']
R['model'] = ("Bloch/Bogoliubov exacto para Dirac en E(t) homogeneo, retroaccion de Maxwell, sin LCFA ni ancho de fuente; estado inicial = vacio adiabatico de orden cero en t=0; "
              "sin colisiones ni transporte espacial (el plasma no se isotropiza); cuadratura en (P, p_perp^2) de paso finito; spin y e+- incluidos")
R['grid'] = dict(dP_over_sqrt_eps0=0.04, J=30, Pa_over_sqrt_eps0=40.0, s_max_over_eps0=14/np.pi, dt_rule="0.2/omega_max", integrator="RK4 de paso fijo (OpenMP)")
save_section('branch_B3', R)
print("\nRESUMEN 05c:", "TODAS PASAN" if ok_all else "HAY FALLOS")
sys.exit(0 if ok_all else 1)
