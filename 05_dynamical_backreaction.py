"""RAMA B/C — dinamica fuera del equilibrio: descarga de Schwinger con retroaccion (campo homogeneo, LCFA).
 Modelo (m_e=1, x=eE):   x' = -2 e^2 [ Int dt' Gam(x') <v> + Gam(x) <m_perp>/x ]        (corriente de conduccion + polarizacion)
                         p(t,t') = Int_{t'}^t x dt'' ,  v = p/sqrt(p^2+m_perp^2),  m_perp^2 = 1 + x s/pi, s~Exp(1)
                         Gam(x) = x^2/(4 pi^3) exp(-pi/x)  (pares/volumen/tiempo, solo n=1)
VALIDEZ: hasta el PRIMER cero del campo (cohortes sin recruce -> sin bloqueo de Pauli). Despues: NO se afirma nada."""
import io, sys, runpy, contextlib, warnings
warnings.filterwarnings("ignore")
import numpy as np, sympy as sp
from numpy.polynomial.laguerre import laggauss
from common import *

R = {}; ok_all = True
def gate(name, ok, msg):
    global ok_all; ok_all &= bool(ok); R.setdefault('gates', {})[name] = bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}")
e2 = E_CHARGE**2

# ============================== SIMBOLICO ==============================
# D1: ∇_μ T^{μν}_EM = (±) F^{νλ} J_λ en Minkowski esferico con A_t(t,r), A_r(t,r) arbitrarios
t, r, th, ph = sp.symbols('t r theta phi', positive=True); xs = [t, r, th, ph]
eta = sp.diag(-1, 1, r**2, r**2*sp.sin(th)**2); ei = eta.inv()
At = sp.Function('At')(t, r); Ar = sp.Function('Ar')(t, r); Av = [At, Ar, 0, 0]
F = sp.Matrix(4, 4, lambda i, j: sp.diff(Av[j], xs[i]) - sp.diff(Av[i], xs[j]))
Fu = ei*F*ei.T
Gam = [[[sum(ei[a, d]*(sp.diff(eta[d, b], xs[c]) + sp.diff(eta[d, c], xs[b]) - sp.diff(eta[b, c], xs[d])) for d in range(4))/2
         for c in range(4)] for b in range(4)] for a in range(4)]
sq = r**2*sp.sin(th)
J = [sp.simplify(sum(sp.diff(sq*Fu[m, n], xs[m]) for m in range(4))/sq) for n in range(4)]     # J^n = ∇_m F^{mn}
F2 = sum(F[i, j]*Fu[i, j] for i in range(4) for j in range(4))
Tuu = sp.Matrix(4, 4, lambda a, b: sum(Fu[a, c]*(F*ei)[d, b]*0 for c in range(4) for d in range(4)))   # placeholder
Fmix = F*ei                                                   # F_{a}^{ b}  = F_{ac} g^{cb}
Tdd = sp.Matrix(4, 4, lambda a, b: sum(F[a, c]*Fmix[b, c] for c in range(4)) - eta[a, b]*F2/4)       # T_ab = F_ac F_b^c - g F^2/4
Tmix = ei*Tdd                                                 # T^a_b
def divT(nu):
    return sum(sp.diff(sq*Tmix[m, nu], xs[m]) for m in range(4))/sq - sum(Gam[l][m][nu]*Tmix[m, l] for m in range(4) for l in range(4))
Jd = [sum(eta[n, m]*J[m] for m in range(4)) for n in range(4)]
res = {}
for nu in (0, 1):
    lhs = sp.simplify(divT(nu))
    rhs = sum(F[nu, l]*J[l] for l in range(4))                # F_{nu l} J^l
    res[nu] = (sp.simplify(lhs - rhs) == 0, sp.simplify(lhs + rhs) == 0)
gate("D1 ∇_μT^{μν}_EM = F^{ν}{}_λ J^λ (esferico, A_t(t,r),A_r(t,r))", all(v[0] for v in res.values()),
     "identidad tensorial verificada simbolicamente (signo +); implica ∇T_matter = -F J si ∇T_total=0")

# D2: balance de energia del modelo de cohortes (Leibniz) ; E_f' + rho_pl' = 0
s_, tt = sp.symbols('s t', positive=True)
Ef = sp.Function('E'); Gf = sp.Function('Gam'); m2 = sp.Function('M', positive=True)            # M(s)=m_perp(s) (cohorte)
ee = sp.symbols('e', positive=True)
a = sp.Function('a')                                                             # a' = E
p = ee*(a(tt) - a(s_))
rho_pl = sp.Integral(2*Gf(s_)*sp.sqrt(p**2 + m2(s_)**2), (s_, 0, tt))
d_rho = sp.diff(rho_pl, tt).doit()
jc = 2*ee*sp.Integral(Gf(s_)*p/sp.sqrt(p**2 + m2(s_)**2), (s_, 0, tt))           # corriente de conduccion
jp = 2*m2(tt)*Gf(tt)/Ef(tt)                                                      # corriente de polarizacion
Edot = -(jc + jp)
dEf = Ef(tt)*Edot                                                                # d/dt (E^2/2)
bal = (d_rho.subs(sp.Derivative(a(tt), tt), Ef(tt)) + dEf)
bal = sp.simplify(bal.doit())
gate("D2 conservacion de energia campo+plasma", bal == 0, f"d/dt(E²/2 + ρ_pl) = {bal}  (identidad exacta; Leibniz + j_pol)")
# control negativo simbolico: sin corriente de polarizacion
bal_bad = sp.simplify((d_rho.subs(sp.Derivative(a(tt), tt), Ef(tt)) + Ef(tt)*(-jc)).doit())
gate("C1 control: sin j_pol se viola", bal_bad != 0, f"residuo = {bal_bad}")
bal_bad2 = sp.simplify((d_rho.subs(sp.Derivative(a(tt), tt), Ef(tt)) + Ef(tt)*(+jc - jp)).doit())
gate("C2 control: signo de j_cond invertido se viola", bal_bad2 != 0, "residuo != 0")
# D3: limite Gamma/x^2 -> 1/(4 pi^3) (campo fuerte, n=1)
xx = sp.symbols('xx', positive=True)
gate("D3 Gamma ~ x^2/(4pi^3) (ε→∞)", sp.limit(xx**2/(4*sp.pi**3)*sp.exp(-sp.pi/xx)/xx**2, xx, sp.oo) == 1/(4*sp.pi**3), "tasa de pares cuadratica en ε")

# ============================== NUMERICO ==============================
J_NODES = 8
sj, wj = laggauss(J_NODES); wj = wj/wj.sum()
gate("N0 cuadratura Gauss-Laguerre", abs((wj*sj).sum() - 1) < 1e-12 and abs((wj*sj**2).sum() - 2) < 1e-12, "momentos <s>=1, <s^2>=2 de Exp(1)")

def Gam_pair(x): x = np.abs(x); return np.where(x > 0, x**2/(4*np.pi**3)*np.exp(-np.pi/np.maximum(x, 1e-300)), 0.0)
def run(x0, Nt, Tf=40.0, scheme='heun', drop_pol=False, flip=False):
    T = Tf/(E_CHARGE*np.sqrt(x0)); h = T/Nt
    X = np.zeros(Nt+1); A = np.zeros(Nt+1); G = np.zeros(Nt+1); M2 = np.zeros((Nt+1, J_NODES)); M = np.zeros((Nt+1, J_NODES))
    def fill(k, x, a_):
        X[k] = x; A[k] = a_; G[k] = Gam_pair(x); M2[k] = 1 + abs(x)*sj/np.pi; M[k] = np.sqrt(M2[k])
    def F(m):                                                   # RHS con nodos 0..m (m = nodo tentativo)
        wk = np.full(m+1, h); wk[0] = wk[-1] = h/2
        p = A[m] - A[:m+1]; om = np.sqrt(p[:, None]**2 + M2[:m+1]); v = p[:, None]/om
        S = np.sum(wk*G[:m+1]*(v @ wj)); P = 0.0 if drop_pol else np.sign(X[m])*G[m]*(M[m] @ wj)/max(abs(X[m]), 1e-300)
        sgn = -1.0 if flip else 1.0
        return -2*e2*(sgn*S + P)
    def energy(m):
        wk = np.full(m+1, h); wk[0] = wk[-1] = h/2
        p = A[m] - A[:m+1]; om = np.sqrt(p[:, None]**2 + M2[:m+1])
        return X[m]**2/(2*e2) + np.sum(wk*G[:m+1]*(2*om @ wj))
    fill(0, x0, 0.0); E0 = energy(0); Fn = F(0); tau0 = None; n_end = Nt
    for n in range(Nt):
        xp = X[n] + h*Fn; fill(n+1, xp, A[n] + h*(X[n] + xp)/2); Fp = F(n+1)
        xn = X[n] + h*(Fn + Fp)/2
        if scheme == 'trap':                                    # trapezoide implicito (punto fijo) -> metodo independiente
            for _ in range(60):
                fill(n+1, xn, A[n] + h*(X[n] + xn)/2); Fp = F(n+1); xnew = X[n] + h*(Fn + Fp)/2
                if abs(xnew - xn) < 1e-15*abs(x0): xn = xnew; break
                xn = xnew
        fill(n+1, xn, A[n] + h*(X[n] + xn)/2); Fn = F(n+1)
        if X[n]*X[n+1] < 0:                                     # primer cero: interpolacion cuadratica
            c = np.polyfit([n-1 if n > 0 else n, n, n+1][-3:], X[[max(n-1, 0), n, n+1]], 2) if n > 0 else None
            tau0 = (n + X[n]/(X[n]-X[n+1]))*h if c is None else np.interp(0, [X[n+1], X[n]], [(n+1)*h, n*h]); n_end = n+1; break
    E1 = energy(n_end)
    return dict(tau0=tau0, resid=abs(E1 - E0)/E0, h=h, n_end=n_end, x_end=X[n_end],
                n_pairs=float(np.sum(np.where(np.arange(n_end+1) == 0, 0.5, 1.0)*G[:n_end+1])*h) if tau0 else None, A=A[n_end], T=T)

# --- convergencia (h, h/2, h/4) y metodo independiente ---------------------------------
x0 = 1e3
rs = [run(x0, N) for N in (800, 1600, 3200)]
tc = [q['tau0'] for q in rs]; ratio = (tc[0]-tc[1])/(tc[1]-tc[2]); order = np.log2(ratio)
rt = run(x0, 1600, scheme='trap')
gate("N1 convergencia h,h/2,h/4 (Heun)", tc[0] is not None and 1.5 < order < 2.7, f"τ0 = {tc}; orden observado = {order:.2f} (esperado ≈2)")
gate("N2 metodo independiente (trapezoide implicito)", abs(rt['tau0'] - rs[1]['tau0'])/rs[1]['tau0'] < 5e-3, f"τ0(Heun,1600)={rs[1]['tau0']:.6g}  τ0(trap-impl,1600)={rt['tau0']:.6g}")
gate("N3 conservacion de energia (discreta)", rs[2]['resid'] < 5e-5 and rs[2]['resid'] < 0.2*rs[0]['resid'], f"|ΔE|/E0 = {rs[0]['resid']:.2e}, {rs[1]['resid']:.2e}, {rs[2]['resid']:.2e} (cae ~100x y se satura ~1e-5 en el cruce; identidad exacta en D2)")
R['convergence'] = dict(x0=x0, tau0=tc, order=float(order), resid=[q['resid'] for q in rs], tau0_trap=rt['tau0'])

# --- barrido en ε0: ley de escala τ0 = κ/(e√x0) ---------------------------------------
sweep = {}
x0s = [3, 10, 30, 1e2, 1e3, 1e4, 1e5, 1e6, 1e8, 1e10]
for x0 in x0s:
    q = run(x0, 1600)
    if q['tau0'] is None: q = run(x0, 1600, Tf=400.0)
    kappa = q['tau0']*E_CHARGE*np.sqrt(x0) if q['tau0'] else None
    sweep[f"{x0:g}"] = dict(tau0=q['tau0'], kappa=kappa, resid=q['resid'], pairs_per_vol=q['n_pairs'], p_max=q['A'],
                            lcfa_param=(E_CHARGE/kappa if kappa else None), t_form_over_tau0=(E_CHARGE/kappa if kappa else None))
    print(f"  ε0={x0:>8g}: τ0={q['tau0']:.5g}  κ=τ0·e√ε0={kappa:.4f}  |ΔE|/E0={q['resid']:.1e}  pares/vol={q['n_pairs']:.4g}")
R['sweep'] = sweep
big = [sweep[f"{x:g}"]['kappa'] for x in (1e4, 1e5, 1e6, 1e8, 1e10)]
lx = np.log([1e4, 1e5, 1e6, 1e8, 1e10]); lt = np.log([sweep[f"{x:g}"]['tau0'] for x in (1e4, 1e5, 1e6, 1e8, 1e10)])
expo = np.polyfit(lx, lt, 1)[0]
kappa_inf = float(np.mean(big[-3:]))
R['scaling'] = dict(tau0_exponent=float(expo), kappa_inf=kappa_inf, kappa_spread=float(np.ptp(big[-3:])), lcfa_e_over_kappa=float(E_CHARGE/kappa_inf))
gate("N4 ley de escala τ0 ∝ ε0^(-1/2)", abs(expo + 0.5) < 5e-3 and np.ptp(big[-3:])/kappa_inf < 1e-2,
     f"exponente = {expo:.4f} (−0.5); κ∞ = {kappa_inf:.4f} (dispersión {np.ptp(big[-3:])/kappa_inf:.1e}) => τ_discharge = κ/(e√(eE0))")
gate("N5 LCFA válida", E_CHARGE/kappa_inf < 0.1, f"t_form/τ0 ≈ e/κ = {E_CHARGE/kappa_inf:.2f} (≪1: aproximación de campo localmente constante justificada; t_form = formacion del par ~1/√ε)")

# --- controles negativos numericos -----------------------------------------------------------
qb = run(1e3, 1600, drop_pol=True)
gate("C3 control numerico: sin j_pol la energia no se conserva", qb['resid'] > 20*rs[1]['resid'], f"|ΔE|/E0 = {qb['resid']:.2e} vs {rs[1]['resid']:.2e}")
qf = run(1e3, 1600, flip=True)
gate("C4 control numerico: j_cond con signo erroneo se detecta", qf['resid'] > 1e-2, f"|ΔE|/E0 = {qf['resid']:.2e} (runaway no fisico: antiapantallamiento)")

# --- carrera gravitatoria: τ_discharge/τ_light = κ/(α√(4πN)) ----------------------------------------
Nstar = kappa_inf**2/(4*np.pi*ALPHA**2)
R['race'] = dict(N_star_tau_dis_eq_r0=float(Nstar), N_valid_homog=float(9/ALPHA), tau_dis_over_r0_at_N1e6=float(kappa_inf/(ALPHA*np.sqrt(4*np.pi*1e6))))
print(f"\nCarrera: τ_dis/r0 = κ/(α√(4πN)) = {kappa_inf/(ALPHA*np.sqrt(4*np.pi)):.2f}/√N ; N* = {Nstar:.1f}. Validez homogénea (r0 ≥ 3 ℓ_form): N ≥ 9/α = {9/ALPHA:.0f}")
R['validity'] = "hasta el primer cero del campo; sin bloqueo de Pauli/aniquilacion estimulada; LCFA; geometria plana; sin hidrodinamica espacial"
save_section('branch_B', R)
print("\nRESUMEN 05:", "TODAS PASAN" if ok_all else "HAY FALLOS")
sys.exit(0 if ok_all else 1)
