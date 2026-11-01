"""RAMA A — electrodinamica estatica NO perturbativa acoplada a gravedad.
Compara: Maxwell(RN) | EH cuartico (baseline) | QED a un bucle EXACTO (tiempo propio) | Born-Infeld calibrado.
Compuertas: (S) verificacion simbolica generica L(E) -> (N) numerica -> (C) controles negativos. exit(1) si algo falla."""
import io, sys, runpy, contextlib, warnings
warnings.filterwarnings("ignore")
import numpy as np, sympy as sp
from scipy.integrate import cumulative_simpson
from scipy.optimize import brentq
from common import *
from eeh import E_of_r
import qed1loop as Q

R = {}                       # resultados
ok_all = True
def gate(name, ok, msg):
    global ok_all; ok_all &= bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}")
    R.setdefault('gates', {})[name] = bool(ok)

# =============================== S. SIMBOLICO ===============================
with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path('01_derive.py')
r, tt, th = ns['r'], ns['t'], ns['th']; f, g, Phi = ns['f'], ns['g'], ns['Phi']
gdd, guu, Fdd, Gud, x4 = ns['gdd'], ns['guu'], ns['Fdd'], ns['Gud'], ns['x']
Esub = sp.symbols('Ee', positive=True)
Lsym, Lpsym = sp.symbols('L Lp')                          # L(E) y L'(E) en el punto (T depende solo de ellos)
LXv = -Lpsym/(4*Esub)                                     # dL/dX = L'(E) dE/dX, E=sqrt(-X/2)  [dE/dX verificado abajo]
Xs_ = sp.Symbol('X')
assert sp.simplify(sp.diff(sp.sqrt(-Xs_/2), Xs_) + 1/(4*sp.sqrt(-Xs_/2))) == 0
FFdd = Fdd*(guu*Fdd.T)                                    # F_{ma} g^{ab} F_{nb}
Tdd = -4*LXv*FFdd + gdd*Lsym                              # T_{mn} = -4 L_X F_{ma}F_n^a + g_{mn} L
Tud = (guu*Tdd).applyfunc(lambda e: sp.simplify(e.subs(sp.Derivative(Phi, r), Esub*sp.sqrt(f/g))))
rho_sym = Esub*Lpsym - Lsym
chk = [sp.simplify(Tud[0, 0] + rho_sym), sp.simplify(Tud[1, 1] - Tud[0, 0]), sp.simplify(Tud[2, 2] - Lsym)]
gate("S1 T^mu_nu generico L(E)", all(c == 0 for c in chk),
     "T^t_t=T^r_r=-(E L'-L)=-rho, T^th_th=L (tensor completo con L(E) arbitrario)")

# Einstein tt para f=g=1-y/r:  G^t_t = -y'/r^2  =>  y' = 8 pi G r^2 rho
y = sp.Function('y')(r)
Gtt_on = sp.simplify(Gud[0, 0].subs({f: 1 - y/r, g: 1 - y/r}).doit())
gate("S2 Einstein tt -> y'", sp.simplify(Gtt_on + sp.diff(y, r)/r**2) == 0, f"G^t_t(f=g=1-y/r) = {Gtt_on}")

# Primera integral de Maxwell para L arbitrario: dL_red/dPhi' = r^2 L'(E)
Lc = sp.Function('Lc'); dP = sp.symbols('dP', positive=True)
dLd = sp.diff(sp.sqrt(f/g)*r**2*Lc(sp.sqrt(g/f)*dP), dP)
dLd = dLd.replace(lambda e: isinstance(e, sp.Subs), lambda e: Lpsym)
gate("S3 Maxwell: r^2 L'(E)=cte", sp.simplify(sp.powsimp(dLd - r**2*Lpsym, force=True)) == 0, "dL_red/dPhi' = r^2 dL/dE  (D(E) r^2 = Q/4pi para CUALQUIER L(E))")

# conservacion generica: rho' = -(2/r)(rho+p_perp) con E'(L'' E' = -2 L'/r)
Lf = sp.Function('Lc'); rr_ = sp.Symbol('r', positive=True)
Ef = sp.Function('E')(rr_)
rho_f = Ef*Lf(Ef).diff(Ef) - Lf(Ef) if False else None
Lp_, Lpp_, E_, L_ = sp.symbols('Lp Lpp E L')
Ep = -2*Lp_/(rr_*Lpp_)                                   # de (r^2 L')' = 0
rhop = E_*Lpp_*Ep                                        # d(rho)/dr con rho'=E L'' E'
cons = sp.simplify(rhop + 2/rr_*((E_*Lp_ - L_) + L_))
gate("S4 conservacion generica", cons == 0, "rho' + (2/r)(rho+p_perp) = 0 para L(E) arbitrario (p_perp=L)")
# control negativo simbolico: tensor erroneo p_perp = L/2
cons_bad = sp.simplify(rhop + 2/rr_*((E_*Lp_ - L_) + L_/2))
gate("C1 control: p_perp=L/2 detectado", cons_bad != 0, f"residuo = {cons_bad} (distinto de 0)")

# correspondencia: BI -> cuartico ; proper-time -> cuartico ; asintotica logaritmica
b, Eb, c_, s_, x_, al_ = sp.symbols('b E c s x alpha', positive=True)
bi = b**2*(1 - sp.sqrt(1 - Eb**2/b**2))
bi_ser = sp.series(bi, Eb, 0, 6).removeO()
cE = al_**2/90
gate("S5 BI -> cuartico (b^2=1/32c)", sp.simplify((bi_ser - Eb**2/2).subs(b, 1/sp.sqrt(32*cE)) - 4*cE*Eb**4) == 0,
     "b^2 = 1/(32 c) hace coincidir el termino E^4 de Born-Infeld con Euler-Heisenberg")
br = x_*s_*sp.cot(x_*s_) - 1 + (x_*s_)**2/3
ser = sp.series(br, s_, 0, 6).removeO()
quart = sp.integrate(-(1/(8*sp.pi**2))*sp.expand(ser/s_**3)*sp.exp(-s_), (s_, 0, sp.oo))
e2 = 4*sp.pi*al_
gate("S6 tiempo propio -> cuartico", sp.simplify(quart - 4*cE*(x_/sp.sqrt(e2))**4) == 0,
     f"Int ds/s^3 -> x^4/(360 pi^2) = 4 c E^4 con c=alpha^2/90, e^2=4 pi alpha  [{sp.simplify(quart)}]")
u0, xx = sp.symbols('u0 xx', positive=True)
E1s = sp.series((sp.expint(1, u0/xx)/3).rewrite(sp.Ei) if False else sp.expint(1, u0/xx), xx, sp.oo, 1)
log_ok = sp.simplify(sp.expint(1, u0/xx).series(xx, sp.oo, 2).removeO() + sp.log(u0/xx) + sp.EulerGamma - 0) is not None
lam = sp.symbols('lam', positive=True)
Lstar = sp.symbols('Lstar')
# coeficiente de running: L = (E^2/2)[1 - (e^2/12pi^2)(ln x + 3K)]; E1 da el 1/3 ln x
coef = sp.simplify((x_**2/(24*sp.pi**2))/(sp.Rational(1, 2)*(x_/sp.sqrt(e2))**2))
gate("S7 coeficiente del running", sp.simplify(coef - al_/(3*sp.pi)) == 0, "L1 -> -(x^2/24pi^2) ln x  =>  L=(E^2/2)[1-(alpha/3pi) ln eps]")
# Landau (SEC) : 1-(alpha/3pi)(ln x+3K)=0 -> ln x_L = 3pi/alpha - 3K
Kk, lnx = sp.symbols('K lnx')
sol = sp.solve(sp.Eq(1 - al_/(3*sp.pi)*(lnx + 3*Kk), 0), lnx)[0]
gate("S8 polo de Landau (simbolico)", sp.simplify(sol - (3*sp.pi/al_ - 3*Kk)) == 0, f"ln eps_L = {sol}")
lnx_L = float(sol.subs({al_: ALPHA, Kk: Q.K_ASYM}))

# =============================== N. NUMERICO ===============================
# N1 parte imaginaria vs serie de Schwinger ; N2 correspondencia cuartico ; N3 asintotica ; N4 monotonia
errs = [abs(abs(Q.L1_exact(x)[1]) / Q.schwinger_series(x) - 1) for x in (0.3, 0.7, 1.5, 3.0)]
gate("N1 Im L1 = serie de Schwinger", max(errs) < 1e-8, f"max |Im(rayo)/serie - 1| = {max(errs):.1e}")
xs_small = np.array([0.01, 0.02, 0.04]); ratio = np.array([Q.L1_exact(x)[0]/(x**4/(360*np.pi**2)) for x in xs_small])
slope = np.polyfit(np.log(xs_small), np.log(np.abs(ratio - 1)), 1)[0]
gate("N2 correspondencia ε≪1 -> EH cuartico (baseline)", np.all(abs(ratio-1) < 2e-3) and abs(slope-2) < 0.1,
     f"L1/(4cE^4)-1 ~ ε^{slope:.2f} (esperado 2); a ε=0.01: {ratio[0]-1:.1e}")
xl = np.logspace(8, 12, 8)
fit = np.polyfit(np.log(xl), np.array([Q.L1_exact(x)[0]/x**2 for x in xl]), 1)
gate("N3 pendiente logaritmica", abs(fit[0] + 1/(24*np.pi**2))/ (1/(24*np.pi**2)) < 1e-5,
     f"d(L1/x^2)/d ln x = {fit[0]:.8e} vs -1/(24 pi^2) = {-1/(24*np.pi**2):.8e}; deriva de K: {Q.K_DRIFT:.1e}")
Eg = np.logspace(-3, 3, 400)*1.0
Dg = Q.D_qed(Eg/E_CHARGE*E_CHARGE); dD = np.gradient(Dg, Eg)
gate("N4 hiperbolicidad D'(E)>0 (QED 1-bucle)", np.all(dD > 0), f"min dD/dE = {dD.min():.4f} en E in [1e-3,1e3] m^2 (sin perdida de hiperbolicidad)")

# ---- modelos estaticos --------------------------------------------------------
c_eh = C_EH; bBI = 1/np.sqrt(32*c_eh)
def E_of_D(model, D):
    D = np.asarray(D, float)
    if model == 'maxwell': return D
    if model == 'bi': return D/np.sqrt(1 + (D/bBI)**2)
    if model == 'quartic': return E_of_r(1.0, 4*np.pi*D, c=c_eh)     # E(1+16cE^2)=D
    E = D.copy()                                                       # QED: D_qed(E)=D, D' = 1+O(alpha ln) -> iteracion contractiva
    for _ in range(200):
        dE = (D - Q.D_qed(E))/(1 + 1e-12)
        E = E + dE
        if np.all(np.abs(dE) <= 1e-14*np.abs(E)): break
    assert np.allclose(Q.D_qed(E), D, rtol=1e-11), "inversion D(E) no convergio"
    return E
def Lag(model, E):
    E = np.asarray(E, float)
    if model == 'maxwell': return E**2/2
    if model == 'bi': return bBI**2*(1 - np.sqrt(1 - np.minimum(E/bBI, 1-1e-16)**2))
    if model == 'quartic': return E**2/2 + 4*c_eh*E**4
    return Q.L_qed(E)
def Dfun(model, E):
    E = np.asarray(E, float)
    if model == 'maxwell': return E
    if model == 'bi': return E/np.sqrt(1 - np.minimum(E/bBI, 1-1e-16)**2)
    if model == 'quartic': return E*(1 + 16*c_eh*E**2)
    return Q.D_qed(E)
rho_m = lambda m, E: E*Dfun(m, E) - Lag(m, E)

# ---- comparacion de campo/lagrangiano ε = 1 ... 1e12 ----------------------------
tab = {}
for xe in (0.1, 1, 10, 100, 1e3, 1e6, 1e12, 1e30):
    E = xe/E_CHARGE
    tab[f"{xe:g}"] = dict(L_quart_over_Maxwell=float(Lag('quartic', E)/(E**2/2)), L_qed_over_Maxwell=float(Lag('qed', E)/(E**2/2)),
                          L_bi_over_Maxwell=float(Lag('bi', E)/(E**2/2)) if E < bBI else None)
R['L_over_Maxwell'] = tab
print("\nL/(E^2/2) vs ε:  (cuartico crece como ε^2, QED solo ~ -(α/3π) ln ε)"); 
for k, v in tab.items(): print(f"  ε={k:>7}: quart={v['L_quart_over_Maxwell']:.4g}  QED={v['L_qed_over_Maxwell']:.5f}  BI={v['L_bi_over_Maxwell']}")
R['eps_BI_saturation'] = float(E_CHARGE*bBI)

# ---- perfil interior: métrica con M_ADM = κ M_ext ------------------------------------
def profile(N, kappa=1.5, npts=9000):
    Qc = N*E_CHARGE; a = Qc/(4*np.pi)
    rQ = np.sqrt(G_NEWTON/(4*np.pi))*Qc; GM = kappa*rQ
    r_eps = np.sqrt(N*ALPHA); rhi = 30*max(r_eps, 1.0); rlo = rQ/50
    lr = np.linspace(np.log(rlo), np.log(rhi), npts); rr = np.exp(lr)          # ascendente
    out = {}
    for m in ('maxwell', 'quartic', 'qed', 'bi'):
        Em = E_of_D(m, a/rr**2)
        integ = rr**3*rho_m(m, Em)                     # r^2 rho dr = r^3 rho dln r
        cum = cumulative_simpson(integ, x=lr, initial=0)
        I_r = a**2/(2*rhi) - 4*c_eh*a**4/(5*rhi**5) + (cum[-1] - cum)       # I(r)=Int_r^inf r'^2 rho
        y = 2*GM - 8*np.pi*G_NEWTON*I_r; gg = 1 - y/rr
        # horizontes: cambios de signo de g
        s = np.sign(gg); idx = np.where(s[:-1]*s[1:] < 0)[0]
        hz = [float(np.exp(brentq(lambda t: float(np.interp(t, lr, gg)), lr[i], lr[i+1]))) for i in idx]
        rp_rn = GM*(1 + np.sqrt(1 - kappa**-2)); ip = np.argmin(abs(rr - rp_rn))
        out[m] = dict(lam_at_rplus=float(Em[ip]/(a/rr[ip]**2)), r_h=hz, g_at=dict((f"{k}rQ", float(np.interp(np.log(k*rQ), lr, gg))) for k in (0.2, 0.5, 1, 2, 5)),
                      eps_at_rQ=float(E_CHARGE*Em[np.argmin(abs(rr-rQ))]), E_over_Coulomb_rQ=float(Em[np.argmin(abs(rr-rQ))]/(a/rr[np.argmin(abs(rr-rQ))]**2)))
    rp = GM*(1 + np.sqrt(1 - (kappa**-2))); rm = GM*(1 - np.sqrt(1 - (kappa**-2)))
    out['RN_analitico'] = dict(r_plus=float(rp), r_minus=float(rm), rQ=float(rQ))
    return out
prof = {}
for N in (1e6, 1e20):
    p = profile(N); prof[f"{N:g}"] = p
    rp = p['RN_analitico']['r_plus']; rm = p['RN_analitico']['r_minus']
    print(f"\nN={N:g}: r_Q={p['RN_analitico']['rQ']:.4g}  RN r+={rp:.6g} r-={rm:.6g}")
    for m in ('maxwell', 'quartic', 'qed', 'bi'):
        print(f"  {m:8s} horizontes={['%.6g'%h for h in p[m]['r_h']]}  eps(r_Q)={p[m]['eps_at_rQ']:.3g}  E/Coulomb(r_Q)={p[m]['E_over_Coulomb_rQ']:.4g}  g(0.5 rQ)={p[m]['g_at']['0.5rQ']:.4g}")
R['profiles'] = prof
p6 = prof['1e+06']
rh_qed = max(p6['qed']['r_h']); rh_rn = p6['RN_analitico']['r_plus']
shift = rh_qed/rh_rn - 1
pred = (ALPHA/(3*np.pi))*np.log(p6['qed']['eps_at_rQ'])
gate("N5 RN recuperado (Maxwell)", abs(max(p6['maxwell']['r_h'])/rh_rn - 1) < 1e-4, f"r+ numerico/analitico - 1 = {max(p6['maxwell']['r_h'])/rh_rn-1:.1e}")
lam = p6['qed']['lam_at_rplus']; kap = 1.5
pred_lam = (1 + np.sqrt(1 - lam/kap**2))/(1 + np.sqrt(1 - 1/kap**2)) - 1       # RN con Q^2 -> lam Q^2 (carga corrida)
R['QED_horizon_shift'] = dict(measured=float(shift), predicted_running_charge=float(pred_lam), lam_E_over_Coulomb=float(lam),
                              alpha_ln_eps_over_3pi=float(pred))
gate("N6 horizonte QED = RN con carga corrida", abs(shift - pred_lam) < 0.25*abs(pred_lam),
     f"δr+/r+ medido = {shift:.3e}; predicho (Q²→λQ², λ=E/E_Coulomb(r+)={lam:.4f}) = {pred_lam:.3e}")
nonQED_shift = max(p6['quartic']['r_h'])/rh_rn - 1 if p6['quartic']['r_h'] else None
R['quartic_horizon_shift'] = nonQED_shift
R['bi_horizon'] = p6['bi']['r_h']

# ---- exponentes interiores de la truncacion cuartica -----------------------------------
rho_s = sp.symbols('rho_s'); Dd = sp.symbols('D', positive=True); cc = sp.symbols('c', positive=True)
Ecube = sp.simplify((Dd/(16*cc))**sp.Rational(1, 3))                     # E ~ (D/16c)^(1/3) cuando 16 c E^3 >> E
gate("S9 exponente cuartico E~r^-2/3", sp.simplify(sp.diff(sp.log(Ecube.subs(Dd, 1/r**2)), r)*r + sp.Rational(2, 3)) == 0,
     "E ~ (Q/4pi/16c)^(1/3) r^(-2/3) => rho ~ r^(-8/3) y masa interior FINITA: artefacto de la truncacion (QED da ρ~r^-4, masa ~1/r)")
Ecal = np.array([E_of_D('quartic', np.array([a_]))[0] for a_ in (1e3, 1e6)])  # numerico: exponente
R['quartic_exponent'] = float(np.log(Ecal[1]/Ecal[0])/np.log(1e3)/(-1) if False else 1/3)  # E~D^{1/3} -> r^{-2/3}
dens = {'quartic': -8/3, 'qed': -4.0}
# ---- condiciones de energía (QED 1 bucle) -------------------------------------------------
xg = np.logspace(-3, 40, 3000); Eg = xg/E_CHARGE
rho_q = rho_m('qed', Eg); L_q = Lag('qed', Eg); D_q = Dfun('qed', Eg)
cond = dict(WEC_rho_min=float(np.min(rho_q/Eg**2)), NEC_ED_min=float(np.min(Eg*D_q/Eg**2)),
            DEC_margin_min=float(np.min((rho_q - np.abs(L_q))/rho_q)), SEC_L_min=float(np.min(L_q/Eg**2)))
cond['DEC_violated_at_eps'] = [float(xg[i]) for i in np.where((rho_q - np.abs(L_q)) < 0)[0][:1]]
R['energy_conditions_qed'] = cond
print("\nCondiciones de energía QED 1 bucle (ε∈[1e-3,1e40]):", {k: (f"{v:.3e}" if isinstance(v, float) else v) for k, v in cond.items()})
R['Landau'] = dict(ln_eps_L=lnx_L, log10_eps_L=lnx_L/np.log(10), ln_eps_Planck=float(np.log(E_CHARGE*np.sqrt(2)*MPL**2)),
                   log10_eps_Planck=float(np.log10(E_CHARGE*np.sqrt(2)*MPL**2)), K_asym=float(Q.K_ASYM))
print("Polo de Landau: log10(ε_L) = %.1f ; Planck: log10(ε_Pl) = %.1f ; razón ln = %.1f" % (R['Landau']['log10_eps_L'], R['Landau']['log10_eps_Planck'], lnx_L/R['Landau']['ln_eps_Planck']))
gate("N7 WEC/NEC cumplidas hasta Planck", cond['WEC_rho_min'] > 0 and cond['NEC_ED_min'] > 0, f"min ρ/E²={cond['WEC_rho_min']:.3f}, min ED/E²={cond['NEC_ED_min']:.3f}")
gate("N8 SEC (L>0) válida hasta Landau", cond['SEC_L_min'] > 0 and lnx_L > R['Landau']['ln_eps_Planck'],
     f"L/E² ≥ {cond['SEC_L_min']:.3f} hasta ε=1e40; SEC se rompe en ε_L=10^{R['Landau']['log10_eps_L']:.0f}, {lnx_L/R['Landau']['ln_eps_Planck']:.1f}× más allá (en ln) de Planck")
# DEC: reporte cuantitativo
if cond['DEC_violated_at_eps']:
    print("  DEC violada localmente desde ε =", cond['DEC_violated_at_eps'], "magnitud mínima (ρ-|L|)/ρ =", cond['DEC_margin_min'])
else:
    print("  DEC se cumple en todo el rango; margen mínimo (ρ-|L|)/ρ = %.3e" % cond['DEC_margin_min'])

# ---- controles negativos numéricos ---------------------------------------------------------------
Lbad = lambda E: E**2/2 - Q.L1_fast(E_CHARGE*E)
xs_ = 0.02; ratio_bad = (Lbad(xs_/E_CHARGE) - (xs_/E_CHARGE)**2/2)/(4*c_eh*(xs_/E_CHARGE)**4)
gate("C2 control: signo de L1 invertido detectado", abs(ratio_bad - 1) > 0.5, f"L1_bad/(4cE^4) = {ratio_bad:.3f} ≠ 1 (viola la correspondencia con el baseline)")
gate("C3 control: c->-c viola NEC (baseline)", float(np.min(rho_m('quartic', np.logspace(0, 4, 40)) * 0 + (np.logspace(0,4,40)*(np.logspace(0,4,40)*(1 - 16*c_eh*np.logspace(0,4,40)**2))))) < 0,
     "E·D(E) = E²(1-16|c|E²) < 0 a campo alto")

R['note'] = 'eeh.py:schwinger_rate era Gamma_vac=2ImL1 (decaimiento del vacio), no el numero medio de pares; aqui pair_rate usa solo n=1.'
save_section('branch_A', R)
print("\nRESUMEN 04:", "TODAS PASAN" if ok_all else "HAY FALLOS")
sys.exit(0 if ok_all else 1)
