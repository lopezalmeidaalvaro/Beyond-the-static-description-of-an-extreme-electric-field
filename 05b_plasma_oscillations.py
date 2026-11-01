"""EXTENSION DE 05: cinetica cuantica (limite de Markov) con BLOQUEO DE PAULI a traves de multiples cruces por cero de E(t).
  d_t f(P,p_perp,t) = (1-2f) S(P,p_perp,t),   S = |x| exp(-pi m_perp^2/|x|) G_sigma(p),  p = P + Phi(t),  Phi' = x = eE
  dx/dt = -2 e^2 [ Sum W f p/omega  +  Sum W omega (1-2f) S / x ]            (conduccion + polarizacion CON el factor (1-2f))
  rho_pl = 2 Sum W f omega ;  identidad exacta  d/dt( x^2/(2e^2) + rho_pl ) = 0.
MODELO (declarado): campo homogeneo, LCFA, fuente de Markov (tasa de Schwinger n=1 repartida en un ancho sigma en p). La ecuacion de Markov con
(1-2f) satura en f<1/2 desde el vacio: NO reproduce la ocupacion exacta f=exp(-pi m^2/eps)>1/2 de la QKE no markoviana (limitacion cuantificada abajo)."""
import sys, warnings, time
warnings.filterwarnings("ignore")
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from common import *

R = {}; ok_all = True
def gate(name, ok, msg):
    global ok_all; ok_all &= bool(ok); R.setdefault('gates', {})[name] = bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}")
e2 = E_CHARGE**2

# ================================ SIMBOLICO ================================
P_, t_ = sp.symbols('P t', real=True)
f_ = sp.Function('f')(P_, t_); S_ = sp.Function('S')(P_, t_); Phi_ = sp.Function('Phi')(t_); x_ = sp.Function('x')(t_)
m2_, e_ = sp.symbols('m2 e', positive=True)
p_ = P_ + Phi_; om_ = sp.sqrt(p_**2 + m2_)
rho_int = sp.Integral(2*f_*om_, (P_, -sp.oo, sp.oo))                           # rho_pl (un modo p_perp; el resto es un peso constante)
drho = sp.diff(rho_int, t_)                                                    # Leibniz (limites fijos)
drho_f = drho.function if hasattr(drho, 'function') else drho
drho_f = drho_f.subs(sp.Derivative(f_, t_), (1 - 2*f_)*S_).subs(sp.Derivative(Phi_, t_), x_)
# campo: x' = -2 e^2 [ Int f p/om + Int om (1-2f) S / x ]  =>  d/dt(x^2/2e^2) = x x'/e^2   (integrando)
field_int = -2*x_*f_*p_/om_ - 2*om_*(1 - 2*f_)*S_
tot = sp.simplify(drho_f + field_int)
gate("K1_energia_con_(1-2f)", tot == 0, f"d/dt(E²/2+ρ_pl) = ∫dP [2ω(1−2f)S + 2x f p/ω − 2x f p/ω − 2ω(1−2f)S] = {tot}  (identidad exacta con el factor estadístico)")
bad_pol = sp.simplify(drho_f + (-2*x_*f_*p_/om_ - 2*om_*S_))                      # j_pol SIN (1-2f) pero f' CON (1-2f)
gate("C5_control_sin_(1-2f)_en_jpol", bad_pol != 0, f"residuo = {sp.simplify(bad_pol)} ≠ 0: la identidad se rompe si j_pol omite (1−2f)")
Sg, ts = sp.symbols('Sg', nonnegative=True), sp.symbols('s', real=True)
fc = (1 - sp.exp(-2*Sg))/2
ode_res = sp.simplify(sp.diff(fc, Sg) - (1 - 2*fc))                              # f(Sigma): df/dSigma = 1-2f  con Sigma'=S
ok_b1 = ode_res == 0 and sp.simplify(fc.subs(Sg, 0)) == 0
ok_b2 = sp.solveset(fc < 0, Sg, domain=sp.Interval(0, sp.oo)) == sp.EmptySet and sp.solveset(fc >= sp.Rational(1, 2), Sg, domain=sp.Interval(0, sp.oo)) == sp.EmptySet
gate("K2_bloqueo_exacto", ok_b1 and ok_b2, "f=(1−e^{−2∫S})/2 resuelve f'=(1−2f)S, f(0)=0 y 0 ≤ f < 1/2 para S≥0 (el bloqueo de Pauli de Markov acota f por 1/2)")
f0s = sp.symbols('f0', positive=True)
R['single_crossing'] = dict(rule="f_after=(1-exp(-2 f0))/2, f0=exp(-pi m_perp^2/eps)", ratio_at_f0_0p5=float(((1 - np.exp(-1.0))/2)/0.5), ratio_at_f0_0p9=float(((1 - np.exp(-1.8))/2)/0.9))

# ================================ NUMERICO ================================
def make(x0, sigf=0.25, dpf=1/3, J=10, Pfac=52.0, block=True, polblock=True, flip=False, wrongsign=False):
    sc = np.sqrt(x0); sig = sigf*sc; dp = sig*dpf
    Pmax = Pfac*sc + 5*sig; P = np.arange(-Pmax, Pmax + dp/2, dp); NP = len(P)
    pm = 4*np.sqrt(x0/np.pi)
    xg, wg = np.polynomial.legendre.leggauss(J); pp = 0.5*pm*(xg + 1); wp = 0.5*pm*wg
    W = (dp/(2*np.pi))*(2*pp*wp/(2*np.pi))[None, :]*np.ones((NP, 1))              # medida (incluye espin 2)
    m2 = 1 + pp**2
    def rhs(t, y):
        x, Phi = y[0], y[1]; f = y[2:].reshape(NP, J)
        p = (P + Phi)[:, None]; om = np.sqrt(p*p + m2[None, :]); ax = abs(x)
        A = (np.exp(-np.pi*m2/ax)[None, :]*np.exp(-0.5*(p/sig)**2)/(sig*np.sqrt(2*np.pi))) if ax > 0 else np.zeros((NP, J))
        S = ax*A
        fac = (1 + 2*f) if wrongsign else ((1 - 2*f) if block else 1.0)
        fdot = fac*S
        polfac = (1 + 2*f) if wrongsign else ((1 - 2*f) if (block and polblock) else 1.0)
        cond = np.sum(W*f*p/om); pol = np.sign(x)*np.sum(W*om*polfac*A)
        return np.concatenate(([-2*e2*((-1 if flip else 1)*cond + pol), x], fdot.ravel()))
    def energy(y):
        x, Phi = y[0], y[1]; f = y[2:].reshape(NP, J); p = (P + Phi)[:, None]; om = np.sqrt(p*p + m2[None, :])
        return x**2/(2*e2) + 2*np.sum(W*f*om)
    def npairs(y): return float(np.sum(W*y[2:].reshape(NP, J)))
    return dict(rhs=rhs, energy=energy, npairs=npairs, NP=NP, J=J, P=P, sig=sig, sc=sc, x0=x0)

def run(x0, ncross=10, rtol=1e-9, tmax_factor=60, **kw):
    K = make(x0, **kw); NP, J = K['NP'], K['J']
    y = np.zeros(2 + NP*J); y[0] = x0; E0 = K['energy'](y)
    tg = 13.5/(E_CHARGE*np.sqrt(x0)); t = 0.0; cr = []; Ed = []; fmax = 0.0; fmin = 0.0; Pcov = 0.0
    atol = np.concatenate(([1e-9*x0, 1e-9*K['sc']*30], np.full(NP*J, 1e-12)))
    T = []; X = []; Np = []; tstart = time.time()
    while len(cr) < ncross and t < tmax_factor*tg:
        sol = solve_ivp(K['rhs'], [t, t + tg], y, method='DOP853', rtol=rtol, atol=atol, dense_output=True)
        te = np.linspace(t, t + tg, 160); ys = sol.sol(te); xs = ys[0]
        for i in range(len(te) - 1):
            if xs[i]*xs[i+1] < 0:
                tc = brentq(lambda tt: sol.sol(tt)[0], te[i], te[i+1], xtol=1e-12*tg); cr.append(tc)
                Ed.append(abs(K['energy'](sol.sol(tc)) - E0)/E0)
        T += list(te); X += list(xs); Np += [K['npairs'](ys[:, i]) for i in range(0, len(te), 20)]
        y = sol.y[:, -1]; t = sol.t[-1]; ff = y[2:]; fmax = max(fmax, float(ff.max())); fmin = min(fmin, float(ff.min()))
        Pcov = max(Pcov, abs(y[1]))
        if not np.isfinite(y).all(): break
    ff = y[2:].reshape(NP, J); edge = float(max(abs(ff[:5]).max(), abs(ff[-5:]).max()))
    T = np.array(T); X = np.array(X); cr = np.array(cr)
    amps = [float(np.abs(X[(T > cr[i]) & (T < cr[i+1])]).max()/x0) for i in range(len(cr) - 1)] if len(cr) > 1 else []
    return dict(cr=cr, Ed=Ed, T=T, X=X, amps=amps, fmax=fmax, fmin=fmin, edge_f=edge, Phimax=float(Pcov), Pmax=float(K['P'][-1]), NP=NP, E0=E0,
                npairs_end=K['npairs'](y), runtime=time.time() - tstart, x0=x0)

def analyze(r):
    cr = r['cr']; hp = np.diff(np.concatenate(([0.0], cr)))                      # hp[0]=tiempo al 1er cero; luego semiperiodos
    half = np.diff(cr); w = np.pi/half if len(half) else np.array([])
    out = dict(n_cross=int(len(cr)), tau0=float(cr[0]) if len(cr) else None, half_periods=[float(v) for v in half], amps=r['amps'],
               omega_early=float(w[0]) if len(w) else None, omega_late=float(np.mean(w[-3:])) if len(w) >= 3 else None, Emax_err=float(max(r['Ed'])) if r['Ed'] else None,
               fmax=r['fmax'], fmin=r['fmin'], edge_f=r['edge_f'], Phimax_over_Pmax=r['Phimax']/r['Pmax'], npairs_end=r['npairs_end'], runtime=r['runtime'])
    if len(r['amps']) >= 3:
        tm = 0.5*(cr[:-1] + cr[1:]); out['gamma'] = float(-np.polyfit(tm, np.log(r['amps']), 1)[0])
    return out

sweep = {}; curves = {}
X0S = [30.0, 100.0, 1e3, 1e4]
for x0 in X0S:
    r = run(x0, ncross=10, block=True); a = analyze(r); sweep[f"{x0:g}"] = a; curves[x0] = r
    a['omega_over_esqrtx0'] = a['omega_late']/(E_CHARGE*np.sqrt(x0)) if a.get('omega_late') else None
    a['gamma_over_omega'] = a['gamma']/a['omega_late'] if a.get('gamma') and a.get('omega_late') else None
    a['tau0_e_sqrtx0'] = a['tau0']*E_CHARGE*np.sqrt(x0)
    print(f"  ε0={x0:>6g}: n_cross={a['n_cross']} τ0·e√ε0={a['tau0_e_sqrtx0']:.3f} ω_late/(e√ε0)={a['omega_over_esqrtx0']:.4f} γ/ω={a['gamma_over_omega']} max|ΔE|/E0={a['Emax_err']:.1e} fmax={a['fmax']:.6f} ({a['runtime']:.0f}s)")
R['sweep'] = sweep

# --- sin bloqueo (para cuantificar su efecto y regresion con 05) ---
rnb = run(100.0, ncross=6, block=False); anb = analyze(rnb)
R['no_block_x100'] = dict(tau0=anb['tau0'], fmax=anb['fmax'], omega_early=anb['omega_early'], amps=anb['amps'])
D = load_all(); tau0_05 = D['branch_B']['sweep']['100']['tau0']
gate("K3_regresion_con_05", abs(anb['tau0']/tau0_05 - 1) < 0.01, f"sin bloqueo: τ0 = {anb['tau0']:.5f} vs modelo de cohortes (05) {tau0_05:.5f}: diferencia {abs(anb['tau0']/tau0_05-1):.1e}; con bloqueo τ0 = {sweep['100']['tau0']:.4f} ({sweep['100']['tau0']/anb['tau0']-1:+.1%})")
R['pauli_effect_x100'] = dict(tau0_block=sweep['100']['tau0'], tau0_noblock=anb['tau0'], f_max_noblock=anb['fmax'], f_max_block=sweep['100']['fmax'])
fmx = max(v['fmax'] for v in sweep.values()); fmn = min(v['fmin'] for v in sweep.values())
gate("K4_0<=f<=1_Pauli", fmn >= -1e-12 and fmx <= 0.5 + 1e-9 and anb['fmax'] > 1.0, f"con (1−2f): {fmn:.1e} ≤ f ≤ {fmx:.9f} (<1/2 ≤ 1) en todas las mallas; SIN bloqueo f llega a {anb['fmax']:.2f} (>1: violación de Pauli que el factor evita)")
ncr = min(v['n_cross'] for v in sweep.values()); eerr = max(v['Emax_err'] for v in sweep.values())
gate("K5_energia_5_oscilaciones", ncr >= 10 and eerr < 1e-4, f"≥{ncr} cruces por cero (≥5 oscilaciones completas); máx |ΔE|/E0 = {eerr:.2e} < 1e-4 (DOP853 rtol=1e-9)")
gate("K5b_cobertura_malla", max(v['Phimax_over_Pmax'] for v in sweep.values()) < 0.9 and max(v['edge_f'] for v in sweep.values()) < 1e-6, f"|Φ|max/P_max = {max(v['Phimax_over_Pmax'] for v in sweep.values()):.2f}; f en los bordes ≤ {max(v['edge_f'] for v in sweep.values()):.1e}: la malla contiene toda la dinámica")

# --- convergencia / sensibilidad ---
base = sweep['100']; conv = {}
for name, kw in (('rtol1e-11', dict(rtol=1e-11)), ('dp/1.5', dict(dpf=1/5)), ('J=14', dict(J=14))):
    a = analyze(run(100.0, ncross=6, block=True, **kw)); conv[name] = dict(tau0=a['tau0'], omega_early=a['omega_early'], Emax_err=a['Emax_err'])
dev = max(abs(v['tau0']/base['tau0'] - 1) for v in conv.values()); devw = max(abs(v['omega_early']/base['omega_early'] - 1) for v in conv.values())
conv['max_rel_change_tau0'] = float(dev); conv['max_rel_change_omega'] = float(devw); R['convergence'] = conv
gate("K6_convergencia_numerica", dev < 2e-3 and devw < 3e-3, f"τ0 y ω bajo rtol 1e-9→1e-11, dp÷1.5 y J=10→14: máx. variación {dev:.1e} / {devw:.1e}")
sens = {}
for sf in (0.15, 0.40):
    a = analyze(run(100.0, ncross=6, block=True, sigf=sf)); sens[f"sigma_{sf}"] = dict(tau0=a['tau0'], omega_early=a['omega_early'])
R['sigma_sensitivity'] = sens
sdev = max(abs(v['tau0']/base['tau0'] - 1) for v in sens.values())
R['sigma_sensitivity']['max_rel_change_tau0'] = float(sdev)
print(f"  sensibilidad al ancho de la fuente σ∈{{0.15,0.25,0.40}}√ε0: Δτ0/τ0 ≤ {sdev:.1%} (parámetro de MODELO, no de discretización)")

# --- ley de escala de la frecuencia de plasma y del amortiguamiento ---
lx = np.log(X0S[1:]); lw = np.log([sweep[f"{v:g}"]['omega_late'] for v in X0S[1:]])
expo = float(np.polyfit(lx, lw, 1)[0])
cs = [sweep[f"{v:g}"]['omega_over_esqrtx0'] for v in X0S[1:]]
gs = [sweep[f"{v:g}"]['gamma_over_omega'] for v in X0S[1:]]
R['scaling'] = dict(omega_exponent=expo, omega_over_esqrtx0=[float(c) for c in cs], gamma_over_omega=[float(g) for g in gs], x0s=X0S[1:])
gate("K7_escala_omega_p", abs(expo - 0.5) < 0.03 and (max(cs) - min(cs))/np.mean(cs) < 0.05, f"ω_p ∝ ε0^{expo:.3f} (esperado 1/2); ω_p/(e√ε0) = {cs} (ε0∈[100,10⁴]); γ/ω_p = {[round(g,4) for g in gs]}")

# --- controles negativos numericos ---
cb = analyze(run(100.0, ncross=4, block=True, polblock=False))
gate("C6_control_jpol_sin_(1-2f)", cb['Emax_err'] > 1e-3 and cb['Emax_err'] > 100*eerr, f"j_pol sin (1−2f): |ΔE|/E0 = {cb['Emax_err']:.2e} ≫ {eerr:.1e} (detectado)")
cw = run(100.0, ncross=3, block=True, wrongsign=True, tmax_factor=15); cwa = analyze(cw)
gate("C7_control_(1+2f)", cw['fmax'] > 0.5 + 1e-3, f"con (1+2f) la ocupación llega a f = {cw['fmax']:.3f} > 1/2 (el bloqueo de Pauli se convierte en estimulación; detectado)")
cf = analyze(run(100.0, ncross=3, block=True, flip=True, tmax_factor=15))
gate("C8_control_signo_conduccion", cf['Emax_err'] is None or cf['Emax_err'] > 1e-3 or cf['n_cross'] == 0, f"corriente de conducción con signo erróneo: |ΔE|/E0 = {cf['Emax_err']}, cruces = {cf['n_cross']} (la identidad de energía se rompe)")

# --- figura vectorial ---
try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
    for x0, c in ((100.0, 'C0'), (1e4, 'C1')):
        r = curves[x0]; ts = r['T']*E_CHARGE*np.sqrt(x0)
        ax[0].plot(ts, r['X']/x0, color=c, lw=1.2, label=rf'$\epsilon_0={x0:g}$ (Pauli blocking)')
    ts = rnb['T']*E_CHARGE*np.sqrt(100.0); ax[0].plot(ts, rnb['X']/100.0, 'k--', lw=1.0, label=r'$\epsilon_0=100$, no blocking')
    ax[0].set_xlabel(r'$t\,e\sqrt{\epsilon_0}$'); ax[0].set_ylabel(r'$\epsilon(t)/\epsilon_0$'); ax[0].axhline(0, color='gray', lw=0.5); ax[0].legend(fontsize=8)
    ax[0].set_title('Damped plasma oscillations')
    xx = np.array(X0S[1:]); ax[1].plot(xx, [sweep[f"{v:g}"]['omega_over_esqrtx0'] for v in xx], 'o-', label=r'$\omega_p/(e\sqrt{\epsilon_0})$')
    ax[1].plot(xx, [sweep[f"{v:g}"]['gamma_over_omega'] for v in xx], 's-', label=r'$\gamma/\omega_p$'); ax[1].set_xscale('log'); ax[1].set_xlabel(r'$\epsilon_0$'); ax[1].legend(); ax[1].set_title('Scalings with $\\epsilon_0$')
    fig.tight_layout(); fig.savefig('fig_plasma_oscillations.pdf'); fig.savefig('fig_plasma_oscillations.png', dpi=160); R['figure'] = 'fig_plasma_oscillations.pdf'
except Exception as ex:
    R['figure'] = f'no generada: {ex}'

R['controls'] = dict(C6_energy_err=float(cb['Emax_err']), C7_fmax=float(cw['fmax']), C8_energy_err=cf['Emax_err'], C8_n_cross=int(cf['n_cross']))
R['model'] = "Markov + LCFA + fuente gaussiana de ancho σ=0.25√ε0 + factor (1-2f); campo homogeneo; sin recoleccion no markoviana; f<1/2 desde el vacio"
save_section('branch_B2', R)
print("\nRESUMEN 05b:", "TODAS PASAN" if ok_all else "HAY FALLOS")
sys.exit(0 if ok_all else 1)
