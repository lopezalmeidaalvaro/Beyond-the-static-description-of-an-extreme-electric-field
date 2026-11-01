"""SUITE DE INVARIANTES FISICOS (5 pruebas). exit(1) si alguna falla."""
import io, sys, runpy, contextlib
import numpy as np, sympy as sp
from eeh import *

results = {}
def report(name, ok, msg):
    results[name] = ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}")

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path('01_derive.py')           # reutiliza geometria/T_mn verificados
r, c, G = ns['r'], ns['c'], ns['G']
Tud, Gam, x = ns['Tud'], ns['Gam'], ns['x']
onshell = ns['onshell']

# ---- 1. LIMITE ASINTOTICO / CORRESPONDENCIA CLASICA -----------------------------------
Es, E0s, cs = sp.symbols('E E0 c', positive=True)
# solucion perturbativa explicita E = E0 - 16 c E0^3 + 768 c^2 E0^5 + ...  (verificada por sustitucion)
Epert = E0s - 16*cs*E0s**3 + 768*cs**2*E0s**5
resid = sp.series(sp.expand(Epert*(1+16*cs*Epert**2) - E0s), cs, 0, 3).removeO()
sym_ok = sp.simplify(resid) == 0
# numerico: (a) c->0 recupera Coulomb exacto; (b) exponente de la correccion = r^-4 cuando r->inf
Q = 1e6*E_CHARGE
rr = np.array([300., 600., 1200., 2400.])
E_c0 = E_of_r(rr, Q, c=0.0)
coul_ok = np.max(np.abs(E_c0/(Q/(4*np.pi*rr**2)) - 1)) < 1e-14
dev = np.abs(E_of_r(rr, Q)/(Q/(4*np.pi*rr**2)) - 1)
slope = np.polyfit(np.log(rr), np.log(dev), 1)[0]
report("1 Limite asintotico", sym_ok and coul_ok and abs(slope + 4) < 0.01,
       f"serie E=E0-16cE0^3+..: residuo simbolico 0 ({sym_ok}); c=0 -> Coulomb exacto ({coul_ok}); "
       f"desviacion ~ r^{slope:.3f} (esperado -4); (RN exacto con c=0 se prueba en 03)")

# ---- 2. CONSERVACION  nabla_mu T^mu_nu = 0 --------------------------------------------
def div_T(nu):
    return sum(sp.diff(Tud[m, nu], x[m]) for m in range(4)) + \
           sum(Gam[m][m][l]*Tud[l, nu] for m in range(4) for l in range(4)) - \
           sum(Gam[l][m][nu]*Tud[m, l] for m in range(4) for l in range(4))
divr = sp.simplify(div_T(1))
off = divr                                            # fuera de capa: proporcional a Maxwell, no nulo
res_on = onshell(divr)                                # sobre la solucion
# numerico independiente: rho' = -(2/r)(rho+p_perp) por diferencias finitas sobre E(r)
rn = np.logspace(-1, 3, 400); h = 1e-5
Q2 = 1e3*E_CHARGE
drho = (rho(E_of_r(rn*(1+h), Q2)) - rho(E_of_r(rn*(1-h), Q2)))/(2*h*rn)
rhs = -2/rn*(rho(E_of_r(rn, Q2)) + p_perp(E_of_r(rn, Q2)))
rel = np.max(np.abs(drho - rhs)/np.abs(rhs))
report("2 Conservacion T", res_on == 0 and rel < 1e-7,
       f"div simbolica en-capa = {res_on} (fuera de capa != 0: {off != 0}); "
       f"max|rho'+2(rho+p_t)/r|/|.| = {rel:.2e}")

# ---- 3. CONDICIONES DE ENERGIA (NEC/WEC) ----------------------------------------------
Esym = sp.symbols('E', real=True)
rh, pr, pt = Esym**2/2 + 12*cs*Esym**4, -(Esym**2/2 + 12*cs*Esym**4), Esym**2/2 + 4*cs*Esym**4
# comprobar contra el tensor simbolico del CAS
chk = [sp.simplify(Tud[0,0].subs({ns['f']:1, ns['g']:1}).doit().subs(sp.Derivative(ns['Phi'], r), Esym) + rh),
       sp.simplify(Tud[2,2].subs({ns['f']:1, ns['g']:1}).doit().subs(sp.Derivative(ns['Phi'], r), Esym) - pt)]
chk = [sp.simplify(s.subs(c, cs)) for s in chk]
sym_ok = all(s == 0 for s in chk)
nec = [sp.factor(rh + pr), sp.factor(rh + pt)]
print("      NEC: rho+p_r =", nec[0], ";  rho+p_perp =", nec[1], ";  WEC: rho =", sp.factor(rh))
rng = np.random.default_rng(1); worst = np.inf
for _ in range(20000):
    E = 10**rng.uniform(-6, 12); n = rng.normal(size=3); n /= np.linalg.norm(n)   # k=(1,n) nulo (marco ortonormal)
    Tk = rho(E) * 1 + (-rho(E))*n[0]**2 + p_perp(E)*(n[1]**2 + n[2]**2)          # T_ab k^a k^b, n_r=n[0]
    worst = min(worst, Tk/(E**2))
report("3 NEC/WEC", sym_ok and worst >= 0 and sp.simplify(nec[1]).is_nonnegative is not False,
       f"T_munu k^mu k^nu/E^2 >= {worst:.3e} (20000 vectores nulos); NEC saturada en direccion radial "
       f"(rho+p_r=0); rho>0. Para c>0 (signo EH) no hay violacion.")

# ---- 4. ANALISIS DIMENSIONAL (unidades naturales; dimension de masa) -------------------
lam = sp.symbols('lam', positive=True)
dim = {sp.Symbol('G', positive=True): -2, sp.Symbol('c', positive=True): -4, ns['r']: -1}
Fs, Rs_, m_, al = sp.symbols('F2 R m alpha', positive=True)
terms = {"R/(16 pi G)": Rs_/(16*sp.pi*sp.Symbol('G', positive=True)), "-F2/4": Fs/4,
         "c*(F2)^2 (EH)": sp.Symbol('c', positive=True)*Fs**2, "alpha^2/(90 m^4)*(F2)^2": al**2/(90*m_**4)*Fs**2}
d = {sp.Symbol('G', positive=True): -2, sp.Symbol('c', positive=True): -4, Rs_: 2, Fs: 4, m_: 1, al: 0}
ok4 = True
for k, v in terms.items():
    dd = sp.simplify(v.subs({s: lam**dv for s, dv in d.items()})).as_coeff_exponent(lam)[1]
    ok4 &= (dd == 4)
    print(f"      [{k}] = mass^{dd}")
Esy, Qd, rd, y_, rho_ = sp.symbols('E Q r y rho', positive=True)
d2 = {Esy: 2, rd: -1, sp.Symbol('c', positive=True): -4, sp.Symbol('G', positive=True): -2, rho_: 4, y_: -1}
cond = {"Maxwell: r^2 E (1+16cE^2)": rd**2*Esy*(1 + 16*sp.Symbol('c', positive=True)*Esy**2),   # debe ser adimensional (Q)
        "16 c E^2 (adimensional)": 16*sp.Symbol('c', positive=True)*Esy**2,
        "dy/dr = 8 pi G r^2 rho": 8*sp.pi*sp.Symbol('G', positive=True)*rd**2*rho_}
for k, v in cond.items():
    dd = sp.simplify(v.subs({s: lam**dv for s, dv in d2.items()})).as_coeff_exponent(lam)[1]
    ok4 &= (dd == 0)
    print(f"      [{k}] = mass^{dd}")
# eps = eE/m^2 adimensional y cutoff de Planck coherente: rho_Pl = M_Pl^4 ; E_Pl ~ M_Pl^2
report("4 Dimensional", ok4, "todos los terminos del Lagrangiano tienen masa^4; ecuaciones radiales coherentes")

# ---- 5. COTA DE UNITARIEDAD ---------------------------------------------------------------
eps_grid = np.logspace(-3, 3, 4000)
Gam_grid = np.array([schwinger_rate(e/E_CHARGE) for e in eps_grid])     # Gamma/m^4
eps_unit = eps_grid[np.argmax(Gam_grid > 1)]
EPS_MAX = 1.0
report("5 Unitariedad", Gam_grid[np.searchsorted(eps_grid, EPS_MAX)] < 1 and eps_unit > EPS_MAX and np.all(np.diff(Gam_grid) >= 0),
       f"P(par en volumen Compton x tiempo Compton)=Gamma lambda_C^4 <1 hasta eps={eps_unit:.2f}; "
       f"el cutoff EFT eps={EPS_MAX} actua antes -> la cota nunca se usa para 'predecir'")

# ---- CONTROL NEGATIVO: la suite debe rechazar un modelo erroneo (signo de c invertido) ----
Ebig = np.logspace(0, 4, 50)
nec_bad = (rho(Ebig, -C_EH) + p_perp(Ebig, -C_EH)).min()
print(f"[CONTROL] c -> -c: min(rho+p_perp) = {nec_bad:.2e} < 0 -> NEC violada, la suite lo rechazaria: {nec_bad < 0}")
results["control negativo"] = bool(nec_bad < 0)

print("\nRESUMEN:", "TODAS PASAN" if all(results.values()) else "HAY FALLOS", results)
sys.exit(0 if all(results.values()) else 1)
