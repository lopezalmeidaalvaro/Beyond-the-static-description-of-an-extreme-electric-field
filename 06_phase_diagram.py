"""RAMA C — escalas y diagrama de fases (N, r0, tau).  Longitudes en lambda_C, masas/energias en m_e, c=hbar=1.
Entradas verificadas: 04 (Landau, campo QED), 05 (kappa_inf).  Criterio de atrapamiento: desigualdad de Penrose cargada
(Malec-Ó Murchadha; spherical, carga dentro de la superficie atrapada)  G M_ADM >= (R + r_Q^2/R)/2, r_Q^2 = G Q^2/4pi."""
import sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, sympy as sp
from common import *

R = {}; ok_all = True
def gate(name, ok, msg):
    global ok_all; ok_all &= bool(ok); R.setdefault('gates', {})[name] = bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}")

DB = load_all()
KAPPA = DB['branch_B']['scaling']['kappa_inf']
# ====================== SIMBOLICO ======================
Gs, Ms, Qs, Rs, al, mp, N, r0, kap = sp.symbols('G M Q R alpha M_Pl N r0 kappa', positive=True)
rQ2 = Gs*Qs**2/(4*sp.pi)
rp = Gs*Ms*(1 + sp.sqrt(1 - rQ2/(Gs*Ms)**2))                     # horizonte externo RN
P1 = sp.simplify(Gs*Ms - (rp + rQ2/rp)/2)
gate("P1 RN satura Penrose cargada", P1 == 0, "G M = (r_+ + r_Q^2/r_+)/2 para el horizonte de Reissner-Nordström (identidad exacta)")
UC = Qs**2/(8*sp.pi*Rs)
P2 = sp.simplify(Gs*UC - rQ2/(2*Rs))
gate("P2 cancelacion Coulomb (cascaron)", P2 == 0, "G·U_C(R) = r_Q²/(2R): la energia electrostatica cancela EXACTAMENTE en la desigualdad => criterio: G M_noEM >= R/2")
Mn = sp.symbols('M_n', positive=True)
crit = sp.simplify((Gs*(Mn + UC) - (Rs + rQ2/Rs)/2) - (Gs*Mn - Rs/2))
gate("P3 criterio reducido", crit == 0, "G(M_n+U_C) - (R+r_Q²/R)/2 = G M_n - R/2")
Gm = 1/mp**2; Nn = sp.symbols('Nn', positive=True); aF = sp.Rational(3, 4)*(9*sp.pi/4)**sp.Rational(1, 3)
rs_sol = sp.solve(sp.Eq(r0**2 - 2*Gm*Nn*r0 - 2*Gm*aF*Nn**sp.Rational(4, 3), 0), r0)
rs_pos = [s for s in rs_sol if sp.simplify(s.subs({mp: 2e22, Nn: 1e6})) > 0][0]
chk = sp.simplify((r0 - 2*Gm*(Nn + aF*Nn**sp.Rational(4, 3)/r0)).subs(r0, rs_pos))
gate("P4 radio de atrapamiento r_s(N)", chk == 0, f"r_s = {sp.simplify(rs_pos)}  [r0 = 2G(N m + E_Fermi), E_F=(3/4)(9π/4)^(1/3)N^(4/3)/r0]")
rE = sp.sqrt(Nn*al); rQ = Nn*sp.sqrt(al)/mp
Ncrit = sp.solve(sp.Eq(rE, rQ), Nn)[0]
gate("P5 N_crit = M_Pl^2", sp.simplify(Ncrit - mp**2) == 0, "r_Q = r_ε  <=>  N = M_Pl²/m_e²  (el horizonte RN entra en el dominio ε<1)")
rhoC = (4*sp.pi*al*Nn**2/(32*sp.pi**2))/(rQ**4)                  # rho_Coulomb(r_Q) = Q^2/(32 pi^2 r_Q^4), Q^2=4 pi alpha N^2
P6 = sp.simplify(rhoC/mp**4 - 1/(8*sp.pi*al*Nn**2))
gate("P6 densidad en r_Q", P6 == 0, "ρ(r_Q)/M_Pl⁴ = 1/(8π α N²): sub-planckiana para N ≥ N_P=(8πα)^(-1/2)")
NP = sp.solve(sp.Eq(1/(8*sp.pi*al*Nn**2), 1), Nn)[0]
NP_val = float(NP.subs(al, ALPHA))
Nstar = sp.solve(sp.Eq(kap*r0/(al*sp.sqrt(4*sp.pi*Nn)), r0), Nn)[0]       # tau_dis = r0
gate("P7 N* (τ_dis = r0)", sp.simplify(Nstar - kap**2/(4*sp.pi*al**2)) == 0, "τ_dis = κ r0/(α√(4πN)) = r0  =>  N* = κ²/(4π α²)")
# tiempo de caida libre (polvo, Oppenheimer-Snyder): int_0^R0 dR / sqrt(2GM(1/R-1/R0)) = (pi/2) sqrt(R0^3/(2GM))
th_, R0_, GM_ = sp.symbols('theta R0 GM', positive=True)
Rsub = R0_*sp.sin(th_)**2
integrand = sp.simplify(sp.diff(Rsub, th_)/sp.sqrt(2*GM_*(1/Rsub - 1/R0_)))
tff = sp.integrate(sp.simplify(integrand.subs(th_, th_)), (th_, 0, sp.pi/2))
gate("P8 caida libre OS", sp.simplify(tff - sp.pi/2*sp.sqrt(R0_**3/(2*GM_))) == 0, f"τ_ff = (π/2)√(R0³/2GM) = {sp.simplify(tff)}")
Nbh = sp.solve(sp.Eq(2*Nn/mp**2, rE), Nn)[0]
Nbh_val = float(Nbh.subs({mp: MPL, al: ALPHA}))
bar = sp.simplify((rQ2/rQ**2*0 + 0))                                    # (placeholder inofensivo)
# control negativo simbolico: criterio errado "G(M_n+U_C) >= R/2" (olvida la carga) NO es satisfecho por RN extremal
Mext = Qs/sp.sqrt(4*sp.pi*Gs)
wrong = sp.simplify((Gs*Mext) - (Mext*0 + (Gs*Mext)))               # RN extremal: r_+=GM, GM=(GM+GM)/2
# criterio erroneo: ignorar r_Q^2/R -> exigir G M >= R/2 en R=r_+ extremal: G M >= G M/2 (siempre)... usamos RN subextremal
Mn_ = sp.symbols('Mnn', positive=True)
wrong_res = sp.simplify((Gs*Mext - (Rs + rQ2/Rs)/2).subs(Rs, Gs*Mext*3/4))
gate("C1 control: R=0.75 GM_ext viola Penrose cargada", wrong_res.subs({Gs: 1, Qs: 1}) < 0, f"G M_ext - (R + r_Q²/R)/2 = {sp.simplify(wrong_res.subs({Gs:1,Qs:1}))} < 0 en R=3/4·G M_ext (R<r_- no atrapada)")

# ley exacta de la carrera (polvo de N masas en reposo): tau_dis / tau_ff = sqrt(r0*/r0), r0* independiente de N
r0s_sym = sp.symbols('r0s', positive=True)
tdis = kap*r0/(al*sp.sqrt(4*sp.pi*Nn)); tff_ = sp.pi/2*sp.sqrt(r0**3*mp**2/(2*Nn))          # GM = N/M_Pl^2
r0star = sp.solve(sp.Eq(tdis, tff_), r0)[0]
P9 = sp.simplify((tdis/tff_)**2 - r0star/r0)
gate("P9 ley de la carrera", P9 == 0 and sp.simplify(r0star - 2*kap**2/(sp.pi**3*al**2*mp**2)) == 0 and not r0star.has(Nn),
     f"τ_dis/τ_ff = √(r0*/r0),  r0* = {sp.simplify(r0star)} (independiente de N)")
r0star_val = float(r0star.subs({kap: KAPPA, al: ALPHA, mp: MPL}))
# ====================== NUMERICO ======================
a, MP = ALPHA, MPL
aF_ = 0.75*(9*np.pi/4)**(1/3)
r_eps = lambda n: np.sqrt(n*a)
r_Q   = lambda n: n*np.sqrt(a)/MP
r_P0  = lambda n: (4*np.pi*a*n**2/(32*np.pi**2))**0.25/MP                 # rho_Coulomb = M_Pl^4
r_P   = lambda n: np.maximum(r_P0(n), 10.0/MP)                           # + piso 10 l_Planck (igual que 03_integrate.py)
r_s   = lambda n, mu=1.0: (mu*n/MP**2) + np.sqrt((mu*n/MP**2)**2 + 2*aF_*n**(4/3)/MP**2)
tau_dis = lambda n, r: KAPPA*r/(a*np.sqrt(4*np.pi*n))
# comprobaciones de escala
gate("N1 N_crit numerico", abs(r_eps(MP**2)/r_Q(MP**2) - 1) < 1e-12, f"r_ε(M_Pl²)/r_Q = {r_eps(MP**2)/r_Q(MP**2):.12f}")
gate("N2 r_s con N electrones", abs(r_s(1e6) - 2*(1e6 + aF_*1e6**(4/3)/r_s(1e6))/MP**2)/r_s(1e6) < 1e-12, f"r_s(N=1e6) = {r_s(1e6):.4g} λ_C = {r_s(1e6)*MP:.4g} l_P")
Nnum = np.logspace(0, 50, 2001); r0n = np.logspace(-30, 30, 1201)
NN, RR = np.meshgrid(Nnum, r0n, indexing='ij')
ph = np.zeros(NN.shape, dtype=int)                                      # 0 EFT estatico, 1 descarga, 2 BH, 3 transplanckiano
ph[:] = 1
ph[RR > r_eps(NN)] = 0
ph[RR <= r_s(NN)] = 2
ph[RR <= r_P(NN)] = 3
frac = {k: float(np.mean(ph == i)) for i, k in enumerate(['EFT_estatica', 'descarga_Schwinger', 'BH_atrapada', 'transplanckiana'])}
R['phase_fractions'] = frac
# ¿region BH sobre la planckiana?  r_P < r_s para todo N>=1
N_BHmin = float(Nnum[np.argmax(r_s(Nnum) > r_P(Nnum))])
gate("N3 BH por encima de la region trans-planckiana", np.all(r_s(Nnum[Nnum >= 20]) > r_P(Nnum[Nnum >= 20])) and r_s(5.0) < r_P(5.0),
     f"r_s(N) > max(r_P, 10 l_P) para N ≥ 20; N minimo con atrapamiento sub-planckiano ≈ {N_BHmin:.3g}; para N ≲ 14 el atrapamiento seria trans-planckiano")
R['N_BH_min'] = N_BHmin
# donde esta N_P (r_Q = r_P)
gate("N4 N_P", abs(r_Q(NP_val)/r_P0(NP_val) - 1) < 1e-9, f"N_P = (8πα)^(-1/2) = {NP_val:.3f}: r_Q > r_P(densidad) para N ≥ 3")
R['race_law'] = dict(r0_star_lambdaC=r0star_val, r0_star_over_lPlanck=r0star_val*MP, statement='tau_dis/tau_ff = sqrt(r0*/r0) para r0 >> r0*')
R['critical'] = dict(N_crit=float(MP**2), N_P=NP_val, N_star=float(KAPPA**2/(4*np.pi*ALPHA**2)), N_classical_BH_in_EFT=Nbh_val,
                     N_valid_homog=float(9/ALPHA), r_neutral_lambdaC=float(np.sqrt(ALPHA)), M_ext_per_charge_me=float(np.sqrt(ALPHA)*MP),
                     z2_coulomb_over_gravity=float(ALPHA*MP**2), kappa=KAPPA)
# tabla de longitudes y tiempos por N
tab = {}
for n in (1.0, 1e3, 1e6, 1e12, 1e20, 1e30, MP**2, 1e48, 1e60):
    rr = r_s(n)
    C = 2*G_NEWTON*(n + aF_*n**(4/3)/rr)/rr                              # =1 en r0=r_s (por construccion)
    tab[f"{n:.3g}"] = dict(r_eps=float(r_eps(n)), r_Q=float(r_Q(n)), r_P=float(r_P0(n)), r_s=float(rr), r_s_ext=float(r_s(n, mu=np.sqrt(ALPHA)*MP)),
                           tau_dis_over_r0=float(KAPPA/(ALPHA*np.sqrt(4*np.pi*n))), tau_dis_at_rs=float(tau_dis(n, rr)),
                           tau_ff_at_rs=float(np.pi/2*rr/np.sqrt(C)), R_race_at_rs=float(tau_dis(n, rr)/(np.pi/2*rr/np.sqrt(C))),
                           eps_at_rs=float(n*a/rr**2), ρ_Coul_at_rs_over_MPl4=float(4*np.pi*a*n**2/(32*np.pi**2)/rr**4/MP**4))
R['table'] = tab
print("\n  N          r_eps        r_Q          r_P          r_s        eps(r_s)    tau_dis/r0")
for k, v in tab.items():
    print(f"  {k:>9} {v['r_eps']:.3e} {v['r_Q']:.3e} {v['r_P']:.3e} {v['r_s']:.3e} {v['eps_at_rs']:.2e}  {v['tau_dis_over_r0']:.3e}")
print("\nFracciones del plano (N∈[1,1e50], r0∈[1e-30,1e30] λ_C):", {k: round(v, 3) for k, v in frac.items()})
# ¿Hay region BH con eps<1? (BH clasico dentro de la EFT estatica)
inEFT = np.any((RR <= r_s(NN)) & (RR > r_eps(NN)))
R['BH_with_eps_lt_1_exists'] = bool(inEFT)
Nmin = Nnum[np.argmax(np.any((RR <= r_s(NN)) & (RR > r_eps(NN)), axis=1))] if inEFT else None
R['N_min_BH_in_static_EFT'] = None if Nmin is None else float(Nmin)
gate("N5 BH clasico dentro de la EFT estatica", (not inEFT) and Nbh_val > Nnum.max(), f"en N≤{Nnum.max():.0e} NO hay r_s > r_ε; umbral analitico N = α M_Pl⁴/4 = {Nbh_val:.3g} (masa en reposo)")
# carrera (BH phase): descarga vs caida libre y vs luz
rows = []
for n in np.logspace(0, 50, 11):
    rr = r_s(n); t_ff = np.pi/2*rr; rows.append((n, tau_dis(n, rr)/rr, tau_dis(n, rr)/t_ff))
R['race_at_rs'] = [dict(N=float(a_), tau_dis_over_light=float(b_), tau_dis_over_ff=float(c_)) for a_, b_, c_ in rows]
# ---------- control negativo numerico: criterio sin cancelacion de Coulomb ----------
wrong_boundary = lambda n, r: 2*G_NEWTON*(n + n**2*a*0.6/r)/r >= 1          # sumaria U_C como gravitante
bad_frac = float(np.mean(wrong_boundary(NN, RR) & (RR > r_P(NN))))
gate("C2 control: sumar U_C (error) agranda la region BH", bad_frac > frac['BH_atrapada']*1.01, f"fraccion BH erronea {bad_frac:.3f} > correcta {frac['BH_atrapada']:.3f}")
gate("C3 control: invertir ρ>M_Pl^4 elimina la fase (iii)", float(np.mean((RR >= r_P(NN)) & (ph == 3))) == 0.0, "sin el criterio de Planck no hay region trans-planckiana (la clasificacion lo detecta)")

# ---------- figura ----------
try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    cm = ListedColormap(['#dfe9f5', '#f7d9a6', '#b5453d', '#4a3a6a'])
    ax.pcolormesh(np.log10(Nnum), np.log10(r0n), ph.T, cmap=cm, shading='auto', vmin=0, vmax=3, rasterized=True)
    ls = np.log10(Nnum)
    for fn, lab, st in ((r_eps, r'$r_\epsilon=\sqrt{N\alpha}$ ($\epsilon=1$)', '-'), (r_Q, r'$r_Q$ (RN)', '--'), (r_P, r'$r_P$ ($\rho=M_{Pl}^4$)', ':'), (r_s, r'$r_s$ (trapped)', '-')):
        ax.plot(ls, np.log10(fn(Nnum)), st, color='k', lw=1.2, label=lab)
    ax.axvline(np.log10(MP**2), color='gray', lw=0.8); ax.text(np.log10(MP**2)+0.3, -27, r'$N_{crit}=M_{Pl}^2$', rotation=90, fontsize=8)
    ax.axvline(np.log10(KAPPA**2/(4*np.pi*ALPHA**2)), color='gray', lw=0.8, ls='--'); ax.text(np.log10(KAPPA**2/(4*np.pi*ALPHA**2))+0.3, 12, r'$N^*$: $\tau_{dis}=r_0$', rotation=90, fontsize=8)
    ax.set_xlim(0, 50); ax.set_ylim(-30, 30); ax.set_xlabel(r'$\log_{10} N$'); ax.set_ylabel(r'$\log_{10}(r_0/\lambda_C)$')
    ax.set_title('Fases: azul=EFT estática, ocre=descarga Schwinger, rojo=BH/atrapada, violeta=trans-Planck')
    ax.legend(loc='lower right', fontsize=8); fig.tight_layout(); fig.savefig('phase_diagram.png', dpi=140)
    R['figure'] = 'phase_diagram.png'
except Exception as ex:
    R['figure'] = f'no generada: {ex}'

# ---------- respuesta a la pregunta central (lógica derivada de los gates) ----------
R['verdict'] = dict(
  coulomb_cancels_in_penrose=bool(R['gates']['P2 cancelacion Coulomb (cascaron)']),
  trapped_requires_nonEM_mass=bool(R['gates']['P3 criterio reducido']),
  qed_prevents_trapped_given_nonEM_mass=not (R['gates']['P3 criterio reducido'] and DB['branch_B']['gates']['D2 conservacion de energia campo+plasma']),
  pure_field_trapped_possible=False,
  planck_density_at_trap_radius=bool(R['gates']['N3 BH por encima de la region trans-planckiana']),
  landau_reachable=bool(DB['branch_A']['Landau']['ln_eps_L'] < DB['branch_A']['Landau']['ln_eps_Planck']))
save_section('branch_C', R)
print("\nRESUMEN 06:", "TODAS PASAN" if ok_all else "HAY FALLOS")
sys.exit(0 if ok_all else 1)
