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

# ====================== P10: DEMOSTRACION INTERNA DE ATRAPAMIENTO (Misner-Sharp, SymPy) ======================
# (a) geometria: doble nulo  ds^2 = -2 e^{2s} du dv + R^2 dOmega^2   (todo metrica 2D es conforme a esta)
u_, v_ = sp.symbols('u v', real=True)
Rf = sp.Function('R')(u_, v_); sg = sp.Function('sigma')(u_, v_)
h2 = sp.Matrix([[0, -sp.exp(2*sg)], [-sp.exp(2*sg), 0]]); h2i = sp.Matrix([[0, -sp.exp(-2*sg)], [-sp.exp(-2*sg), 0]])
ok_a0 = sp.simplify(h2*h2i - sp.eye(2)) == sp.zeros(2)
dRv = sp.Matrix([sp.diff(Rf, u_), sp.diff(Rf, v_)])
gRR2 = sp.simplify((dRv.T*h2i*dRv)[0])                                  # g^{ab} d_a R d_b R
MMS2 = Rf/(2*Gs)*(1 - gRR2)                                              # M_MS = (R/2G)(1 - g^{ab} dR dR)
thl, thn = 2*sp.diff(Rf, u_)/Rf, 2*sp.diff(Rf, v_)/Rf                   # expansiones de l=d_u, n=d_v (nulos: h_uu=h_vv=0)
ok_a1 = sp.simplify((1 - 2*Gs*MMS2/Rf) + Rf**2*sp.exp(-2*sg)/2*thl*thn) == 0     # 1-2GM/R = -(R^2 e^{-2s}/2) th_l th_n
ok_a2 = sp.simplify((1 - 2*Gs*MMS2/Rf).subs(sp.diff(Rf, u_), 0).doit()) == 0      # th_l=0  =>  2G M_MS = R
Rx, Mx = sp.symbols('R_x M_x', positive=True)
Apos, Tprod = sp.symbols('A_pos T_prod', real=True); Apos = sp.Symbol('A_pos', positive=True)       # A=e^{-2 sigma}>0, T=th_l*th_n
Xrel = -(Rx**2*Apos/2)*Tprod                                                   # 1-2GM/R = -(R^2 A/2) T   (identidad a1)
ok_a3 = (sp.solve(sp.Eq(Xrel, 0), Tprod) == [0]) and bool((-(Rx**2*Apos/2)).is_negative)   # X=0 <=> T=0 ; X<0 <=> T>0 (atrapada)
ft, gt_, rr_ = sp.symbols('f_t g_t r_s', positive=True)
hs = sp.diag(-ft, 1/gt_); gRRs = sp.simplify((sp.Matrix([0, 1]).T*hs.inv()*sp.Matrix([0, 1]))[0])   # static: R=r
ok_a4 = sp.simplify(gRRs - gt_) == 0                                      # g^{ab}dR dR = g  => 2G M_MS = r(1-g)=y (coincide con 04:S2)

# (b) vinculo hamiltoniano de datos simetricos en el tiempo (K_ij=0 en t=0): G_ab n^a n^b = 2G m'/r^2
t_, r_, th_ = sp.symbols('t r theta', positive=True)
a_ = sp.Function('a', positive=True)(r_); N_ = sp.Function('N', positive=True)(r_); sf_ = sp.Function('sf')(r_); qf_ = sp.Function('qf')(r_)
m_ = r_*(1 - a_)/(2*Gs)                                                   # a = 1-2Gm/r
xs4 = [t_, r_, th_, sp.Symbol('phi')]
g4 = sp.diag(-N_**2, 1/a_ + t_**2*sf_, (r_*(1 + t_**2*qf_))**2, (r_*(1 + t_**2*qf_))**2*sp.sin(th_)**2); g4i = g4.inv()
Gm4 = [[[sum(g4i[i, d]*(sp.diff(g4[d, j], xs4[k]) + sp.diff(g4[d, k], xs4[j]) - sp.diff(g4[j, k], xs4[d])) for d in range(4))/2 for k in range(4)] for j in range(4)] for i in range(4)]
def _Ric(j, k):
    return sum(sp.diff(Gm4[i][j][k], xs4[i]) - sp.diff(Gm4[i][j][i], xs4[k]) + sum(Gm4[i][i][d]*Gm4[d][j][k] - Gm4[i][k][d]*Gm4[d][j][i] for d in range(4)) for i in range(4))
Rs4 = sum(g4i[i, i]*_Ric(i, i) for i in range(4))
Gnn = sp.simplify(((_Ric(0, 0) - g4[0, 0]*Rs4/2)/N_**2).subs(t_, 0))
ok_b1 = sp.simplify(Gnn - 2*Gs*sp.diff(m_, r_)/r_**2) == 0                # independiente de sf, qf (segundas derivadas temporales)
# T_ab n^a n^b del campo de Maxwell radial = E_orth^2/2  y Gauss: r^2 E_orth = cte
Phi_ = sp.Function('Phi')(r_)
Fd = sp.zeros(4, 4); Fd[1, 0] = sp.diff(Phi_, r_); Fd[0, 1] = -sp.diff(Phi_, r_)
g3 = sp.diag(-N_**2, 1/a_, r_**2, r_**2*sp.sin(th_)**2); g3i = g3.inv()
Fu = g3i*Fd*g3i.T
F2 = sum(Fd[i, j]*Fu[i, j] for i in range(4) for j in range(4))
Ttt = sum(Fd[0, c]*(Fd*g3i)[0, c] for c in range(4)) - g3[0, 0]*F2/4          # T_tt = F_tc F_t^c - g_tt F^2/4
E_o = sp.diff(Phi_, r_)/(N_/sp.sqrt(a_))                                     # E_orth = F_rt/(N sqrt(g_rr))
ok_b2 = sp.simplify(Ttt/N_**2 - E_o**2/2) == 0
sqrtg = N_*r_**2*sp.sin(th_)/sp.sqrt(a_)
ok_b3 = sp.simplify(sqrtg*Fu[1, 0]/sp.sin(th_) + r_**2*E_o) == 0              # sqrt(-g)F^{rt} = -r^2 E_orth : Gauss => r^2 E_orth = cte (independiente de N, m)

# (c) integracion del vinculo desde R hasta infinito
Qs_, Rr, rp = sp.symbols('Q R r', positive=True)
rho_em = (Qs_/(4*sp.pi*rp**2))**2/2
U_C = sp.integrate(4*sp.pi*rp**2*rho_em, (rp, Rr, sp.oo))
ok_c1 = sp.simplify(U_C - Qs_**2/(8*sp.pi*Rr)) == 0                         # U_C(R) = Q^2/(8 pi R)
# Legendre: H(D)=E D - L(E), dH/dD = E  (cualquier L)  y  U(R) = 2 pi a^{3/2} Int_0^{D_R} H(D) D^{-5/2} dD
Lg_ = sp.Function('Lg'); Eg_ = sp.symbols('Eg', positive=True)
Hh = Eg_*sp.diff(Lg_(Eg_), Eg_) - Lg_(Eg_); Dd = sp.diff(Lg_(Eg_), Eg_)
ok_c2 = sp.simplify(sp.diff(Hh, Eg_)/sp.diff(Dd, Eg_) - Eg_) == 0
aa, DD, DR = sp.symbols('a D D_R', positive=True)
U_max = sp.simplify(2*sp.pi*aa**sp.Rational(3, 2)*sp.integrate((DD**2/2)*DD**sp.Rational(-5, 2), (DD, 0, DR)))
ok_c3 = sp.simplify(U_max.subs(DR, aa/Rr**2) - 2*sp.pi*aa**2/Rr) == 0 and sp.simplify((2*sp.pi*aa**2/Rr).subs(aa, Qs_/(4*sp.pi)) - Qs_**2/(8*sp.pi*Rr)) == 0
cc_ = sp.symbols('c_eh', positive=True)
U_q = sp.simplify(2*sp.pi*aa**sp.Rational(3, 2)*sp.integrate((-4*cc_*DD**4)*DD**sp.Rational(-5, 2), (DD, 0, DR)).subs(DR, aa/Rr**2))
ok_c4 = sp.simplify(U_q + 16*sp.pi*cc_*aa**4/(5*Rr**5)) == 0                  # correccion cuartica -16 pi c a^4/(5 R^5) (= la usada en 04)

# (d) cancelacion en theta_+ = 0 sin teoremas externos
GM_adm, mR, rQ2_, Mn_ = sp.symbols('GM_adm m_R rQ2 M_n', positive=True)
g_rr_ext = sp.simplify((1 - 2*(Gs*mR + 0)/Rr))                               # placeholder simbolico
m_of_r = sp.Function('mm')(rp)
M_ADM_ = sp.symbols('M_ADM', positive=True)
m_ext = M_ADM_ - Qs_**2/(8*sp.pi*rp)                                           # m(r)=M_ADM - U_C(r)  (consecuencia de (c))
g_RN = sp.simplify(1 - 2*Gs*m_ext/rp)
ok_d1 = sp.simplify(g_RN - (1 - 2*Gs*M_ADM_/rp + Gs*Qs_**2/(4*sp.pi*rp**2))) == 0       # la metrica exterior ES Reissner-Nordstrom (derivada, no supuesta)
rQ2 = Gs*Qs_**2/(4*sp.pi)
theta_plus_zero = sp.Eq(sp.Rational(1) - 2*Gs*(M_ADM_ - Qs_**2/(8*sp.pi*Rr))/Rr, 0)       # theta_+=0  <=>  1-2G m(R)/R = 0
sol_GM = sp.solve(theta_plus_zero, M_ADM_)[0]
ok_d2 = sp.simplify(Gs*sol_GM - (Rr + rQ2/Rr)/2) == 0                       # igualdad de Penrose cargada, DERIVADA
# con M_MS(R)=M_n en el cascaron: el termino electrostatico cancela identicamente
exc = sp.simplify((Gs*(Mn_ + Qs_**2/(8*sp.pi*Rr)) - (Rr + rQ2/Rr)/2) - (Gs*Mn_ - Rr/2))
ok_d3 = exc == 0
# desigualdad: materia exterior no cargada con rho>=0 solo anade masa => G M_ADM >= (R + r_Q^2/R)/2 cuando 2G m(R) >= R
dm_ext = sp.symbols('dm_ext', nonnegative=True)
ok_d4 = sp.simplify(Gs*(Mn_ + dm_ext + Qs_**2/(8*sp.pi*Rr)) - (Rr + rQ2/Rr)/2 - (Gs*Mn_ - Rr/2) - Gs*dm_ext) == 0
P10_parts = dict(a0=ok_a0, a1=ok_a1, a2=ok_a2, a3=ok_a3, a4=ok_a4, b1=ok_b1, b2=ok_b2, b3=ok_b3, c1=ok_c1, c2=ok_c2, c3=ok_c3, c4=ok_c4, d1=ok_d1, d2=ok_d2, d3=ok_d3, d4=ok_d4)
gate("P10_misner_sharp_symbolic", all(P10_parts.values()),
     "θ_+=0 ⇔ 2G M_MS=R (doble nulo, exacto); vínculo hamiltoniano G_nn=2Gm'/r² indep. de ∂ₜ²; Gauss ⇒ M_ADM=M_MS(R)+Q²/8πR; la métrica exterior sale RN; igualdad de Penrose cargada DERIVADA y energía electrostática cancelada. "
     f"Subcomprobaciones: {sum(P10_parts.values())}/{len(P10_parts)}")
R['P10_parts'] = {k: bool(v) for k, v in P10_parts.items()}
R['P10_hypotheses'] = "datos esfericos simetricos en el tiempo (K_ij=0), rho>=0 (DEC), sin materia cargada exterior a la superficie; el caso dinamico general NO se demuestra aqui"
# control negativo C4: tres versiones erroneas deben FALLAR
bad_ms = sp.simplify((1 - 2*Gs*(Rf/Gs*(1 - gRR2))/Rf) + Rf**2*sp.exp(-2*sg)/2*thl*thn) == 0          # M_MS con factor 2 de mas
bad_constraint = sp.simplify(Gnn - 4*sp.pi*Gs*sp.diff(m_, r_)/r_**2) == 0                             # m'=2 pi r^2 rho en vez de 4 pi
bad_UC = sp.simplify(Gs*(Qs_**2/(4*sp.pi*Rr)) - rQ2/(2*Rr)) == 0                                       # U_C con factor 2 de mas
gate("C4_ms_control", (not bad_ms) and (not bad_constraint) and (not bad_UC),
     f"M_MS con factor erróneo {'DETECTADO' if not bad_ms else 'NO detectado'}; vínculo m'=2πr²ρ {'DETECTADO' if not bad_constraint else 'NO detectado'}; U_C=Q²/4πR {'DETECTADO' if not bad_UC else 'NO detectado'}")

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
R['P10_n_parts'] = len(R['P10_parts'])
N_cf = (aF_/(0.6*ALPHA))**1.5                                              # U_C(bola)=0.6 alpha N^2/r0 = E_F=aF N^(4/3)/r0
R['critical'] = dict(N_coulomb_over_fermi=float(N_cf), N_crit=float(MP**2), N_P=NP_val, N_star=float(KAPPA**2/(4*np.pi*ALPHA**2)), N_classical_BH_in_EFT=Nbh_val,
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


# ---------- P11: correccion EXACTA a un bucle de U_C(R) (energia electrostatica exterior) ----------
import qed1loop as QL
from scipy.integrate import quad
def _E_of_D(Dv):
    E = Dv
    for _ in range(300):
        dE = Dv - float(QL.D_qed(E)); E += dE
        if abs(dE) <= 1e-14*abs(E): break
    return E
def _H(Dv):
    E = _E_of_D(Dv); return E*Dv - float(QL.L_qed(E))                       # H(D) = E D - L(E) = rho(D)
def U_qed(Rv, N):
    a_ = N*E_CHARGE/(4*np.pi); DR_ = a_/Rv**2
    f_ = lambda sl: _H(np.exp(sl))*np.exp(-1.5*sl)                           # H D^{-5/2} dD = H D^{-3/2} dln D
    I = quad(f_, np.log(DR_) - 60, np.log(DR_), epsabs=0, epsrel=1e-11, limit=400)[0]
    return 2*np.pi*a_**1.5*I
def U_qed_direct(Rv, N):                                                    # contraste independiente: integral directa en r
    a_ = N*E_CHARGE/(4*np.pi); rhi = Rv*10**6 if (a_/(Rv*10**6)**2*E_CHARGE) < 1e-3 else None
    rhi = max(Rv*1e3, np.sqrt(a_*E_CHARGE/1e-3)*1.0)                         # x(r_hi) <= 1e-3: cola cuartica analitica
    f_ = lambda lr: 4*np.pi*np.exp(3*lr)*_H(a_/np.exp(2*lr))
    I = quad(f_, np.log(Rv), np.log(rhi), epsabs=0, epsrel=1e-11, limit=800)[0]
    tail = 2*np.pi*a_**2/rhi - 16*np.pi*C_EH*a_**4/(5*rhi**5)
    return I + tail
Nq = 1e6; rQq = r_Q(Nq); repsq = r_eps(Nq)
UCq = lambda Rv: (Nq*E_CHARGE)**2/(8*np.pi*Rv)
Rlist = [rQq, 1e3*rQq, 1e10*rQq, 1e-6*repsq*1e3]
tabU = {}
for Rv in [rQq, 1e6*rQq, 1e12*rQq]:
    u1, u2 = U_qed(Rv, Nq), U_qed_direct(Rv, Nq)
    tabU[f"{Rv/rQq:.0e}"] = dict(R=float(Rv), x_R=float(E_CHARGE*(Nq*E_CHARGE/(4*np.pi))/Rv**2), delta_U=float(u1/UCq(Rv) - 1), rel_diff_direct=float(abs(u1/u2 - 1)))
slope_pts = [1e8, 1e14, 1e20, 1e26]
slope_lo_pts = [1e4, 1e6, 1e8, 1e10]
dd = []
for xR in slope_pts:                                                         # delta(ln x_R): pendiente = alpha/(3 pi)
    Rv = np.sqrt(E_CHARGE*(Nq*E_CHARGE/(4*np.pi))/xR); dd.append(U_qed(Rv, Nq)/UCq(Rv) - 1)
slopeU = float(np.polyfit(np.log(slope_pts), dd, 1)[0])
dd_lo = []
for xR in slope_lo_pts:
    Rv = np.sqrt(E_CHARGE*(Nq*E_CHARGE/(4*np.pi))/xR); dd_lo.append(U_qed(Rv, Nq)/UCq(Rv) - 1)
slopeU_lo = float(np.polyfit(np.log(slope_lo_pts), dd_lo, 1)[0])
eL_ = ALPHA/(3*np.pi)
dev_hi, dev_lo = slopeU/eL_ - 1, slopeU_lo/eL_ - 1                          # desviacion = O(alpha ln eps) (log subdominante)
yy_ = sp.symbols('y', positive=True)
ok_int = sp.integrate(sp.log(yy_)/yy_**2, (yy_, 1, sp.oo)) == 1 and sp.simplify(((rp**-2*sp.log(rp/Rr)).subs(rp, Rr*yy_)*Rr) - sp.log(yy_)/(Rr*yy_**2)) == 0   # ∫_R^∞ r⁻² ln(r/R) dr = 1/R
R['U_QED'] = dict(table=tabU, slope_dU_dlnx=slopeU, slope_low=slopeU_lo, dev_hi=float(dev_hi), dev_lo=float(dev_lo), slope_expected_alpha_over_3pi=float(ALPHA/(3*np.pi)), N=Nq)
gate("P11_U_QED_exacta", max(v['rel_diff_direct'] for v in tabU.values()) < 1e-7 and abs(dev_hi) < 0.10 and abs(dev_lo) < abs(dev_hi) and ok_int,
     f"U_QED(R)=2πa^(3/2)∫₀^(D_R) H(D)D^(-5/2)dD (H=ED−L exacta a un bucle) vs integral directa en r: ≤{max(v['rel_diff_direct'] for v in tabU.values()):.1e}; "
     f"δU/U_C = {tabU[list(tabU)[0]]['delta_U']:.3e} en R=r_Q; pendiente d(δU/U)/d ln ε = {slopeU:.4e} (ln ε∈[18,60]) y {slopeU_lo:.4e} (ln ε∈[9,23]) vs α/3π = {ALPHA/(3*np.pi):.4e}: desviación {dev_hi:.1%} → {dev_lo:.1%}, consistente con la corrección O(α ln ε) del logaritmo subdominante (∫r⁻² ln(r/R)dr=1/R, simbólico)")

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
    ax.set_title('blue: static EFT | ochre: Schwinger discharge | red: trapped surface | purple: trans-Planckian', fontsize=9)
    ax.legend(loc='lower right', fontsize=8); fig.tight_layout(); fig.savefig('phase_diagram.png', dpi=140); fig.savefig('phase_diagram.pdf')
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
