"""SOLVER NUMERICO: integra desde r>>lambda_C hacia el centro el sistema reducido
   y(r)=2Gm(r),  g=f=1-y/r,  Phi'=-E,  r^2 E(1+16cE^2)=Q/4pi
con compuertas de parada de la EFT.  NUNCA extrapola mas alla del primer criterio violado."""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp, quad
from eeh import *

EPS_WARN, EPS_MAX = 0.1, 1.0        # parametro de pequeñez eps=eE/m^2: aviso / ruptura (polo de Schwinger ~ eps~1)
ELL_MIN = 1.0                       # escala de variacion del campo >= lambda_C (expansion derivativa)
RHO_PL, R_PL = MPL**4, 1.0/MPL      # densidad de energia de Planck, longitud de Planck

def run(N, kappa=None, M_over_m=None, c=C_EH, label="", verbose=True):
    Q = N*E_CHARGE; a = Q/(4*np.pi)
    r_eps1 = np.sqrt(N*ALPHA)                       # E(r)=Q/4pi r^2 -> eps=1 en r=sqrt(N alpha) (Coulomb)
    rQ = np.sqrt(G_NEWTON/(4*np.pi))*Q              # radio de RN
    GM = kappa*rQ if kappa is not None else G_NEWTON*M_over_m
    r_max = 1e3*max(r_eps1, 1.0)
    # condicion de contorno asintotica (cola EH incluida a primer orden)
    y0 = 2*GM - 8*np.pi*G_NEWTON*(a**2/(2*r_max) - 4*c*a**4/(5*r_max**5))
    P0 = a/r_max - 16*c*a**3/(7*r_max**7)

    def rhs(s, u):                                   # s = ln r (descendente)
        r = np.exp(s); E = E_of_r(r, Q, c)
        return [r*8*np.pi*G_NEWTON*r**2*rho(E, c), -r*(-E)]   # dy/ds, dPhi/ds  (Phi'=-E)
    def ev(fun, name):
        def h(s, u): return fun(np.exp(s), E_of_r(np.exp(s), Q, c), u)
        h.terminal = True; h.direction = -1; h.name = name; return h
    crit = [ev(lambda r,E,u: np.log(EPS_MAX/eps_of(E)), "eps=eE/m^2 >= 1 (serie EH / Schwinger no perturbativo)"),
            ev(lambda r,E,u: np.log(ell_field(r,E,c)/ELL_MIN), "escala de variacion del campo < lambda_C (expansion derivativa)"),
            ev(lambda r,E,u: np.log(RHO_PL/rho(E,c)), "densidad de energia >= M_Pl^4 (gravedad cuantica)"),
            ev(lambda r,E,u: np.log(r/R_PL/10), "r < 10 l_Planck"),
            ev(lambda r,E,u: -np.log(max(schwinger_rate(E), 1e-300)), "tasa de creacion de pares Gamma*lambda_C^4 >= 1 (unitariedad)")]
    hor = lambda s, u: 1 - u[0]/np.exp(s); hor.terminal = False
    warn = lambda s, u: np.log(EPS_WARN/eps_of(E_of_r(np.exp(s), Q, c))); warn.terminal = False
    sol = solve_ivp(rhs, [np.log(r_max), np.log(R_PL)], [y0, P0], events=crit + [hor, warn],
                    method='DOP853', rtol=1e-12, atol=0.0 if y0 else 1e-300, dense_output=True)
    stops = [(sol.t_events[i][0], crit[i].name) for i in range(len(crit)) if len(sol.t_events[i])]
    s_stop, why = (min(stops, key=lambda z: -z[0]) if stops else (sol.t[-1], "ninguno (alcanzado r_min)"))
    r_stop = np.exp(s_stop); E_stop = E_of_r(r_stop, Q, c)
    horiz = sorted(np.exp(sol.t_events[len(crit)]).tolist(), reverse=True)
    r_warn = np.exp(sol.t_events[len(crit)+1][0]) if len(sol.t_events[len(crit)+1]) else None

    # --- verificaciones numericas del propio solver -------------------------------------
    rq = np.exp(np.linspace(np.log(r_max), s_stop, 7))            # cuadratura independiente
    yq = lambda r: y0 + 8*np.pi*G_NEWTON*quad(lambda z: np.exp(3*z)*rho(E_of_r(np.exp(z),Q,c),c), np.log(r_max), np.log(r),
                                               epsabs=0, epsrel=1e-13, limit=500)[0]
    ydiff = np.max([abs(sol.sol(np.log(r))[0]-yq(r))/abs(yq(r)) for r in rq[1:]])
    # g EH vs RN exacto (misma M,Q): y_RN = 2GM - G Q^2/(4 pi r)
    y_RN = lambda r: 2*GM - G_NEWTON*Q**2/(4*np.pi*r)
    gEH = 1 - sol.sol(s_stop)[0]/r_stop; gRN = 1 - y_RN(r_stop)/r_stop
    out = dict(label=label, N=N, r_eps1=r_eps1, rQ=rQ, r_stop=r_stop, why=why, eps_stop=float(eps_of(E_stop)),
               rho_over_Pl=float(rho(E_stop,c)/RHO_PL), horizons=horiz, r_warn=r_warn,
               dEH_E=float(16*c*E_stop**2), g_EH=float(gEH), g_RN=float(gRN), ode_vs_quad=float(ydiff))
    if verbose:
        print(f"\n--- {label}  (N={N:.3g}, Q=N e) ---")
        print(f"  aviso eps>0.1 en r = {r_warn if r_warn is None else f'{r_warn:.4g}'} lambda_C ;  radio RN r_Q = {rQ:.3g} lambda_C ; "
              f"eps=1 (Coulomb) en r = {r_eps1:.4g} lambda_C")
        print(f"  PARADA en r = {r_stop:.5g} lambda_C   [{why}]")
        print(f"  eps(r_stop)={out['eps_stop']:.3g}  rho/M_Pl^4={out['rho_over_Pl']:.2e}  "
              f"correccion EH relativa 16cE^2={out['dEH_E']:.2e}")
        print(f"  horizontes (g=0) encontrados antes de la parada: {horiz if horiz else 'ninguno'}")
        print(f"  g_EH(r_stop)={gEH:.12g}   g_RN(r_stop)={gRN:.12g}   |g_EH-g_RN|={abs(gEH-gRN):.2e}")
        print(f"  ODE vs cuadratura (y): {ydiff:.1e}")
        print(f"  DIAGNOSTICO: ruptura de la EFT en r={r_stop:.4g} lambda_C -> NO se extrapola mas adentro "
              f"(no se afirma singularidad ni geometria interior)")
    return out

if __name__ == "__main__":
    res = []
    # (0) control: c=0 debe reproducir Reissner-Nordstrom exactamente (limite hbar->0)
    kappa = 1.5
    r0 = run(1e48, kappa=kappa, c=0.0, label="CONTROL c=0  (RN exacto)")
    rQ = r0["rQ"]; GM = kappa*rQ
    gexact = 1 - 2*GM/r0["r_stop"] + rQ**2/r0["r_stop"]**2
    print(f"  RN analitico g={gexact:.12g}; |g_num-g_RN,exacto| = {abs(r0['g_EH']-gexact):.2e}")
    assert abs(r0['g_EH']-gexact) < 1e-9, "El control c=0 no reproduce RN: solver roto"
    res.append(r0)
    res.append(run(1,   M_over_m=1.0, label="A: electron (Q=e, M=m_e)"))
    res.append(run(1e6, M_over_m=1e6, label="B: cuerpo cargado N=1e6"))
    res.append(run(1e20, M_over_m=1e20, label="C: N=1e20"))
    res.append(run(1e48, kappa=kappa, label="D: RN macroscopico, N=1e48, GM=1.5 r_Q (horizontes en regimen valido)"))
    json.dump(res, open('results.json', 'w'), indent=1, default=float)
