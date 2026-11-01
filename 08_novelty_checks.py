"""Comprobaciones numericas que alimentan NOVELTY_AUDIT.md (todo reproducible; nada de memoria).
 a) polo de Landau en GeV vs el valor citado en la literatura (~1e277 GeV, extracto de arXiv:1812.08673)
 b) umbral de Ritus-Narozhny: alpha*chi^(2/3)=1  ->  chi=alpha^(-3/2); para una particula en reposo chi=eps
 c) validez de la expansion derivativa en los perfiles de 04 (escala de variacion del campo vs lambda_C)
 d) contraste cuantitativo del tiempo de descarga con el dato de Ruffini et al. (extracto: L~10 hbar/mc a E0=10 Ec; ~1e7 a 0.15 Ec)
 e) DEC: umbral exacto y asintota (anomalia de traza)"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from common import *
D = load_all(); A = D['branch_A']; R = {}
M_E_GEV = 0.51099895e-3
# a)
lnL = A['Landau']['ln_eps_L']; log10_mu = np.log10(M_E_GEV) + 0.5*lnL/np.log(10)
R['landau_scale_GeV_log10'] = float(log10_mu); R['landau_literature_GeV_log10'] = 277.0
# b)
R['ritus_narozhny_chi_threshold'] = float(ALPHA**-1.5); R['ritus_narozhny_eps_at_rest'] = float(ALPHA**-1.5)
R['eps_sweep_max_in_05'] = 1e10; R['orders_above_RN_threshold_in_05'] = float(np.log10(1e10/ALPHA**-1.5))
R['eps_planck_over_RN'] = float(10**A['Landau']['log10_eps_Planck']/ALPHA**-1.5)
# c)  ell = |E/E'| ~ r/2 (Coulomb);  EFT derivativa requiere ell >= lambda_C (criterio del baseline en 03)
for key, N in (('1e+06', 1e6), ('1e+20', 1e20)):
    rQ = A['profiles'][key]['RN_analitico']['rQ']
    R[f'ell_over_lambdaC_at_rQ_N{key}'] = float(rQ/2)
R['ell_criterion_violated_in_04_profiles'] = bool(R['ell_over_lambdaC_at_rQ_N1e+06'] < 1.0)
# d)
src = open('05_dynamical_backreaction.py').read()
pre = src[:src.index("# --- convergencia")].replace('gate("D1', '(lambda *a,**k: None)("D1')
body = pre.split("# ============================== SIMBOLICO")[0] + "\nJ_NODES=8\n" + pre[pre.index("sj, wj = laggauss"):]
import io, contextlib
ns = {}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(body, 'x', 'exec'), ns)
run = ns['run']
cmp_ = {}
for x0, Tf, Nts in ((10.0, 40.0, (1500, 3000)), (0.15, 3.5e6, (3000, 6000))):
    t = [run(x0, Nt, Tf=Tf)['tau0'] for Nt in Nts]
    cmp_[f"{x0:g}"] = dict(tau0_first_zero=t, rel_change_with_Nt=float(abs(t[1]/t[0]-1)))
R['ruffini_crosscheck'] = dict(ours=cmp_, literature_extract={'10': 10.0, '0.15': 1e7}, caveat='Definicion de "longitud de oscilacion" no verificada (solo extracto); es consistencia de orden de magnitud, NO validacion')
R['ruffini_ratio_ours'] = float(cmp_['0.15']['tau0_first_zero'][1]/cmp_['10']['tau0_first_zero'][1]); R['ruffini_ratio_literature'] = 1e7/10.0
# e)
R['dec_onset_exact'] = A['energy_conditions_qed']['DEC_onset_eps_exact']
R['dec_margin_asymptote'] = A['energy_conditions_qed']['DEC_margin_asymptote_alpha_over_3pi']
R['trace_anomaly_relation'] = 'ρ−|L| = x³ d(L1/x²)/dx / (...) ; asintota −α/(3π) = coeficiente β de un bucle'
save_section('audit_checks', R)
for k, v in R.items(): print(f"{k}: {v}")
