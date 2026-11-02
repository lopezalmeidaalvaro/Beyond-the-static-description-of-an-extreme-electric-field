"""Compuerta global + generador de DISCOVERY_REPORT.md.
 - exige 100% de gates (simbolicos, numericos, controles negativos) en 04/05/06 y que el baseline siga vigente;
 - renderiza el informe SOLO con valores leidos de discovery_results.json / results.json (KeyError => falla)."""
import json, re, sys
D = json.load(open('discovery_results.json')); BASE = json.load(open('results.json'))
allg = {}
for br in ('branch_A', 'branch_B', 'branch_B2', 'branch_B3', 'branch_C', 'branch_C2'):
    for k, v in D[br]['gates'].items(): allg[f"{br}:{k}"] = v
bad = [k for k, v in allg.items() if not v]
neg = [k for k in allg if re.match(r'branch_(A|B|B2|B3|C|C2):C', k)]
if bad: print("GATES FALLIDOS:", bad); sys.exit(1)
if len(neg) < 24: print("faltan controles negativos"); sys.exit(1)
D['meta'] = dict(n_gates=len(allg), n_negative_controls=len(neg), all_pass=True)
json.dump(D, open('discovery_results.json', 'w'), indent=1)

def get(path):
    if path.startswith('base.'):
        cur = BASE[4]; path = path[5:]
    else: cur = D
    parts = path.split('.'); i = 0
    while i < len(parts):
        if isinstance(cur, list): cur = cur[int(parts[i])]; i += 1; continue
        for j in range(len(parts), i, -1):                         # coincidencia mas larga (claves con puntos, p.ej. "1.0")
            key = '.'.join(parts[i:j])
            if key in cur: cur = cur[key]; i = j; break
        else: raise KeyError(f"{'.'.join(parts[:i])} -> {'.'.join(parts[i:])}")
    return cur
def sub(m):
    path, fmt = (m.group(1).split('|') + [''])[:2]
    v = get(path)
    return format(v, fmt) if fmt else str(v)
A, B, C = D['branch_A'], D['branch_B'], D['branch_C']
rt = C['race_at_rs']; N_race = rt[1]['N'] * rt[1]['tau_dis_over_ff']**2     # tau_dis/tau_ff ∝ N^(-1/2) en el regimen de masa en reposo
D['branch_C']['critical']['N_race_at_trap_boundary'] = N_race
json.dump(D, open('discovery_results.json', 'w'), indent=1)

tab = "| N | r_ε=√(Nα) | r_Q (RN) | r_s (atrapamiento) | ε(r_s) | τ_dis/r0 |\n|---|---|---|---|---|---|\n"
for k, v in C['table'].items():
    tab += f"| {k} | {v['r_eps']:.3g} | {v['r_Q']:.3g} | {v['r_s']:.3g} | {v['eps_at_rs']:.2g} | {v['tau_dis_over_r0']:.3g} |\n"
SC = {'branch_A': '04_nonperturbative_qed.py', 'branch_B': '05_dynamical_backreaction.py', 'branch_B2': '05b_plasma_oscillations.py', 'branch_B3': '05c_quantum_vlasov_exact.py', 'branch_C': '06_phase_diagram.py', 'branch_C2': '06b_dynamical_trapping_1d.py'}
trace = "| Gate | Script | Estado |\n|---|---|---|\n" + "".join(f"| {k.split(':')[1]} | {SC[k.split(':')[0]]} | PASS |\n" for k in allg)

B3 = D['branch_B3']
def _f(v, f): return '—' if v is None else format(v, f)
tabx = "| ε₀ | f_max exacto | f_max Markov | τ₀·e√ε₀ exacto | Markov | ω_p/(e√ε₀) exacto | Markov | γ/ω_p exacto | Markov | ⟨w∥⟩_pl | cruces |\n|---|---|---|---|---|---|---|---|---|---|---|\n"
for k, a in B3['exact_sweep'].items():
    m = B3['markov_sweep'][k]
    tabx += f"| {float(k):g} | {a['fmax']:.5f} | {m['fmax']:.5f} | {a['tau0_e_sqrtx0']:.3f} | {_f(m.get('tau0_e_sqrtx0'), '.3f')} | {_f(a['omega_over_esqrtx0'], '.4f')} | {_f(m.get('omega_over_esqrtx0'), '.4f')} | {_f(a.get('gamma_over_omega'), '.4f')} | {_f(m.get('gamma_over_omega'), '.4f')} | {_f(a.get('w_pl_avg_last_period'), '.3f')} | {a['n_cross']} |\n"

T = open('report_template.md').read()
T = T.replace('@@TABLA_ESCALAS@@', tab).replace('@@TRAZABILIDAD@@', trace).replace('@@TABLA_EXACTA@@', tabx).replace('@@N_RACE@@', f"{N_race:.3g}")
out = re.sub(r'⟦([^⟧]+)⟧', sub, T)
open('DISCOVERY_REPORT.md', 'w').write(out)
print(f"OK: {len(allg)} gates PASS ({len(neg)} controles negativos). Informe: DISCOVERY_REPORT.md ({len(out)} caracteres)")
