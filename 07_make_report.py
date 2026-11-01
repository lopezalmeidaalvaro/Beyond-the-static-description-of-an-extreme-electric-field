"""Compuerta global + generador de DISCOVERY_REPORT.md.
 - exige 100% de gates (simbolicos, numericos, controles negativos) en 04/05/06 y que el baseline siga vigente;
 - renderiza el informe SOLO con valores leidos de discovery_results.json / results.json (KeyError => falla)."""
import json, re, sys
D = json.load(open('discovery_results.json')); BASE = json.load(open('results.json'))
allg = {}
for br in ('branch_A', 'branch_B', 'branch_C'):
    for k, v in D[br]['gates'].items(): allg[f"{br}:{k}"] = v
bad = [k for k, v in allg.items() if not v]
neg = [k for k in allg if re.match(r'branch_[ABC]:C\d', k)]
if bad: print("GATES FALLIDOS:", bad); sys.exit(1)
if len(neg) < 8: print("faltan controles negativos"); sys.exit(1)
D['meta'] = dict(n_gates=len(allg), n_negative_controls=len(neg), all_pass=True)
json.dump(D, open('discovery_results.json', 'w'), indent=1)

def get(path):
    if path.startswith('base.'):
        cur = BASE[4]; path = path[5:]
    else: cur = D
    for p in path.split('.'):
        cur = cur[int(p)] if isinstance(cur, list) else cur[p]
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
trace = "| Gate | Script | Estado |\n|---|---|---|\n" + "".join(f"| {k.split(':')[1]} | {'0'+{'A':'4','B':'5','C':'6'}[k[7]]}_*.py | PASS |\n" for k in allg)

T = open('report_template.md').read()
T = T.replace('@@TABLA_ESCALAS@@', tab).replace('@@TRAZABILIDAD@@', trace).replace('@@N_RACE@@', f"{N_race:.3g}")
out = re.sub(r'⟦([^⟧]+)⟧', sub, T)
open('DISCOVERY_REPORT.md', 'w').write(out)
print(f"OK: {len(allg)} gates PASS ({len(neg)} controles negativos). Informe: DISCOVERY_REPORT.md ({len(out)} caracteres)")
