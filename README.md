# Experimento: campo eléctrico extremo más allá de la EFT estática (Einstein–Maxwell–Euler–Heisenberg)

Pipeline "el LLM propone, el CAS ejecuta": SymPy deriva y verifica, una suite de invariantes actúa de compuerta
y los solvers numéricos se detienen donde la EFT deja de ser válida.

```bash
pip install -r requirements.txt
./run_all.sh        # ejecuta 01..07; falla si algún gate o control negativo falla
```

| Script | Contenido |
|---|---|
| `01_derive.py`, `02_invariants.py`, `03_integrate.py` | Línea base: EEH cuártico estático, 5 invariantes, integración con cortes de EFT |
| `04_nonperturbative_qed.py` | QED a un bucle (tiempo propio exacto), perfil interior vs cuártico / Born–Infeld |
| `05_dynamical_backreaction.py` | Descarga de Schwinger con retroacción (LCFA, válida hasta el primer cero del campo) |
| `06_phase_diagram.py` | Penrose cargada, leyes de escala y diagrama de fases (`phase_diagram.png`) |
| `07_make_report.py` | Compuerta global (50 gates, 10 controles negativos) y genera `DISCOVERY_REPORT.md` |

Ver `DISCOVERY_REPORT.md` (resultados y límites), `SYSTEM_PROMPT.md` (reglas del agente) y `discovery_results.json` (todas las métricas).
