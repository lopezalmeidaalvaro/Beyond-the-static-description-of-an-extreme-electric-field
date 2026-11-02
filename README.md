# Experimento: campo eléctrico extremo más allá de la EFT estática (Einstein–Maxwell–Euler–Heisenberg)

Pipeline "el LLM propone, el CAS ejecuta": SymPy deriva y verifica, una suite de invariantes actúa de compuerta
y los solvers numéricos se detienen donde la EFT deja de ser válida.

```bash
pip install -r requirements.txt
# 05c compila qvlasov_kernel.c con gcc -fopenmp (necesita gcc); el barrido completo tarda ~20 min en 4 núcleos
# para compilar paper.pdf (Debian/Ubuntu):
apt-get install -y texlive-latex-base texlive-latex-recommended texlive-latex-extra texlive-publishers texlive-science texlive-fonts-recommended texlive-bibtex-extra
./run_all.sh        # ejecuta 01..10; falla si algún gate o control negativo falla
```

| Script | Contenido |
|---|---|
| `01_derive.py`, `02_invariants.py`, `03_integrate.py` | Línea base: EEH cuártico estático, 5 invariantes, integración con cortes de EFT |
| `04_nonperturbative_qed.py` | QED a un bucle (tiempo propio exacto), perfil interior vs cuártico / Born–Infeld |
| `05_dynamical_backreaction.py` | Descarga de Schwinger con retroacción (LCFA, válida hasta el primer cero del campo) |
| `06_phase_diagram.py` | Penrose cargada, leyes de escala y diagrama de fases (`phase_diagram.png`) |
| `08_novelty_checks.py`, `09_make_audit.py` | Comprobaciones numéricas de la auditoría y generación de `NOVELTY_AUDIT.md` / `novelty_audit.json` |
| `05b_plasma_oscillations.py` | Cinética de Markov con bloqueo de Pauli a través de ≥5 oscilaciones de plasma (extiende 05) |
| `05c_quantum_vlasov_exact.py` | Cinética cuántica exacta (Bloch/Bogoliubov + Maxwell): gates Q0–Q6, Landau–Zener, barrido ε₀∈[4,10³], comparación con Markov, ecuación de estado (motor: `qvlasov_engine.py` + `qvlasov_kernel.c`, compilado con gcc/OpenMP; `markov_engine.py` = copia de 05b) |
| `10_build_paper.py` | Genera `paper.tex`/`references.bib` (REVTeX 4-2, PRD) con cifras leídas de `discovery_results.json` y compila `paper.pdf` |
| `06b_dynamical_trapping_1d.py` | Misner–Sharp dinámico con K_ij≠0 (SymPy: D_MS1–D_MS5, P12, controles) + solver 1+1D polar-areal (`dyn1d.py`) de la descarga de una bola cargada: balance, convergencia, canales fuga/inflow, desplazamiento del umbral y mapa de fases corregido |
| `07_make_report.py` | Compuerta global (67 gates, 15 controles negativos) y genera `DISCOVERY_REPORT.md` |

Ver `DISCOVERY_REPORT.md` (resultados y límites), `SYSTEM_PROMPT.md` (reglas del agente) y `discovery_results.json` (todas las métricas).
