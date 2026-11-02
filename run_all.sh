#!/usr/bin/env bash
set -eo pipefail
cd "$(dirname "$0")"
rm -f discovery_results.json
python3 01_derive.py > out_01_derive.txt
python3 02_invariants.py > out_02_invariants.txt
python3 03_integrate.py > out_03_integrate.txt
python3 04_nonperturbative_qed.py > out_04_nonperturbative_qed.txt
python3 05_dynamical_backreaction.py > out_05_dynamical_backreaction.txt
python3 05b_plasma_oscillations.py > out_05b_plasma_oscillations.txt
python3 05c_quantum_vlasov_exact.py > out_05c_quantum_vlasov_exact.txt
python3 06_phase_diagram.py > out_06_phase_diagram.txt
python3 06b_dynamical_trapping_1d.py > out_06b_dynamical_trapping_1d.txt
python3 07_make_report.py | tee out_07_report.txt
python3 08_novelty_checks.py > out_08_novelty_checks.txt
python3 09_make_audit.py | tee out_09_audit.txt
python3 10_build_paper.py | tee out_10_paper.txt
echo "OK: baseline + ramas A/B/C verificadas. Ver out_*.txt, results.json, discovery_results.json, phase_diagram.png, DISCOVERY_REPORT.md, NOVELTY_AUDIT.md, novelty_audit.json, paper.tex, references.bib, paper.pdf"
