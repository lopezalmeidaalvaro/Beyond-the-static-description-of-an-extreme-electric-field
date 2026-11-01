"""Constantes en unidades naturales hbar=c=1, Heaviside-Lorentz (eps0=1), m_e=1.
Unidad de longitud = longitud de Compton reducida lambda_C = hbar/(m_e c).
Unidad de carga: Q adimensional;  Q = N*e  con e = sqrt(4*pi*alpha)."""
import math

ALPHA = 1 / 137.035999084
E_CHARGE = math.sqrt(4 * math.pi * ALPHA)       # carga elemental (HL)
C_EH = ALPHA**2 / 90.0                          # alpha^2/(90 m^4), m_e=1
MPL = 1.220890e19 / 0.51099895e-3               # masa de Planck (no reducida) en unidades de m_e
G_NEWTON = 1.0 / MPL**2                         # G = 1/M_Pl^2

# ---------------- persistencia de resultados (discovery_results.json) ----------------
import json, os
_DR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'discovery_results.json')
def save_section(name, data):
    db = json.load(open(_DR)) if os.path.exists(_DR) else {}
    db[name] = data
    json.dump(db, open(_DR, 'w'), indent=1, default=float)
def load_all():
    return json.load(open(_DR)) if os.path.exists(_DR) else {}
