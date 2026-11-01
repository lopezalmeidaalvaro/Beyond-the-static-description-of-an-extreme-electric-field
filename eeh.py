"""Biblioteca numerica del sistema reducido (deducido y verificado en 01_derive.py)."""
import numpy as np
from scipy.special import spence
from common import *

def E_of_r(r, Q, c=C_EH):
    """Raiz unica de r^2 E (1+16 c E^2) = Q/4pi  (Newton monotono; h(E) convexa creciente)."""
    E0 = Q / (4*np.pi*np.asarray(r, float)**2)
    E = E0.copy() if isinstance(E0, np.ndarray) else E0
    for _ in range(200):
        h = E*(1 + 16*c*E**2) - E0
        dE = h / (1 + 48*c*E**2)
        E = E - dE
        if np.all(np.abs(dE) <= 1e-15*np.abs(E)):
            break
    return E

rho    = lambda E, c=C_EH: E**2/2 + 12*c*E**4
p_perp = lambda E, c=C_EH: E**2/2 + 4*c*E**4
eps_of = lambda E: E_CHARGE*np.abs(E)                 # eps = e E/m^2  (m_e=1)

def schwinger_rate(E):
    """Gamma/m^4 = (eE)^2/(4 pi^3) sum_n n^-2 exp(-n pi/eps) = eps^2/(4pi^3) Li2(exp(-pi/eps))."""
    e = eps_of(E)
    return e**2/(4*np.pi**3) * spence(1 - np.exp(-np.pi/e))   # spence(1-z)=Li2(z)

def ell_field(r, E, c=C_EH):
    """Escala de variacion |E/E'| (E' de la integral primera) -> criterio de expansion derivativa."""
    return r*(1 + 48*c*E**2)/(2*(1 + 16*c*E**2))
