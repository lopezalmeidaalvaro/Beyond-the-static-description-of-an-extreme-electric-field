"""Lagrangiano de Euler-Heisenberg a un bucle, QED de espinores, campo ELECTRICO puro, tiempo propio EXACTO (m_e=1).
import warnings; warnings.filterwarnings("ignore")
   L1(x) = -(1/8pi^2) Int_0^inf ds/s^3 e^{-s} [ x s cot(x s) - 1 + (x s)^2/3 ],   x = eE = epsilon.
Los polos de cot estan sobre el eje real: Re L1 = parte real del integral sobre el rayo s=rho*exp(-i pi/4);
|Im L1| = mitad de la tasa de decaimiento del vacio (se verifica contra la serie de Schwinger)."""
import numpy as np, mpmath as mp
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from common import ALPHA, E_CHARGE

PHI = np.pi/4
_NB = 28
_COEF = [(-1)**n * 2**(2*n) * float(mp.bernoulli(2*n)) / float(mp.factorial(2*n)) for n in range(_NB)]  # u cot u = sum c_n u^{2n}

def _g_over_u3(u):
    """g(u)/u^3 con g = u cot u - 1 + u^2/3 ; serie para |u|<1 (evita cancelacion)."""
    if abs(u) < 1.0:
        return sum(_COEF[n]*u**(2*n-3) for n in range(2, _NB))
    ex = np.exp(-2j*u)                                  # Im u<0 -> |ex| pequeño (estable)
    cot = 1j*(1 + ex)/(1 - ex)
    return (u*cot - 1 + u*u/3)/u**3

def _I_J(x):
    """I=Int du u^-3 g e^{-u/x},  J=Int du u^-2 g e^{-u/x} sobre el rayo u=rho e^{-i phi}; devuelve (I,J) complejos."""
    ph = np.exp(-1j*PHI)
    def f(t, k):
        rho = np.exp(t); u = rho*ph
        v = _g_over_u3(u)*np.exp(-u/x)*ph*rho
        return (v if k == 0 else v*u)
    tmax = np.log(60*x/np.cos(PHI) + 60)
    out = []
    for k in (0, 1):
        re = quad(lambda t: f(t, k).real, np.log(1e-8), tmax, limit=800, epsabs=0, epsrel=1e-11)[0]
        im = quad(lambda t: f(t, k).imag, np.log(1e-8), tmax, limit=800, epsabs=0, epsrel=1e-11)[0]
        out.append(re + 1j*im)
    return out

def L1_exact(x):
    I, J = _I_J(x)
    return -(x**2/(8*np.pi**2))*I.real, -(x**2/(8*np.pi**2))*I.imag     # (Re L1, Im L1) [signo de Im: ver tests]

def dL1_exact(x):
    I, J = _I_J(x)
    return -(x/(4*np.pi**2))*I.real - J.real/(8*np.pi**2)

def schwinger_series(x, nmax=400):
    n = np.arange(1, nmax+1)
    return x**2/(8*np.pi**3)*np.sum(np.exp(-n*np.pi/x)/n**2)           # Im L1 (espinor: todos los terminos positivos)

# --- tabla + extension asintotica logaritmica ----------------------------------------------
_X = np.logspace(-3, 12, 181)
def _build():
    F1 = np.array([L1_exact(x)[0]/x**2 for x in _X]); F2 = np.array([dL1_exact(x)/x for x in _X])
    return np.log(_X), F1, F2
_lx, _F1, _F2 = _build()
_s1, _s2 = CubicSpline(_lx, _F1), CubicSpline(_lx, _F2)
# I(x) = (1/3) ln x + K + o(1)  ->  K desde x grande
def _K(x): return (_I_J(x)[0].real - np.log(x)/3)
K_ASYM = _K(1e12)
K_DRIFT = abs(_K(1e12) - _K(1e10))
def L1_fast(x):
    x = np.asarray(x, float); lx = np.log(x)
    asym = -(1/(8*np.pi**2))*(lx/3 + K_ASYM)
    return x**2*np.where(lx > _lx[-1], asym, _s1(np.minimum(lx, _lx[-1])))
def dL1_fast(x):
    x = np.asarray(x, float); lx = np.log(x)
    asym = -(1/(12*np.pi**2))*(lx + 3*K_ASYM + 0.5)
    return x*np.where(lx > _lx[-1], asym, _s2(np.minimum(lx, _lx[-1])))
# funciones en E (m=1)
L_qed = lambda E: E**2/2 + L1_fast(E_CHARGE*np.abs(E))
D_qed = lambda E: E + E_CHARGE*dL1_fast(E_CHARGE*np.abs(E))
def pair_rate(x):
    """TASA MEDIA DE PARES por volumen y tiempo (por especie): solo el termino n=1 (f=exp(-pi m_perp^2/x))."""
    return x**2/(4*np.pi**3)*np.exp(-np.pi/x)
def vacuum_decay_rate(x):
    """2 Im L1 = Gamma_vac = x^2/(4 pi^3) sum_n n^-2 exp(-n pi/x)  (era la 'Gamma' de eeh.py; NO es el numero medio de pares)."""
    n = np.arange(1, 600)[:, None]
    return (np.asarray(x)**2/(4*np.pi**3)*np.sum(np.exp(-n*np.pi/np.asarray(x))/n**2, axis=0))
