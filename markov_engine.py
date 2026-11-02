"""Motor de Markov con bloqueo de Pauli (copia literal de las funciones make/run/analyze de 05b_plasma_oscillations.py, sin ejecucion al importar).
Se usa en 05c para comparar Markov vs cinetica exacta con los MISMOS epsilon_0. Equivalencia con 05b verificada por el gate M0 de 05c."""
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from common import *
e2 = E_CHARGE**2

def make(x0, sigf=0.25, dpf=1/3, J=10, Pfac=52.0, block=True, polblock=True, flip=False, wrongsign=False):
    sc = np.sqrt(x0); sig = sigf*sc; dp = sig*dpf
    Pmax = Pfac*sc + 5*sig; P = np.arange(-Pmax, Pmax + dp/2, dp); NP = len(P)
    pm = 4*np.sqrt(x0/np.pi)
    xg, wg = np.polynomial.legendre.leggauss(J); pp = 0.5*pm*(xg + 1); wp = 0.5*pm*wg
    W = (dp/(2*np.pi))*(2*pp*wp/(2*np.pi))[None, :]*np.ones((NP, 1))              # medida (incluye espin 2)
    m2 = 1 + pp**2
    def rhs(t, y):
        x, Phi = y[0], y[1]; f = y[2:].reshape(NP, J)
        p = (P + Phi)[:, None]; om = np.sqrt(p*p + m2[None, :]); ax = abs(x)
        A = (np.exp(-np.pi*m2/ax)[None, :]*np.exp(-0.5*(p/sig)**2)/(sig*np.sqrt(2*np.pi))) if ax > 0 else np.zeros((NP, J))
        S = ax*A
        fac = (1 + 2*f) if wrongsign else ((1 - 2*f) if block else 1.0)
        fdot = fac*S
        polfac = (1 + 2*f) if wrongsign else ((1 - 2*f) if (block and polblock) else 1.0)
        cond = np.sum(W*f*p/om); pol = np.sign(x)*np.sum(W*om*polfac*A)
        return np.concatenate(([-2*e2*((-1 if flip else 1)*cond + pol), x], fdot.ravel()))
    def energy(y):
        x, Phi = y[0], y[1]; f = y[2:].reshape(NP, J); p = (P + Phi)[:, None]; om = np.sqrt(p*p + m2[None, :])
        return x**2/(2*e2) + 2*np.sum(W*f*om)
    def npairs(y): return float(np.sum(W*y[2:].reshape(NP, J)))
    return dict(rhs=rhs, energy=energy, npairs=npairs, NP=NP, J=J, P=P, sig=sig, sc=sc, x0=x0)

def run(x0, ncross=10, rtol=1e-9, tmax_factor=60, **kw):
    K = make(x0, **kw); NP, J = K['NP'], K['J']
    y = np.zeros(2 + NP*J); y[0] = x0; E0 = K['energy'](y)
    tg = 13.5/(E_CHARGE*np.sqrt(x0)); t = 0.0; cr = []; Ed = []; fmax = 0.0; fmin = 0.0; Pcov = 0.0
    atol = np.concatenate(([1e-9*x0, 1e-9*K['sc']*30], np.full(NP*J, 1e-12)))
    T = []; X = []; Np = []; tstart = time.time()
    while len(cr) < ncross and t < tmax_factor*tg:
        sol = solve_ivp(K['rhs'], [t, t + tg], y, method='DOP853', rtol=rtol, atol=atol, dense_output=True)
        te = np.linspace(t, t + tg, 160); ys = sol.sol(te); xs = ys[0]
        for i in range(len(te) - 1):
            if xs[i]*xs[i+1] < 0:
                tc = brentq(lambda tt: sol.sol(tt)[0], te[i], te[i+1], xtol=1e-12*tg); cr.append(tc)
                Ed.append(abs(K['energy'](sol.sol(tc)) - E0)/E0)
        T += list(te); X += list(xs); Np += [K['npairs'](ys[:, i]) for i in range(0, len(te), 20)]
        y = sol.y[:, -1]; t = sol.t[-1]; ff = y[2:]; fmax = max(fmax, float(ff.max())); fmin = min(fmin, float(ff.min()))
        Pcov = max(Pcov, abs(y[1]))
        if not np.isfinite(y).all(): break
    ff = y[2:].reshape(NP, J); edge = float(max(abs(ff[:5]).max(), abs(ff[-5:]).max()))
    T = np.array(T); X = np.array(X); cr = np.array(cr)
    amps = [float(np.abs(X[(T > cr[i]) & (T < cr[i+1])]).max()/x0) for i in range(len(cr) - 1)] if len(cr) > 1 else []
    return dict(cr=cr, Ed=Ed, T=T, X=X, amps=amps, fmax=fmax, fmin=fmin, edge_f=edge, Phimax=float(Pcov), Pmax=float(K['P'][-1]), NP=NP, E0=E0,
                npairs_end=K['npairs'](y), runtime=time.time() - tstart, x0=x0)

def analyze(r):
    cr = r['cr']; hp = np.diff(np.concatenate(([0.0], cr)))                      # hp[0]=tiempo al 1er cero; luego semiperiodos
    half = np.diff(cr); w = np.pi/half if len(half) else np.array([])
    out = dict(n_cross=int(len(cr)), tau0=float(cr[0]) if len(cr) else None, half_periods=[float(v) for v in half], amps=r['amps'],
               omega_early=float(w[0]) if len(w) else None, omega_late=float(np.mean(w[-3:])) if len(w) >= 3 else None, Emax_err=float(max(r['Ed'])) if r['Ed'] else None,
               fmax=r['fmax'], fmin=r['fmin'], edge_f=r['edge_f'], Phimax_over_Pmax=r['Phimax']/r['Pmax'], npairs_end=r['npairs_end'], runtime=r['runtime'])
    if len(r['amps']) >= 3:
        tm = 0.5*(cr[:-1] + cr[1:]); out['gamma'] = float(-np.polyfit(tm, np.log(r['amps']), 1)[0])
    return out
