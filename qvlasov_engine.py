"""Motor numerico del sistema exacto de Bloch/Bogoliubov (nucleo C via ctypes + malla + analisis). Importado por 05c_quantum_vlasov_exact.py."""
import os, time, ctypes, subprocess
import numpy as np
from common import *
e2 = E_CHARGE**2

# ================================ 2. COMPILACION DEL NUCLEO C ================================
SRC, LIB = 'qvlasov_kernel.c', 'libqvlasov.so'
if not os.path.exists(LIB) or os.path.getmtime(LIB) < os.path.getmtime(SRC):
    subprocess.check_call(['gcc', '-O3', '-ffast-math', '-fopenmp', '-shared', '-fPIC', '-o', LIB, SRC, '-lm'])
lib = ctypes.CDLL(os.path.abspath(LIB)); dp_ = ctypes.POINTER(ctypes.c_double)
lib.integrate.restype = ctypes.c_long
lib.integrate.argtypes = [ctypes.c_long] + [dp_]*3 + [ctypes.c_double]*3 + [dp_]*4 + [ctypes.c_double, ctypes.c_long, ctypes.c_long, dp_, ctypes.c_long, ctypes.c_double]
def ptr(a): return a.ctypes.data_as(dp_)
def kernel(P, s, W, state, f, u, v, dt, ns, every, clam=1.0, cpol=1.0, e2_=e2, pur_ommax=1e300):
    maxrec = ns//every + 2; out = np.zeros(10*maxrec)
    nrec = lib.integrate(len(P), ptr(P), ptr(s), ptr(W), e2_, clam, cpol, ptr(state), ptr(f), ptr(u), ptr(v), dt, ns, every, ptr(out), maxrec, pur_ommax)
    return out[:10*nrec].reshape(nrec, 10)


# ================================ 4. MOTOR EXACTO ================================
def ex_run(x0, dPf=0.04, J=30, Paf=40.0, ncross=10, chunkU=40.0, TmaxU=600.0, dtf=0.2, clam=1.0, cpol=1.0, init='vac', smax_c=14.0, autogrow=True, tag=''):
    sq = np.sqrt(x0); tauU = 1.0/(E_CHARGE*sq)
    while True:
        dP = dPf*sq; Pa = Paf*sq; Pb = 6*sq + 12*x0**(1/3) + 10
        P1 = np.arange(-Pa, Pb + dP/2, dP); smax = smax_c*x0/np.pi
        xg, wg = np.polynomial.legendre.leggauss(J); s1 = 0.5*smax*(xg + 1); ws = 0.5*smax*wg
        P = np.repeat(P1, J); s = np.tile(s1, len(P1)); W = (dP/(2*np.pi))*np.tile(ws/(2*np.pi), len(P1))
        N = len(P); om_max = np.hypot(Pa, np.sqrt(1 + smax)); dt = dtf/om_max
        f = np.zeros(N); u = np.zeros(N); v = np.zeros(N); st = np.array([x0, 0.0])
        if init == 'dressed':                                                        # vacio superadiabatico de 1er orden: R || (0, lam, 2 om)
            om0 = np.sqrt(P**2 + 1 + s); lam0 = x0*np.sqrt(1 + s)/om0**2; phi = np.arctan(lam0/(2*om0)); f[:] = (1 - np.cos(phi))/2; v[:] = np.sin(phi)
        every = max(1, int(0.02*tauU/dt)); ns_chunk = int(chunkU*tauU/dt)
        rec = []; toff = 0.0; ncr = 0; tstart = time.time(); fmax = 0.0; fmin = 0.0
        while toff < TmaxU*tauU and ncr < ncross:
            o = kernel(P, s, W, st, f, u, v, dt, ns_chunk, every, clam=clam, cpol=cpol, pur_ommax=1e300)
            o = o.copy(); o[:, 0] += toff
            rec.append(o if not rec else o[1:]); toff = o[-1, 0]
            allo = np.vstack(rec); ncr = int(np.sum(allo[:-1, 1]*allo[1:, 1] < 0))
            if os.environ.get('QV_VERBOSE'): print(f'   [ex_run x0={x0:g}] t={toff/tauU:.0f}tauU cross={ncr} Phimax/Pa={allo[:,2].max()/Pa:.2f} N={N} wall={time.time()-tstart:.0f}s', flush=True)
            if autogrow and (allo[:, 2].max() > 0.9*Pa or -allo[:, 2].min() > 0.9*Pb): break                  # cobertura insuficiente: abortar y ampliar
        o = np.vstack(rec); runtime = time.time() - tstart
        PhiMax = float(o[:, 2].max()); PhiMin = float(o[:, 2].min())
        edge_left = float(f[:5*J].max()); edge_right = float(f[-5*J:].max())
        cov = dict(Phimax_over_Pa=PhiMax/Pa, minusPhimin_over_Pb=max(-PhiMin, 0)/Pb, edge_f_left=edge_left, edge_f_right=edge_right)
        if autogrow and (cov['Phimax_over_Pa'] > 0.92 or cov['minusPhimin_over_Pb'] > 0.92):                    # ampliar la malla y repetir
            Paf *= max(1.3, 1.15*PhiMax/(0.9*Pa)); continue
        purf = np.abs((1 - 2*f)**2 + u**2 + v**2 - 1.0); pur_w = float(np.sum(W*purf)/np.sum(W))
        return dict(o=o, N=N, dt=dt, pur_w=pur_w, pur_max_end=float(purf.max()), runtime=runtime, cov=cov, Paf=Paf, x0=x0, tauU=tauU, f=f, P=P, s=s, dPf=dPf, J=J)

def crossings(o):
    t = o[:, 0]; x = o[:, 1]; idx = np.where(x[:-1]*x[1:] < 0)[0]
    return np.array([t[i] - x[i]*(t[i+1] - t[i])/(x[i+1] - x[i]) for i in idx]), idx

def ex_analyze(r):
    o = r['o']; x0 = r['x0']; t = o[:, 0]; x = o[:, 1]; cr, idx = crossings(o); k = E_CHARGE*np.sqrt(x0)
    half = np.diff(cr)
    amps = [float(np.abs(x[idx[i]:idx[i+1] + 1]).max()/x0) for i in range(len(idx) - 1)]
    out = dict(n_cross=int(len(cr)), tau0=float(cr[0]), tau0_e_sqrtx0=float(cr[0]*k), half_periods_u=[float(h*k) for h in half], amps=amps,
               omega_early=float(np.pi/half[0]) if len(half) else None, omega_late=float(np.mean(np.pi/half[-3:])) if len(half) >= 3 else None,
               fmax=float(o[:, 6].max()), fmin=float(o[:, 7].min()), Emax_err=float(np.max(np.abs(o[:, 9] - o[0, 9]))/o[0, 9]), purity_dev=float(o[:, 8].max()),
               runtime=r['runtime'], N_modes=int(r['N']), cov=r['cov'], dt=r['dt'], pur_w=r['pur_w'], pur_max_end=r['pur_max_end'])
    out['omega_over_esqrtx0'] = out['omega_late']/k if out.get('omega_late') else None
    if len(amps) >= 3:
        tm = 0.5*(cr[:-1] + cr[1:]); out['gamma'] = float(-np.polyfit(tm, np.log(amps), 1)[0]); out['gamma_over_omega'] = out['gamma']/out['omega_late']
    # ecuacion de estado
    E2 = x**2/(2*e2); rho = o[:, 3]; Pp = o[:, 4]; m = rho > 1e-12*E2[0]
    wpl = np.full(len(t), np.nan); wpl[m] = Pp[m]/rho[m]; wtot = (Pp - E2)/(rho + E2)
    out['wtot_t0'] = float(wtot[0])
    if len(cr) >= 3:
        k0, k1 = idx[-3] + 1, idx[-1] + 1                                                      # un periodo completo entre cruces k y k+2
        out['w_pl_avg_last_period'] = float(np.trapezoid(Pp[k0:k1], t[k0:k1])/np.trapezoid(rho[k0:k1], t[k0:k1]))
        out['w_tot_avg_last_period'] = float(np.trapezoid(Pp[k0:k1] - E2[k0:k1], t[k0:k1])/np.trapezoid(rho[k0:k1] + E2[k0:k1], t[k0:k1]))
        out['plasma_energy_fraction_end'] = float(rho[-1]/(rho[-1] + E2[-1]))
    i0 = int(np.argmax(wtot > 0)) if (wtot > 0).any() else None
    out['t_wtot_cross_0_u'] = float(t[i0]*k) if i0 else None
    out['t_wpl_cross_third_u'] = float(t[np.argmax(wpl > 1/3)]*k) if np.nanmax(wpl) > 1/3 else None
    out['wpl_max'] = float(np.nanmax(wpl))
    return out, dict(t=t*k, x=x/x0, wtot=wtot, wpl=wpl, E2=E2, rho=rho)

