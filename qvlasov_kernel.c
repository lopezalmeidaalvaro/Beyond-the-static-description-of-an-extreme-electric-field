/* Nucleo C (OpenMP) del sistema exacto de Bloch/Bogoliubov acoplado a Maxwell (homogeneo, QED escalar->espinorial).
   Modo i: (P_i, s_i=p_perp^2, W_i).  p=P+Phi, om=sqrt(p^2+1+s), eperp=sqrt(1+s), lam=c_lam*x*eperp/om^2
     f'=lam u/2 ; u'=lam(1-2f)-2 om v ; v'=2 om u ; Phi'=x ; x'=-2 e^2 sum W [ (p/om) f + c_pol (eperp/(2 om)) u ]
   c_lam=1 (estandar, reproduce Schwinger) ; c_pol=1 consistente con c_lam=1 (conserva la energia) ; flags para controles.
   rho_pl = 2 sum W om f ; P_par = sum W [2 p^2 f/om + p eperp u/om]. RK4 de paso fijo. */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
typedef struct { double x, Phi; } Glob;

static void stage(long N, const double *P, const double *s, const double *W, double e2, double clam, double cpol,
                  double x, double Phi, const double *f, const double *u, const double *v,
                  double *kf, double *ku, double *kv, double *dx)
{
    double acc = 0.0;
    #pragma omp parallel for reduction(+:acc) schedule(static)
    for (long i = 0; i < N; i++) {
        double p = P[i] + Phi, ep2 = 1.0 + s[i], om2 = p*p + ep2, om = sqrt(om2), ep = sqrt(ep2);
        double lam = clam * x * ep / om2;
        kf[i] = 0.5*lam*u[i];
        ku[i] = lam*(1.0 - 2.0*f[i]) - 2.0*om*v[i];
        kv[i] = 2.0*om*u[i];
        acc += W[i] * ((p/om)*f[i] + cpol*(ep/(2.0*om))*u[i]);
    }
    *dx = -2.0*e2*acc;
}

/* integra nsteps pasos de tamano dt; guarda cada 'every' pasos: t,x,Phi,rho,Ppar,npairs,fmax,fmin,purity_dev,energia en out[10*nrec] */
long integrate(long N, const double *P, const double *s, const double *W, double e2, double clam, double cpol,
               double *state /* [x,Phi] */, double *f, double *u, double *v, double dt, long nsteps, long every, double *out, long maxrec, double pur_ommax)
{
    double *f1 = malloc(N*sizeof(double)), *u1 = malloc(N*sizeof(double)), *v1 = malloc(N*sizeof(double));
    double *kf[4], *ku[4], *kv[4]; for (int k = 0; k < 4; k++) { kf[k] = malloc(N*sizeof(double)); ku[k] = malloc(N*sizeof(double)); kv[k] = malloc(N*sizeof(double)); }
    double x = state[0], Phi = state[1], t = 0.0; long nrec = 0;
    for (long n = 0; n <= nsteps; n++) {
        if (n % every == 0 && nrec < maxrec) {
            double rho = 0, Pp = 0, np_ = 0, fmx = -1e300, fmn = 1e300, pur = 0;
            #pragma omp parallel for reduction(+:rho,Pp,np_) reduction(max:fmx,pur) reduction(min:fmn) schedule(static)
            for (long i = 0; i < N; i++) {
                double p = P[i] + Phi, ep = sqrt(1.0 + s[i]), om = sqrt(p*p + ep*ep);
                rho += 2.0*W[i]*om*f[i]; Pp += W[i]*(2.0*p*p*f[i]/om + p*ep*u[i]/om); np_ += W[i]*f[i];
                if (f[i] > fmx) fmx = f[i]; if (f[i] < fmn) fmn = f[i];
                if (om < pur_ommax) { double d = fabs((1-2*f[i])*(1-2*f[i]) + u[i]*u[i] + v[i]*v[i] - 1.0); if (d > pur) pur = d; }
            }
            double *o = out + 10*nrec; o[0] = t; o[1] = x; o[2] = Phi; o[3] = rho; o[4] = Pp; o[5] = np_; o[6] = fmx; o[7] = fmn; o[8] = pur; o[9] = x*x/(2.0*e2) + rho; nrec++;
        }
        if (n == nsteps) break;
        double dx[4];
        stage(N,P,s,W,e2,clam,cpol, x, Phi, f,u,v, kf[0],ku[0],kv[0], &dx[0]);
        double x1 = x + 0.5*dt*dx[0], Ph1 = Phi + 0.5*dt*x;
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N; i++) { f1[i] = f[i] + 0.5*dt*kf[0][i]; u1[i] = u[i] + 0.5*dt*ku[0][i]; v1[i] = v[i] + 0.5*dt*kv[0][i]; }
        stage(N,P,s,W,e2,clam,cpol, x1, Ph1, f1,u1,v1, kf[1],ku[1],kv[1], &dx[1]);
        double x2 = x + 0.5*dt*dx[1], Ph2 = Phi + 0.5*dt*x1;
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N; i++) { f1[i] = f[i] + 0.5*dt*kf[1][i]; u1[i] = u[i] + 0.5*dt*ku[1][i]; v1[i] = v[i] + 0.5*dt*kv[1][i]; }
        stage(N,P,s,W,e2,clam,cpol, x2, Ph2, f1,u1,v1, kf[2],ku[2],kv[2], &dx[2]);
        double x3 = x + dt*dx[2], Ph3 = Phi + dt*x2;
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N; i++) { f1[i] = f[i] + dt*kf[2][i]; u1[i] = u[i] + dt*ku[2][i]; v1[i] = v[i] + dt*kv[2][i]; }
        stage(N,P,s,W,e2,clam,cpol, x3, Ph3, f1,u1,v1, kf[3],ku[3],kv[3], &dx[3]);
        double xn = x + dt/6.0*(dx[0] + 2*dx[1] + 2*dx[2] + dx[3]);
        double Phn = Phi + dt/6.0*(x + 2*x1 + 2*x2 + x3);
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < N; i++) {
            f[i] += dt/6.0*(kf[0][i] + 2*kf[1][i] + 2*kf[2][i] + kf[3][i]);
            u[i] += dt/6.0*(ku[0][i] + 2*ku[1][i] + 2*ku[2][i] + ku[3][i]);
            v[i] += dt/6.0*(kv[0][i] + 2*kv[1][i] + 2*kv[2][i] + kv[3][i]);
        }
        x = xn; Phi = Phn; t += dt;
    }
    state[0] = x; state[1] = Phi;
    free(f1); free(u1); free(v1); for (int k = 0; k < 4; k++) { free(kf[k]); free(ku[k]); free(kv[k]); }
    return nrec;
}
