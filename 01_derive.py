"""EJECUTOR FORMAL (CAS): deriva y VERIFICA simbolicamente el sistema Einstein-Euler-Heisenberg
con simetria esferica estatica.  Si alguna identidad no se reduce a 0 -> exit(1) (deduccion descartada).
Ansatz: ds^2=-f dt^2+dr^2/g+r^2 dOmega^2,  A=Phi(r)dt,  B=0.
"""
import sys
import sympy as sp

t, r, th, ph = sp.symbols('t r theta phi', positive=True)
G, c = sp.symbols('G c', positive=True)           # c == alpha^2/(90 m^4) (coef. EH)
f, g, Phi = (sp.Function(n)(r) for n in ('f', 'g', 'Phi'))
x = [t, r, th, ph]

gdd = sp.diag(-f, 1/g, r**2, r**2*sp.sin(th)**2)
guu = gdd.inv()
sqrtg = sp.sqrt(f/g) * r**2 * sp.sin(th)          # sqrt(-det g)

# --- geometria -------------------------------------------------------------
Gam = [[[sum(guu[a, d]*(sp.diff(gdd[d, b], x[cc]) + sp.diff(gdd[d, cc], x[b]) - sp.diff(gdd[b, cc], x[d]))
             for d in range(4))/2 for cc in range(4)] for b in range(4)] for a in range(4)]
def Ricci(b, cc):
    return sum(sp.diff(Gam[a][b][cc], x[a]) - sp.diff(Gam[a][b][a], x[cc])
               + sum(Gam[a][a][d]*Gam[d][b][cc] - Gam[a][cc][d]*Gam[d][b][a] for d in range(4))
               for a in range(4))
Ric = sp.Matrix(4, 4, lambda i, j: sp.simplify(Ricci(i, j)))
Rs = sp.simplify(sum(guu[i, j]*Ric[i, j] for i in range(4) for j in range(4)))
Gdd = Ric - Rs*gdd/2
Gud = sp.simplify(guu*Gdd)                        # G^mu_nu

# --- materia: L(X,Y) = -X/4 + c (X^2 + 7/4 Y^2),  X=F_{mn}F^{mn}, Y=F_{mn} Ftilde^{mn} ---
A = [Phi, 0, 0, 0]
Fdd = sp.Matrix(4, 4, lambda i, j: sp.diff(A[j], x[i]) - sp.diff(A[i], x[j]))
Fuu = guu*Fdd*guu.T
X = sp.simplify(sum(Fdd[i, j]*Fuu[i, j] for i in range(4) for j in range(4)))
eps = lambda *idx: sp.LeviCivita(*idx)
Ftil_uu = sp.Matrix(4, 4, lambda a, b: sum(-eps(a, b, cc, d)/sqrtg*Fdd[cc, d]/2
                                           for cc in range(4) for d in range(4)))
Y = sp.simplify(sum(Fdd[i, j]*Ftil_uu[i, j] for i in range(4) for j in range(4)))
Xs, Ys = sp.symbols('X Y')
Lsym = -Xs/4 + c*(Xs**2 + sp.Rational(7, 4)*Ys**2)
LX, LY = sp.diff(Lsym, Xs), sp.diff(Lsym, Ys)
sub = {Xs: X, Ys: Y}
L, LX, LY = Lsym.subs(sub), LX.subs(sub), LY.subs(sub)
FF_dd = Fdd*gdd.inv()*Fdd.T                        # F_{ma}F_n^a  (F_{ma} g^{ab} F_{nb})
Tdd = -4*LX*FF_dd - LY*Y*gdd + gdd*L               # T_{mn} = -2 dL/dg^{mn} + g_{mn} L
Tud = sp.simplify(guu*Tdd)

print("X = F_mn F^mn =", X, "\nY = F Ftilde =", Y)
print("T^t_t =", sp.simplify(Tud[0, 0]), "\nT^r_r =", sp.simplify(Tud[1, 1]),
      "\nT^th_th =", sp.simplify(Tud[2, 2]))

# --- ecuaciones de Einstein ----------------------------------------------------
Ett, Err, Eth = (sp.simplify(Gud[i, i] - 8*sp.pi*G*Tud[i, i]) for i in (0, 1, 2))
print("\n[G^t_t - G^r_r] =", sp.simplify(Gud[0, 0] - Gud[1, 1]), "   (T^t_t - T^r_r =",
      sp.simplify(Tud[0, 0] - Tud[1, 1]), ")")

# --- accion reducida y Euler-Lagrange (verificacion independiente) --------------
Lred = sp.simplify((sqrtg/sp.sin(th))*(Rs/(16*sp.pi*G) + L))
from sympy.calculus.euler import euler_equations
EL = euler_equations(Lred, [f, g, Phi], r)
EL = [sp.simplify(e.lhs) for e in EL]

# --- sustitucion en-capa: f=g=1-y/r, Phi'=E, E(1+16 c E^2) r^2 = Q/4pi (integral primera) ---
y_, E_, r_ = sp.symbols('y E r_', positive=True)
yp = None
rho = E_**2/2 + 12*c*E_**4
yp = 8*sp.pi*G*r**2*rho                            # y' = 8 pi G r^2 rho,  y = 2 G m(r)
Ep = -2*E_*(1 + 16*c*E_**2)/(r*(1 + 48*c*E_**2))   # derivada de E(r) por la integral primera
def D(h):                                          # derivada total en la capa
    return sp.diff(h, r) + sp.diff(h, y_)*yp + sp.diff(h, E_)*Ep

def onshell(expr):
    y = sp.Function('y')(r)
    e = expr.subs({f: 1 - y/r, g: 1 - y/r}).doit()
    # Phi^(k) -> D^{k-1}(E) ; y^(k) -> D^{k-1}(yp)
    dPhi = {1: E_}
    dy = {1: yp}
    for k in range(2, 6):
        dPhi[k] = sp.simplify(D(dPhi[k-1]))
        dy[k] = sp.simplify(D(dy[k-1]))
    for k in range(5, 0, -1):
        e = e.subs(sp.Derivative(y, (r, k)), dy[k]) if k > 1 else e.subs(sp.Derivative(y, r), dy[1])
        e = e.subs(sp.Derivative(Phi, (r, k)), dPhi[k]) if k > 1 else e.subs(sp.Derivative(Phi, r), dPhi[1])
    return sp.simplify(e.subs(y, y_))

checks = {'Einstein tt': Ett, 'Einstein rr': Err, 'Einstein thth': Eth,
          'E-L delta f': EL[0], 'E-L delta g': EL[1], 'E-L delta Phi': EL[2]}
ok = True
print("\n=== Verificacion simbolica en-capa (debe dar 0) ===")
for name, e in checks.items():
    res = onshell(e)
    flag = (res == 0)
    ok &= flag
    print(f"{name:16s} -> {res}   {'OK' if flag else 'FALLA'}")

print("\n=== Sistema reducido deducido ===")
print("Maxwell (1a integral):  r^2 E (1 + 16 c E^2) = Q/(4 pi),   E = sqrt(g/f) Phi'")
print("f = g   (T^t_t = T^r_r  =>  (ln f/g)' = 0)")
print("y' = 8 pi G r^2 rho(E),  rho = E^2/2 + 12 c E^4,  g = 1 - y/r,  y = 2 G m(r)")
print("p_r = -rho,   p_perp = E^2/2 + 4 c E^4")
if not ok:
    print("\nRESULTADO: DEDUCCION DESCARTADA"); sys.exit(1)
print("\nRESULTADO: todas las identidades se reducen a 0 simbolicamente.")
