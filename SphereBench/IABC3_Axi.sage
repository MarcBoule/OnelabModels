# Calculate 3rd order IABC relative permittivities in tables 5 to 8 of paper:
# Improvised Asymptotic Boundary Conditions for Electrostatic Finite Elements
# by David C. Meeker
# in IEEE TRANSACTIONS ON MAGNETICS, VOL. 50, NO. 6, JUNE 2014

var('n','rho','R','d','c11','c21','c12','c22','c13','c23','e0','e1','e2','e3','neum')

neum = 1     # 0 = Dirichlet, 1 = Neumann     
R = 1        # problem is scale independant, with R = 1 we have delta = d/R = d
d = 1/3      # with above, this is now delta
e0 = 1       # permittivities are now relative permittivities

u1(n, rho) = c11*rho^n + c21*rho^(-1 - n)
u2(n, rho) = c12*rho^n + c22*rho^(-1 - n)
u3(n, rho) = c13*rho^n + c23*rho^(-1 - n)

du1(n, rho) = diff(u1(n,rho), rho)
du2(n, rho) = diff(u2(n,rho), rho)
du3(n, rho) = diff(u3(n,rho), rho)

eq1 = u1(n, R) == 1
if neum == 1:
    eq2 = du3(n, R + 3*d) == 0
else:
    eq2 = u3(n, R + 3*d) == 0
eq3 = u1(n, R + d) == u2(n, R + d)
eq4 = e1*du1(n, R + d) == e2*du2(n, R + d)
eq5 = u2(n, R + 2*d) == u3(n, R + 2*d)
eq6 = e2*du2(n, R + 2*d) == e3*du3(n, R + 2*d)

sol1 = solve([eq1, eq2, eq3, eq4, eq5, eq6], c11, c21, c12, c22, c13, c23)

u1(n, rho) = u1(n, rho).subs(sol1)
du1(n, rho) = du1(n, rho).subs(sol1)
eqA = e1*du1(0+neum, R) + (1+neum)*e0*u1(0+neum, R)/R
eqB = e1*du1(1+neum, R) + (2+neum)*e0*u1(1+neum, R)/R
eqC = e1*du1(2+neum, R) + (3+neum)*e0*u1(2+neum, R)/R
eqA = eqA.factor()
eqB = eqB.factor()
eqC = eqC.factor()

dA = eqA.denominator()
fA = eqA.numerator() #+ dA  # already included
dB = eqB.denominator()
fB = eqB.numerator() #+ dB  # already included
dC = eqC.denominator()
fC = eqC.numerator() #+ dC  # already included

R.<e1, e2, e3> = PolynomialRing(QQ)
I = Ideal([R(fA), R(fB), R(fC)])

D = R(dA) * R(dB) * R(dC)
I_sat = I.saturation(Ideal([D]))[0]

G = I_sat.groebner_basis()
sols = I_sat.variety(QQbar)
for i, s in enumerate(sols, 1):
    if (s[e1] < 0 or s[e2] < 0 or s[e3] < 0): continue
    print({e1: CC(s[e1]), e2: CC(s[e2]), e3: CC(s[e3])})
