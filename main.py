from sympy import *
import re

print("Ax'+By'+Cx+Dy=Ee^(at+b)")
print("A1x'+B1y'+C1x+D1y=E1e^(a1t+b1)")

D, y, x, t, c1, c2, c3, c4 = symbols('D y x t c1 c2 c3 c4')
f = Function('x')(t)

A = 1
B = 0
C = -3
d = 1
F = -1
a = 1
b = 1

A1 = 0
B1 = 1
C1 = -1
d1 = -1
F1 = 4
a1 = 1
b1 = 2
'''
print('A')
A = float(input())
print('B')
B = float(input())
print('C')
C = float(input())
print('D')
d = float(input())
print('E')
F = float(input())
print('a')
a = float(input())
print('b')
b = float(input())

print('A1')
A1 = float(input())
print('B1')
B1 = float(input())
print('C1')
C1 = float(input())
print('D1')
d1 = float(input())
print('E1')
F1 = float(input())
print('a1')
a1 = float(input())
print('b1')
b1 = float(input())
'''

def wronskiano(y1, y2, fx):
    w = simplify(Matrix([[y1, Derivative(y1,t)], [y2, Derivative(y2,t)]]).det())
    u1_d = simplify(-y2*fx/w)
    u2_d = simplify(y1*fx/w)

    u1 = integrate(u1_d, t)
    u2 = integrate(u2_d, t)

    u1y1 = simplify(u1*y1)
    u2y2 = simplify(u2*y2)
    result = simplify(u1y1+u2y2)
    return result

det = simplify(Matrix([[A*D+C,B*D+d],[A1*D+C1,B1*D+d1]]).det())
print(det)

expressions = solve(det)
print(expressions)
exp1 = 0
exp2 = 0
if expressions != []:
    if len(expressions) == 1:
        exp1 = expressions[0]
        exp2 = expressions[0]
    elif len(expressions) == 2:
        exp1 = expressions[0]
        exp2 = expressions[1]
y1 = E**(exp1*t)
y2 = E**(exp2*t)

if exp1 == exp2:
    y1 = E**(exp1*t)
    y2 = t*E**(exp2*t)

x = c1*y1+c2*y2
y = c3*y1+c4*y2

if F == 0 and F1 == 0:
    ec1 = A*Derivative(x, t) + B*Derivative(y, t) + C*x + d*y

    C1 = solve(simplify(ec1.subs(c2,0).subs(c4,0)),c1)[0]
    C2 = solve(simplify(ec1.subs(c1,0).subs(c3,0)),c2)[0]
    C3 = solve(simplify(ec1.subs(c2,0).subs(c4,0)),c3)[0]
    C4 = solve(simplify(ec1.subs(c1,0).subs(c3,0)),c4)[0]

    print(C1)
    print(C2)
    print(C3)
    print(C4)

    new_x = x.subs(c1, C1).subs(c2, C2)

    new_y = y.subs(c3, C3).subs(c4, C4)

    print('x=',x)
    print('y=',new_y)

    print()

    print('x=',new_x)
    print('y=',y)
else:
    y = c1*y1+c2*y2
    detx = expand(simplify(Matrix([[F*E**(a*t+b),F1*E**(a1*t+b1)],[B*D+d,B1*D+d1]]).det()))
    dety = expand(simplify(Matrix([[A*D+C,A1*D+C1], [F*E**(a*t+b),F1*E**(a1*t+b1)]]).det()))
    patron = r"D\*exp\(([^()]+)\)"
    cadena_reemplazo = r"Derivative(exp(\1),t)"
    #print(re.sub(patron, cadena_reemplazo, str(detx)))
    detx = simplify(sympify(re.sub(patron, cadena_reemplazo, str(detx))).subs(D, 0))
    dety = simplify(sympify(re.sub(patron, cadena_reemplazo, str(dety))).subs(D, 0))
    #print(detx)
    #print(dety)
    if 'D' not in str(det):
        x = detx/det
        y = dety/det
        print('x=',x)
        print('y=',y)
        exit()
    #print(x)
    #print(y)
    #print(detx)
    #print(dety)
    #ec = solve(Eq(det, detx), D)
    xp = wronskiano(y1, y2, detx)
    x = expand(simplify(x+xp))
    print('x=',x)

    yp = wronskiano(y1, y2, dety)
    y = expand(simplify(y+yp))
    print('y=',y)
    #wronskiano(y1,y2,dety)
