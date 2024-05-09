from sympy import *
import re

text_1 = "Ax'+By'+Cx+Dy=Ee^(at+b)"
text_2 = "A1x'+B1y'+C1x+D1y=E1e^(a1t+b1)"
print(text_1)
print(text_2)

D, y, x, t, c1, c2, c3, c4, X = symbols('D y x t c1 c2 c3 c4 X')
f = Function('x')(t)

'''
A = 1
B=0
C=0
d=-4
F=1
a=0
b=0

A1 = 0
B1=1
C1=1
d1=0
F1=2
a1=0
b1=0
'''
A = 1
B = 0
C = 0
d = -4
F = 1
a = 0
b = 0

A1 = 0
B1 = 1
C1 = 1
d1 = 0
F1 = 2
a1 = 0
b1 = 0
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

expressions = solve(det)
exp1 = 0
exp2 = 0
if expressions != []:
    if len(expressions) == 1:
        exp1 = expressions[0]
        exp2 = expressions[0]
    elif len(expressions) >= 2:
        exp1 = expressions[0]
        exp2 = expressions[1]
y1 = E**(exp1*t)
y2 = E**(exp2*t)

if exp1 == exp2:
    y1 = E**(exp1*t)
    y2 = t*E**(exp2*t)

if "I" in str(exp1):
    y1 = (E**(exp1.subs(I, 0)*t))*cos(exp1.subs(I, t))
if "I" in str(exp2):
    y2 = (E**(exp1.subs(I, 0)*t))*sin(exp2.subs(I, t))

x = c1*y1+c2*y2
y = c3*y1+c4*y2

if F == 0 and F1 == 0:
    ec1 = A*Derivative(x, t) + B*Derivative(y, t) + C*x + d*y

    C1 = solve(simplify(ec1.subs(c2,0).subs(c4,0)),c1)[0]
    C2 = solve(simplify(ec1.subs(c1,0).subs(c3,0)),c2)[0]
    C3 = solve(simplify(ec1.subs(c2,0).subs(c4,0)),c3)[0]
    C4 = solve(simplify(ec1.subs(c1,0).subs(c3,0)),c4)[0]

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
    wronskiano_y = expand(simplify(y+yp))
    ec1 = Eq(A*Derivative(x, t) + B*Derivative(y, t) + C*x + d*y,F*E**(a*t+b))
    solve_x = solve(ec1, X)
    solve_dx = solve(ec1, Derivative(X, t))
    if solve_x != []:
        y = expand(solve_x[0])
    elif solve_dx != []:
        solve_dx = integrate(solve_dx[0], t)+c2
        y = expand(solve_dx)

    if 'X' in str(y) or 'Derivative(X, t)' in str(y) or (solve_x == [] and solve_dx == []):
        ec2 = Eq(A1*Derivative(x, t) + B1*Derivative(X, t) + C1*x + d1*X,F1*E**(a1*t+b1))
        solve_x_1 = solve(ec2, X)
        solve_dx_1 = solve(ec2, Derivative(X, t))
        if solve_x_1 != []:
            y = expand(solve_x_1[0])
        elif solve_dx_1 != []:
            solve_dx_1 = integrate(solve_dx_1[0], t)+c2
            y = expand(solve_dx_1)

    if 'X' in str(y) or 'Derivative(X, t)' in str(y):
        y = wronskiano_y
    print('y=',y)

    print('Sustituir en la parte homogenea de la primera ecuacion:')
    result1 = simplify(A*Derivative(x, t) + B*Derivative(y, t) + C*x + d*y)
    print(simplify(result1))
    if "c1" in str(simplify(result1)):
        c_value = solve(Eq(result1,F*E**(a*t+b)), c1)
        if c_value!= []:
            print('Resultado actualizado:')
            x = expand(x.subs(c1,c_value[0]))
            y = expand(y.subs(c1,c_value[0]))
            print('x=',x)
            print('y=',y)
        result1 = result1.subs(c1, c_value[0])
    if 'c2' in str(simplify(result1)):
        c_value = solve(Eq(result1,F*E**(a*t+b)), c2)
        if c_value != []:
            print('Resultado actualizado:')
            x = expand(x.subs(c2,c_value[0]))
            y = expand(y.subs(c2,c_value[0]))
            print('x=',x)
            print('y=',y)
        result1 = result1.subs(c2, c_value[0])
    print('Sustituir en la parte homogenea de la segunda ecuacion:')
    result2 = simplify(A1*Derivative(x, t) + B1*Derivative(y, t) + C1*x + d1*y)
    print(simplify(result2))
    if 'c1' in str(simplify(result2)):
        c_value = solve(Eq(result2,F1*E**(a1*t+b1)), c1)
        if c_value != []:
            print('Resultado actualizado:')
            x = expand(x.subs(c1,c_value[0]))
            y = expand(y.subs(c1,c_value[0]))
            print('x=',x)
            print('y=',y)
        result2 = result2.subs(c2, c_value[0])
    if 'c2' in str(simplify(result2)):
        c_value = solve(Eq(result2,F1*E**(a1*t+b1)), c2)
        if c_value != []:
            print('Resultado actualizado:')
            x = expand(x.subs(c2,c_value[0]))
            y = expand(y.subs(c2,c_value[0]))
            print('x=',x)
            print('y=',y)
        result2 = result2.subs(c2, c_value[0])

    if expand(result1) == expand(F*E**(a*t+b)) and expand(result2) == expand(F1*E**(a1*t+b1)):
        print('El resultado es exacto')
        print('x =',x)
        print('y =',y)
