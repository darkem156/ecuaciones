from sympy import *
import re
import sys
import copy

text_1 = "Ax'+By'+Cx+Dy=Ee^(at+b)"
text_2 = "A1x'+B1y'+C1x+D1y=E1e^(a1t+b1)"
print(text_1)
print(text_2)

D, y, x, t, c1, c2, c3, c4, X, s = symbols('D y x t c1 c2 c3 c4 X s')
f = Function('x')(t)

A = 1
B = 0
C = 1
d = -1
F = 0
a = 0
b = 0

A1 = 0
B1 = 1
C1 = -2
d1 = 0
F1 = 0
a1 = 0
b1 = 0

x0=0
y0=1

# (A)x + (B)y + (C)z ... = D
# D representa otra derivada
# Ejemplo: x''+x'-x+y'=5,x'+y'+y=1
#system = [
#    [D**2+D-1,D,5],
#    [D,D+1,1]
#]

system = [
    [0,Rational(1,100)*D+5,5,100],
    [0,5,Rational(1,80)*D+5,100],
    [1,-1,-1,0]
    #[D+1,1,0,0],
    #[Rational(-3,4),D+Rational(3,2),-3,0],
    #[Rational(-1,8),Rational(-1,4),D+Rational(1,2),0]
]

# Una fila para cada variable. Del orden mas pequeño al mas grande
# Ejemplo:
# x(0)=0, x'(0)=1, x'''(0)=1
# y(0)=1, y'(0)=1, y'''(0)=0
#conditions = [
#    [0, 1, 1],
#    [1, 1, 0]
#]
conditions = [
    [0,0],
    [0,0],
    [0]
]

def replace_D(text):
    if text == 0 or text == "0":
        return 0
    elements = text.as_ordered_terms()#re.split(r"(?<![a-zA-Z]\()([+\-\s=])", text)
    new_elements = []
    for i in range(0, len(elements)):
        if "D" in str(elements[i]):
            match = re.search(r"D\*\*(\d+)", str(elements[i]))
            if match:
                num = match.group(1)
                new_elements.append("Derivative(" + str(elements[i]).replace("D", "1") + ", t, " + num + ")")
            else:
                new_elements.append("Derivative(" + str(elements[i]).replace("D", "1") + ", t)")
        else:
            new_elements.append(str(elements[i]))
    return sympify("+".join(new_elements))

def wronskiano(expressions, fx):
    if fx == 0:
        return 0
    w_matrix = []
    solutions = []
    for expression in expressions:
        current_array = []
        for i in range(0, len(expressions)):
            current_array.append(simplify(Derivative(expression, t, i)))
        w_matrix.append(current_array)
    w = simplify(Matrix(w_matrix).det())
    for i in range(0, len(expressions)):
        wi = copy.deepcopy(w_matrix)
        for j in range(0, len(expressions)):
            wi[i][j] = 0
        wi[i][-1] = fx
        dui = Matrix(wi).det()/w
        ui = expand(integrate(dui, t))
        el = simplify(ui*expressions[i])
        solutions.append(str(el))
    solutions = expand(simplify(sympify("+".join(solutions))))
    return solutions

#det = simplify(Matrix([[A*D+C,B*D+d],[A1*D+C1,B1*D+d1]]).det())

solutions = []

def get_solution(eq, result):
    solutions = roots(eq)
    expressions = []
    for solution in solutions:
        for i in range(0, solutions[solution]):
            #expressions.append('c'+str(i+1)+'*'+str(t**i*exp(solution*t)))
            if 'I' in str(solution):
                if im(solution).subs(I,1) < 0:
                    expressions.append((E**(solution.subs(I, 0)*t))*cos(im(solution).subs(I, 1)*t))
                else:
                    expressions.append((E**(solution.subs(I, 0)*t))*sin(im(solution).subs(I, 1)*t))
            else:
                expressions.append(t**i*exp(solution*t))
    extras = wronskiano(expressions, result)
    for i in range(0, len(expressions)):
        expressions[i] = 'c' + str(i+1) + '*' + str(expressions[i])
    return sympify("+".join(expressions))+extras

def sort_d(matrix, index):
    i = index if index < len(matrix) else index-1
    new_matrix = []
    d_matrix = []
    for el in matrix:
        if 'D' in str(el[i]):
            d_matrix.append(el)
        else:
            new_matrix.append(el)
    d_matrix = sorted(d_matrix, key=lambda x: simplify(sympify(x[i]).subs(D, 2)), reverse=True)
    return new_matrix+d_matrix

vars = len(system)
def solve_system(system):
    global solutions, vars
    results = [row.pop() for row in system]
    for i in range(0, len(system[0])-len(system)):
        system.append([1 for i in range(0, len(system[0]))])
    det = expand(simplify(Matrix(system).det()))

    current_matrix = copy.deepcopy(system)
    for i in range(0, len(current_matrix)):
        current_matrix[i][0] = results[i]
    current_det = simplify(expand(simplify(replace_D(Matrix(current_matrix).det()))))

    solution = 0
    if 'D' not in str(det):
        solution = expand(current_det)/det
    else:
        are_all_negative = all(['-' in str(term) for term in det.as_ordered_terms()])
        if are_all_negative:
            current_det = -current_det
            det = -det
            print(det, current_det)
        solution = get_solution(det, expand(current_det))
        '''
        eq = simplify(Eq(det*x, current_det))
        if 'D' in str(eq.lhs):
            solution = get_solution(eq.lhs.subs(x,1), eq.rhs.subs(x,1))
        else:
            solution = get_solution(eq.rhs.subs(x,1), eq.lhs.subs(x,1))
        '''
    solutions.append(solution)
    if len(solutions) < vars:
        i = 0
        subsystem = copy.deepcopy(system)
        for index in range(0, len(subsystem)):
            results[index] = results[index] - simplify(replace_D(expand(solution * subsystem[index][i])))
            for j in range(0, len(subsystem[index])):
                if i == j:
                    #subsystem[index][j] = simplify(replace_D(expand(solution * subsystem[index][j])))
                    subsystem[index][j] = 0
            subsystem[index].append(expand(results[index]))
            subsystem[index] = subsystem[index][1:]
        subsystem = sort_d([eq for eq in subsystem if eq[0] != 0], i+2)[0:len(subsystem)-i-1]
        print('subsystem', subsystem)
        print('continue')
        solve_system(subsystem)

if conditions == []:
    solve_system(system)
    print(solutions)
else:
    laplaces = []
    results = []
    for n in range(0, len(system)):
        ec = system[n]
        result = ec.pop()
        result = simplify(laplace_transform(result, t, s, noconds=True))
        results.append(result)
        eq = ''
        for i in range(1, len(ec)+1):
            current_eq = sympify(ec[i-1] * Function('x'+str(i))(t))
            eq = str(eq) + '+' + str(current_eq)
        eq = simplify(replace_D(expand(sympify(eq))))
        laplace_eq = str(laplace_transform(eq, t, s, noconds=True))
        for i in range(1, len(ec)+1):
            laplace_eq = laplace_eq.replace('x'+str(i)+'(0)', str(conditions[i-1][0]))
        laplace_eq = sympify(laplace_eq)
        elements = laplace_eq.as_ordered_terms()
        for i in range(1, len(ec)+1):
            subs_elements = []
            for j in range(0, len(elements)):
                element = elements[j]
                if "Subs" in str(element) and ("x"+str(i)) in str(element):
                    subs_elements.append(element)
                    elements[j] = 0
            for j in range(0, len(subs_elements)):
                elements.append(-conditions[i-1][j+1])
        elements = sympify("+".join([str(element) for element in elements]))
        coeficientes = elements.as_coefficients_dict()
        for coeficiente in coeficientes:
            if 'LaplaceTransform' not in str(coeficiente):
                num = coeficiente * coeficientes[coeficiente]
                results[n] = results[n] - num
                elements = elements - num
        #print(elements, '=', results[n])
        s_matrix = []
        for i in range(1, len(ec)+1):
            s_matrix.append(sympify("+".join([str(element) if "x"+str(i) in str(element) else "0" for element in elements.as_ordered_terms()]).replace("LaplaceTransform(x" + str(i) + "(t), t, s)", "1")))
            #elements = sympify(str(elements).replace('LaplaceTransform(x' + str(i) + '(t), t, s)', "1"))
        laplaces.append(s_matrix)
    det = Matrix(laplaces).det()
    print(laplaces)
    print(results)
    for i in range(0, len(system)):
        current_matrix = copy.deepcopy(laplaces)
        for j in range(0, len(current_matrix)):
            current_matrix[j][i] = results[j]
        current_det = Matrix(current_matrix).det()
        current_var = expand(current_det/det)
        current_var = expand(simplify(sympify(str(inverse_laplace_transform(current_var, s, t)).replace('Heaviside(t)', "1"))))
        print('x'+str(i+1) + ' =', current_var)
exit()
exp1 = 0
exp2 = 0
if expressions != []:
    if "I" not in str(expressions):
        expressions.sort(reverse=False)
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
    y1 = (E**(exp1.subs(I, 0)*t))*cos(im(exp1).subs(I, 1)*t)
if "I" in str(exp2):
    y2 = (E**(exp2.subs(I, 0)*t))*sin(im(exp2).subs(I, 1)*t)

x = c1*y1+c2*y2
y = c1*y1+c2*y2

if x0 != None and y0 != None:
    x = Function('x')(t)
    y = Function('y')(t)
    s, t = symbols('s t')
    ec1 = A*Derivative(x, t) + B*Derivative(y, t) + C*x + d*y-(F*E**(a*t+b))
    ec2 = A1*Derivative(x, t) + B1*Derivative(y, t) + C1*x + d1*y-(F1*E**(a1*t+b1))
    laplace_ec1 = laplace_transform(ec1, t, s, noconds=True)
    laplace_ec1 = sympify(str(laplace_ec1).replace("x(0)", str(x0)).replace("y(0)", str(y0)))
    coeficientes = laplace_ec1.as_coefficients_dict()
    coeficientes_matrix = [0, 0]
    if Integer(1) in coeficientes:
        coeficientes_matrix[0] = -1*coeficientes.pop(Integer(1))
    laplace_ec1 = sympify(str(laplace_ec1).replace(" ", "").replace(str(coeficientes_matrix[0]), "0"))
    laplace_ec2 = laplace_transform(ec2, t, s, noconds=True)
    laplace_ec2 = sympify(str(laplace_ec2).replace("x(0)", str(x0)).replace("y(0)", str(y0)))
    coeficientes = laplace_ec2.as_coefficients_dict()
    if Integer(1) in coeficientes:
        coeficientes_matrix[1] = -1*coeficientes.pop(Integer(1))
    laplace_ec2 = sympify(str(laplace_ec2).replace(" ", "").replace(str(coeficientes_matrix[1]), "0"))

    xs = [
        sympify(str(laplace_ec1).replace("LaplaceTransform(x(t), t, s)", "1").replace("LaplaceTransform(y(t), t, s)", "0")),
        sympify(str(laplace_ec2).replace("LaplaceTransform(x(t), t, s)", "1").replace("LaplaceTransform(y(t), t, s)", "0")),
    ]
    ys = [
        sympify(str(laplace_ec1).replace("LaplaceTransform(x(t), t, s)", "0").replace("LaplaceTransform(y(t), t, s)", "1")),
        sympify(str(laplace_ec2).replace("LaplaceTransform(x(t), t, s)", "0").replace("LaplaceTransform(y(t), t, s)", "1")),
    ]

    det = Matrix([xs, ys]).det()
    print(det)
    detx = Matrix([coeficientes_matrix, ys]).det()
    xt = sympify(str(inverse_laplace_transform(detx/det, s, t)).replace("Heaviside(t)", "1"))
    dety = Matrix([xs, coeficientes_matrix]).det()
    yt = sympify(str(inverse_laplace_transform(dety/det, s, t)).replace("Heaviside(t)", "1"))
    print(xt)
    print(yt)
