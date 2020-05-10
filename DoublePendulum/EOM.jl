using SymPy
using LinearAlgebra

@vars t L1 L2 m1 m2 g x1 x2 x3 x4
theta1 = SymFunction("theta1")(t)
theta2 = SymFunction("theta2")(t)

r1 = [L1*sin(theta1); -L1*cos(theta1); 0]
r2 = r1 + [L2*sin(theta2);  -L2*cos(theta2); 0]

r1dot = diff.(r1,t)
r2dot = diff.(r2,t)

omega1 = [0; 0; diff(theta1,t)]
omega2 = [0; 0; diff(theta2,t)]

v1c = cross(omega1,r1)
v2c = v1c + cross(omega2,r2)

J = [m1*L1^2, m2*L2^2]
T = 1//2*(v1c.T*m1*v1c+omega1.T*J[1]*omega1
    + v2c.T*m2*v2c+omega2.T*J[2]*omega2)

V = 1//2*(m1*g*r1[2] + m2*g*r2[2])

L = T[1]-V

EOM1 = diff(diff(L,diff(theta1,t)),t)-diff(L,theta1)
EOM2 = diff(diff(L,diff(theta2,t)),t)-diff(L,theta2)

dictionary = solve([EOM1; EOM2],[diff(diff(theta1,t),t); diff(diff(theta2,t),t)])

eq1 = dictionary[diff(diff(theta1,t),t)]
eq2 = dictionary[diff(diff(theta2,t),t)]

changeofvariables = Dict(theta1=>x1,theta2=>x2,diff(theta1,t)=>x3,diff(theta2,t)=>x4)

dx3 = eq1.subs(changeofvariables)
sympy.julia_code(dx3) |> println

dx4 = eq2.subs(changeofvariables)
sympy.julia_code(dx4) |> println
