using Plots
using DifferentialEquations

f(t,x) = [4*exp.(.8*t)] - .5*x
x0 = [2]
tf = 2
tol = 1e-8


include("Adaptive_time_step_Runge_Kutta.jl")

t,x = Adaptive_time_step_Runge_Kutta.vRK(f,tf,x0;tol = .1)
plot(t,x[:,1])


##

# J*thetaddot + kt * theta = 0

# k*y + m*yddot = 0

# x1 = y
# x2 = theta
# x3 = ydot
# x4 = thetadot


# THE FOUR EQUATIONS TO SOLVE:
# x1dot = x3
# x2dot = x4
# x3dot = -k*x1/m
# x4dot = -kt*x2/J
m = 1
J = 1
k = 100
kt = 100
f_af(t,x) = [x[3] x[4] -k*x[1]/m -kt*x[2]/J]

f_af(0,[1 2 3 4])

function builtin!(dx,x,p,t)
    p = k, m, kt, J
    dx = [x[3]; x[4]; -k*x[1]/m; -kt*x[2]/J]
end

prob = ODEProblem(builtin!,u₀,tspan,M)
sol = solve(prob)
