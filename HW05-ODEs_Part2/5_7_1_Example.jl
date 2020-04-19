using Plots
using DifferentialEquations

f(t,x) = [4*exp.(.8*t)] - .5*x
x0 = [2]
tf = 10
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

tf = 5.
x0 = [5 1 2 1.]
include("Adaptive_time_step_Runge_Kutta.jl")
t,x = Adaptive_time_step_Runge_Kutta.vRK(f_af,tf,x0)
plot(t,x[:,1],label = "y",color = "green")
    plot!(t,x[:,2],label = "theta",color = "red")
    plot!(t,x[:,3],label = "ydot",color = "purple")
    plot!(t,x[:,4],label = "thetadot",color = "orange")
    title!("My ODE23")




##
function builtin!(dx,x,p,t)
    dx[1] = x[3]
    dx[2] = x[4]
    dx[3] = -100*x[1]/1
    dx[4] = -100*x[2]/1
end


tspan = (0,tf)
prob = ODEProblem(builtin!,x0',tspan)
sol = solve(prob,BS3())
plot(sol,vars=(0,1),label = "y",color = "green")
    plot!(sol,vars=(0,2),label = "theta",color = "red")
    plot!(sol,vars=(0,3),label = "ydot",color = "purple")
    plot!(sol,vars=(0,4),label = "thetadot",color = "orange")
    title!("Julia's version of ODE23")
##
