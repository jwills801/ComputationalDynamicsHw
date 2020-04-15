using Plots
using SymPy
@vars p q theta t

Ham = p^2/2 +cos(q)


# now I'll find the nondimentional momentum of the angle
pdot = -diff(Ham,q) |> subs(q=>theta)
thetadot = diff(Ham,p)

# therefore thetaddot = sin(theta)
# we need to get it in first order form
# thetaddot = x1dot
# thetadot = x2dot

# Now, there are two first order ODEs to solve.
# 1) x1dot = sin(x2)
# 2) x2dot = x1




include("BackwardEuler.jl")


f(x1,x2) = [sin(x2);x1]
h = .01
x0 = [0 1.]
tf = 5

x,t = BackwardEuler.beuler(f,tf,h,x0)
plot(t,x[:,1],label = "qdot")
    plot!(t,x[:,2],label = "q")
    xlabel!("Time")
    ylabel!("State")

##

include("Adams_Bashforth2.jl")
x,t = Adams_Bashforth2.ab2(f,tf,h,x0)
plot(t,x[:,1],label = "qdot")
    plot!(t,x[:,2],label = "q")
    xlabel!("Time")
    ylabel!("State")

##
include("GaussLegendreRungeKutta.jl")

x,t = GaussLegendreRungeKutta.IRK(f,tf,h,x0)
plot(t,x[:,1],label = "qdot")
    plot!(t,x[:,2],label = "q")
    xlabel!("Time")
    ylabel!("State")

fill(1.,502,1)

oftype(x[:,1],fill(1.,502,1))

x[:,1]
