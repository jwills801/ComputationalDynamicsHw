


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








f(x) = [sin(x[2]) x[1]]
h = .01
x0 = [0 1.]
tf = 10

include("BackwardEuler.jl")
x_BE,t = BackwardEuler.beuler(f,tf,h,x0)
plot(t,x_BE[:,1],label = "qdot")
    plot!(t,x_BE[:,2],label = "q")
    xlabel!("Time")
    ylabel!("State")
    title!("Backwards Euler")

##

include("Adams_Bashforth2.jl")
x_AB,t = Adams_Bashforth2.ab2(f,tf,h,x0)
plot(t,x_AB[:,1],label = "qdot")
    plot!(t,x_AB[:,2],label = "q")
    xlabel!("Time")
    ylabel!("State")
    title!("Adams Bashforth")

##
include("GL_RK.jl")
x_GL,t = GL_RK.gl_rk(f,tf,h,x0)
plot(t,x_GL[:,1],label = "qdot")
    plot!(t,x_GL[:,2],label = "q")
    xlabel!("Time")
    ylabel!("State")
    title!("Gauss Legendre Runge Kutta")

##
H_BE = x_BE[:,1].^2/2 + cos.(x_BE[:,2])
H0 = x_BE[1,1].^2/2 + cos.(x_BE[1,2])

H_AB = x_AB[:,1].^2/2 +cos.(x_AB[:,2])

H_GL = x_GL[:,1].^2/2 +cos.(x_GL[:,2])

plot(t,H_BE/H0,label = "Backward Euler")
    plot!(t,H_AB/H0,label = "Adams Bashforth")
    plot!(t,H_GL/H0,label = "Gauss Legendre Runge Kutta")
    xlabel!("Time")
    ylabel!("Normalized Hamiltonian")
    title!("How Energy Changes with Time if Step Size is .01 Seconds")
