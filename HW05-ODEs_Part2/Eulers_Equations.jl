
y0 = [cos(1.1) 0 sin(1.1)]

I1 = 2
I2 = 1
I3 = 2/3

##
a1 = (I2-I3)/I2/I3
a2 = (I3-I1)/I3/I1
a3 = (I1-I2)/I1/I2
f(y) = [a1*y[2]*y[3] a2*y[3]*y[1] a3*y[1]*y[2]]
##
# f(y) = [y[1]/I1 y[2]/I2 y[3]/I3]
h = .001
tf = 11
H0 = 1/2*(y0[1]^2/I1+y0[2]^2/I2+y0[3]^2/I3)


include("MyForwardEuler.jl")
y_FE,t = MyForwardEuler.feuler(f,tf,h,y0)
H_FE = 1/2*(y_FE[:,1].^2/I1+y_FE[:,2].^2/I2+y_FE[:,3].^2/I3)
plot3d(y_FE[:,1],y_FE[:,2],y_FE[:,3])
    title!("Eulers Equations: Forward Euler")


include("MyImplicitMidpoint.jl")
y_IM, t = MyImplicitMidpoint.my_implicit_mid(f,tf,h,y0)
H_IM = 1/2*(y_IM[:,1].^2/I1+y_IM[:,2].^2/I2+y_IM[:,3].^2/I3)
plot3d(y_IM[:,1],y_IM[:,2],y_IM[:,3])
    title!("Eulers Equations: Implicit Midpoint")

include("GL_RK.jl")
y_GL,ti = GL_RK.gl_rk(f,tf,h,y0)
H_GL = 1/2*(y_GL[:,1].^2/I1+y_GL[:,2].^2/I2+y_GL[:,3].^2/I3)/H0
plot3d(y_GL[:,1],y_GL[:,2],y_GL[:,3])
    title!("Eulers Equations: Implicit Runge Kutta")

plot(t,H_FE/H0,label= "Forward Euler")
    plot!(t,H_IM/H0,label = "Implicit Midpoint")
    plot!(t,H_GL/H0,label = "Implicit Runge Kutta")
    ylims!((.9,1.6))
