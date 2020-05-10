using Plots, Printf

L1 = .5
L2 = .5
m1 = 1
m2 = 1
g = 9.81

f(x) = [x[3] x[4] (-L2.*m2.*(L1 + L2.*cos(x[1] - x[2])).*(-2*L1.*(x[3].^2).*sin(x[1] - x[2]) - 4*L1.*x[3].*x[4].*sin(x[1] - x[2]) + 2*L1.*(x[4].^2).*sin(x[1] - x[2]) + g.*sin(x[2])) + (L1.^2 + 2*L1.*L2.*cos(x[1] - x[2]) + 2*L2.^2).*(4*L2.*m2.*(x[4].^2).*sin(x[1] - x[2]) + g.*m1.*sin(x[1]) + g.*m2.*sin(x[1])))./(2*L1.*(m2.*(L1 + L2.*cos(x[1] - x[2])).^2 - (2*m1 + m2).*(L1.^2 + 2*L1.*L2.*cos(x[1] - x[2]) + 2*L2.^2)))                                                                                (L2.*(2*m1 + m2).*(-2*L1.*(x[3].^2).*sin(x[1] - x[2]) - 4*L1.*x[3].*x[4].*sin(x[1] - x[2]) + 2*L1.*(x[4].^2).*sin(x[1] - x[2]) + g.*sin(x[2])) - (L1 + L2.*cos(x[1] - x[2])).*(4*L2.*m2.*(x[4].^2).*sin(x[1] - x[2]) + g.*m1.*sin(x[1]) + g.*m2.*sin(x[1])))./(2*(m2.*(L1 + L2.*cos(x[1] - x[2])).^2 - (2*m1 + m2).*(L1.^2 + 2*L1.*L2.*cos(x[1] - x[2]) + 2*L2.^2)))]
f(x0)
x0 = [3π/4 π/4 1 -2]
include("Adaptive_time_step_Runge_Kutta.jl")
sol,t = Adaptive_time_step_Runge_Kutta.vRK(f,20,x0)

plot(t,sol[:,2],label = "θ2")
        plot!(t,sol[:,1],label = "θ1")


pyplot()

θ1 = sol[:,1]
θ2 = sol[:,2]
x1 = L1*sin.(θ1)
y1 = -L1*cos.(θ1)
x2 = x1 + L2*sin.(θ2)
y2 = y1 - L2*cos.(θ2)

plt = plot( [0,x1[i]], [0,y1[i]] , lw=3, label="",
        aspect_ratio=:equal,
        xaxis="x [m]", yaxis="y [m]",
        xlim=(-L1-L2, +L1+L2),
        ylim=(-L1-L2, +L1+L2));
        plot!(plt, [x1[i], x2[i]], [y1[i], y2[i]], lw=3, label="");


anim = @animate for i = 1:100:length(sol[:,1])
    plt[1] = ( [0,x1[i]], [0,y1[i]] )
    plt[2] = ([x1[i], x2[i]], [y1[i], y2[i]] )
    plot!(plt, title=@sprintf("Time t = %6.3f [s]", t[i]) )
end
gif(anim, "test1.gif", fps=10)
