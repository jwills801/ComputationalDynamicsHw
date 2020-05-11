using Plots, Printf

L1 = .5
L2 = .5
m1 = 1
m2 = 1
g = 9.81

f(x) = [x[3] x[4] (-L2.*m2.*(L1 + L2.*cos(x[1] - x[2])).*(-2*L1.*(x[3].^2).*sin(x[1] - x[2]) - 4*L1.*x[3].*x[4].*sin(x[1] - x[2]) + 2*L1.*(x[4].^2).*sin(x[1] - x[2]) + g.*sin(x[2])) + (L1.^2 + 2*L1.*L2.*cos(x[1] - x[2]) + 2*L2.^2).*(4*L2.*m2.*(x[4].^2).*sin(x[1] - x[2]) + g.*m1.*sin(x[1]) + g.*m2.*sin(x[1])))./(2*L1.*(m2.*(L1 + L2.*cos(x[1] - x[2])).^2 - (2*m1 + m2).*(L1.^2 + 2*L1.*L2.*cos(x[1] - x[2]) + 2*L2.^2)))                                                                                (L2.*(2*m1 + m2).*(-2*L1.*(x[3].^2).*sin(x[1] - x[2]) - 4*L1.*x[3].*x[4].*sin(x[1] - x[2]) + 2*L1.*(x[4].^2).*sin(x[1] - x[2]) + g.*sin(x[2])) - (L1 + L2.*cos(x[1] - x[2])).*(4*L2.*m2.*(x[4].^2).*sin(x[1] - x[2]) + g.*m1.*sin(x[1]) + g.*m2.*sin(x[1])))./(2*(m2.*(L1 + L2.*cos(x[1] - x[2])).^2 - (2*m1 + m2).*(L1.^2 + 2*L1.*L2.*cos(x[1] - x[2]) + 2*L2.^2)))]
f(x0)
x0 = [3π/4 π/4 1 -2]
tf = 20
include("Adaptive_time_step_Runge_Kutta.jl")
sol,t = Adaptive_time_step_Runge_Kutta.vRK(f,tf,x0)

plot(t,sol[:,2],label = "θ2")
        plot!(t,sol[:,1],label = "θ1")


pyplot()

θ1 = sol[:,1]
θ2 = sol[:,2]
x1 = L1*sin.(θ1)
y1 = -L1*cos.(θ1)
x2 = x1 + L2*sin.(θ2)
y2 = y1 - L2*cos.(θ2)



plt = plot( [0,x1[1]], [0,y1[1]] , lw=3, label="",
        aspect_ratio=:equal,
        xaxis="x [m]", yaxis="y [m]",
        xlim=(-L1-L2, +L1+L2),
        ylim=(-L1-L2, +L1+L2));
        plot!(plt, [x1[1], x2[1]], [y1[1], y2[1]], lw=3, label="Runge Kutta Variable Time Step");


anim = @animate for i = 1:100:length(sol[:,1])
    plt[1] = ( [0,x1[i]], [0,y1[i]] )
    plt[2] = ([x1[i], x2[i]], [y1[i], y2[i]] )
    plot!(plt, title=@sprintf("Time t = %6.3f [s]", t[i]) )
end
gif(anim, "RungeKutta.gif", fps=10)

## Now to animate the Adams_Bashforth2 pendulum

include("Adams_Bashforth2.jl")
sol_AB,t_AB = Adams_Bashforth2.ab2(f,tf,.0005,x0)


plot(t_AB,sol_AB[:,1],label = "θ1 AB")
        plot!(t,sol[:,1],label = "θ1 RK")


θ1_AB = sol_AB[:,1]
θ2_AB = sol_AB[:,2]
x1_AB = L1*sin.(θ1_AB)
y1_AB = -L1*cos.(θ1_AB)
x2_AB = x1_AB + L2*sin.(θ2_AB)
y2_AB = y1_AB - L2*cos.(θ2_AB)


plt_AB = plot( [0,x1_AB[1]], [0,y1_AB[1]] , lw=3, label="",
        aspect_ratio=:equal,
        xaxis="x [m]", yaxis="y [m]",
        xlim=(-L1-L2, +L1+L2),
        ylim=(-L1-L2, +L1+L2));
        plot!(plt_AB, [x1_AB[1], x2_AB[1]], [y1_AB[1], y2_AB[1]], lw=3, label="Adams Bashforth 2 Constant Time Step");


anim_AB = @animate for i = 1:100:length(x1_AB)
    plt_AB[1] = ( [0,x1_AB[i]], [0,y1_AB[i]] )
    plt_AB[2] = ([x1_AB[i], x2_AB[i]], [y1_AB[i], y2_AB[i]] )
    plot!(plt_AB, title=@sprintf("Time t = %6.3f [s]", t_AB[i]) )
end
gif(anim_AB, "Adams_Bashforth_2.gif", fps=10)

### Now to compare different initial conditions
# two different initial conditions
x01 = 0.001 .+ x0

sol_AB_plus,t_AB_plus = Adams_Bashforth2.ab2(f,tf,.0005,x01)



θ1_AB_plus = sol_AB_plus[:,1]
θ2_AB_plus = sol_AB_plus[:,2]
x1_AB_plus = L1*sin.(θ1_AB_plus)
y1_AB_plus = -L1*cos.(θ1_AB_plus)
x2_AB_plus = x1_AB_plus + L2*sin.(θ2_AB_plus)
y2_AB_plus = y1_AB_plus - L2*cos.(θ2_AB_plus)




plt_chaos = plot( [0,x1_AB_plus[1]], [0,y1_AB_plus[1]] , lw=3, label="",
        aspect_ratio=:equal,
        xaxis="x [m]", yaxis="y [m]",
        xlim=(-L1-L2, +L1+L2),
        ylim=(-L1-L2, +L1+L2));
        plot!(plt_chaos, [x1_AB_plus[1], x2_AB_plus[1]], [y1_AB_plus[1], y2_AB_plus[1]], lw=3, label="x0 + 0.01")
        plot!(plt_chaos, [0, x1_AB[1]], [0, y1_AB[1]], lw=3, label="")
        plot!(plt_chaos, [x1_AB[1], x2_AB[1]], [y1_AB[1], y2_AB[1]], lw=3, label="x0");


anim_chaos = @animate for i = 1:100:length(x1_AB_plus)
    plt_chaos[1] = ( [0,x1_AB_plus[i]], [0,y1_AB_plus[i]] )
    plt_chaos[2] = ([x1_AB_plus[i], x2_AB_plus[i]], [y1_AB_plus[i], y2_AB_plus[i]] )
            plt_chaos[3] = ( [0,x1_AB[i]], [0,y1_AB[i]] )
            plt_chaos[4] = ([x1_AB[i], x2_AB[i]], [y1_AB[i], y2_AB[i]] )
            plot!(plt_chaos, title=@sprintf("Time t = %6.3f [s]", t_AB[i]) )

    end
gif(anim_chaos, "Chaos.gif", fps=15)






## Now, Ill put both Adams and Runge on the same animation

plt_both = plot( [0,x1[1]], [0,y1[1]] , lw=3, label="",
        aspect_ratio=:equal,
        xaxis="x [m]", yaxis="y [m]",
        xlim=(-L1-L2, +L1+L2),
        ylim=(-L1-L2, +L1+L2));
        plot!(plt_both, [x1[1], x2[1]], [y1[1], y2[1]], lw=3, label="Runge Kutta Variable Time Step")
        plot!(plt_both, [0, x1_AB[1]], [0, y1_AB[1]], lw=3, label="")
        plot!(plt_both, [x1_AB[1], x2_AB[1]], [y1_AB[1], y2_AB[1]], lw=3, label="Adams Bashforth constant time step");


anim_both = @animate for i = 1:100:length(x1)
    plt_both[1] = ( [0,x1[i]], [0,y1[i]] )
    plt_both[2] = ([x1[i], x2[i]], [y1[i], y2[i]] )
            plt_both[3] = ( [0,x1_AB[i]], [0,y1_AB[i]] )
            plt_both[4] = ([x1_AB[i], x2_AB[i]], [y1_AB[i], y2_AB[i]] )
            plot!(plt_both, title=@sprintf("Time t_AB = %6.3f [s]      Time = t_RK = %6.3f [s] ", t_AB[i],t[i]) )

end
gif(anim_both, "Two_Integrators.gif", fps=10)
