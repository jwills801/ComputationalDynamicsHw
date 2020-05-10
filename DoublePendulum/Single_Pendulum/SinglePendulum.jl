using Plots

g= 9.81
l =1
fp1(x) = [x[2] -g/l*sin(x[1])]

x0 = [π/2 0]

include("../Adaptive_time_step_Runge_Kutta.jl")
sol,t = Adaptive_time_step_Runge_Kutta.vRK(fp1,10,x0)

pyplot()

plot(t,sol[:,1])

θ = sol[:,1]
x = sin.(θ)
y = -cos.(θ)

plt = plot( [0,x[i]], [0,y[i]] , lw=3, label="",
        aspect_ratio=:equal,
        xaxis="x [m]", yaxis="y [m]",
        xlim=(-l, +l),
        ylim=(-l, +l));


anim = @animate for i = 1:100:length(sol[:,1])
    plt[1] = ( [0,x[i]], [0,y[i]] )
    plot!(plt, title=@sprintf("Time t = %6.3f [s]", t[i]) )
end
gif(anim, "test1.gif", fps=10)
