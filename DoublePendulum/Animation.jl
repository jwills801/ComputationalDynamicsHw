module Animation
using Plots, Printf
function animate_double_pen(f,tf,x0; h = .0005, title = "Adams Bashforth Constant Time Step of 0.0005 seconds")

    include("Adams_Bashforth2.jl")
    sol,t = Adams_Bashforth2.ab2(f,tf,x0)
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
    gif(anim, "test1.gif", fps=10)
end

end  # module Animation
