module BackwardEuler

using LinearAlgebra

export beuler

function beuler(f, tf, h, x0; tol = 1e-4, iterMax = 500 )

    # x(i+1) = x(i) + h*f(i+1)

    time = 0:h:tf
    n = length(time)

    x = x0

    for i = 1:n-1
        # x(i+1) = x(i) + h*f( x(i+1) )

        # Predict via forward Euler x[i+1]
        y = x[i,:] + h*f( x[i,1], x[i,2] )

        flag = 0
        iter = 0
        while flag == 0
            iter += 1
            y = x[i,:] + h*f(y[1], y[2])

            residual = norm( y - x[i,:] - h*f(y[1], y[2]))

            if residual <= tol
                flag = 1
            elseif iter >= iterMax
                flag = -1
                error("Error: failed to converge")
            end
        end

        x = vcat(x,y')

    end

    return x, time

end

end
