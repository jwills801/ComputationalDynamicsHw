module MyImplicitMidpoint

using LinearAlgebra

function my_implicit_mid(f,tf,h,x0)
    x = x0
    time = 0:h:tf
    n = length(time)

    for i=1:n-1
    xnext_guess = x[1,:]'

        for j = 1:500
            xnext = x[i,:]' + h*f(1/2*(x[i,:]' + xnext_guess))
            if norm(xnext-xnext_guess) <= .001
                x = vcat(x,xnext)
                break
            end
            xnext_guess = xnext
        end # j for loop
    end #i for loop

    return x,time


end # end my_implicit_mid

end  # moduleMyImplicitMidpoint
