module Adaptive_time_step_Runge_Kutta

using LinearAlgebra

export vRk

function vRK(f,tf,x0;h0 = 1,tol = 1e-8)

    c1 = 1/2
    c2 = 3/4
    c3 = 1

    a1 = 1/2
    a21 = 0
    a22 = 3/4
    a31 = 2/9
    a32 = 1/3
    a33 = 4/9

    b11 = 2/9
    b12 = 1/3
    b13 = 4/9
    b14 = 0
    b21 = 7/24
    b22 = 1/4
    b23 = 1/3
    b24 = 1/8

    t = [0]
    x = x0
    h = h0
    flag = 0
    i = 0
    error_m1 = 100
    while flag == 0
        i += 1
        for j= 1:50
            k1 = f(t[i],x[i,:]')
            k2 = f(t[i]+c1*h,x[i,:]'+h*k1*a1)
            k3 = f(t[i] + c2*h, x[i,:]' + h*k1*a21 + h*k2*a22)

            xnext = x[i,:]' + h*(b11*k1 + b12*k2 + b13*k3)

            k4 = f(t[i]+c3*h,xnext)

            znext = x[i,:]' + h*(b21*k1 + b22*k2 + b23*k3 + b24*k4)

            error = norm(xnext-znext)


            if error <= tol
                t = vcat(t,t.+h)
                h = .9*h*min(max(error/error_m1,.3),2)
                x = vcat(x,xnext)
                break
            else
                h = h/2
            end #if error
            if j == 49
                println("did not converge")
                break
            end

            error_m1 = error
        end # j for loop

        if t[i] >= tf
            flag = 1
        end
        if i >= 100
            flag = -1
        end



    end # while loop

    return t, x

end #function vRK

end  # moduleAdaptive_time_step_Runge_Kutta
