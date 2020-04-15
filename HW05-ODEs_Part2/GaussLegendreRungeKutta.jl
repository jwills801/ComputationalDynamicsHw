module GaussLegendreRungeKutta

using LinearAlgebra
using SymPy

export IRK


function IRK(f, tf, h, x0; tol = 1e-4, iterMax = 500 )

    # x(i+1) = x(i) + h*f(i+1)

    time = 0:h:tf
    n = length(time)

    x = x0
    x_float = x0

    c1 = 1/2 - sqrt(3)/6
    c2 = 1/2 + sqrt(3)/6
    a11 = 1/4
    a12 = 1/4 - sqrt(3)/6
    a21 = 1/4 + sqrt(3)/6
    a22 = 1/4
    b1 = 1/2
    b2 = 1/2


    # k1 = f(x[1,1]+h*a11*k1 + h*a12*k2,x[1,2]+h*a11*k1 + h*a12*k2)
    # f(x1,x2) = [sin(x2);x1]
    # k1 =[sin(x[1,2]+h*a11*k1 + h*a12*k2);x[1,1]+h*a11*k1 + h*a12*k2]
    for i = 1:n-1
        k11 = Symbol("k11")
        k12 = Symbol("k12")
        k21 = Symbol("k21")
        k22 = Symbol("k22")
        k1 = [k11,k12]
        k2 = [k21,k22]
        eqns1_2 = k1 - [sin(x[i,2]+h*a11*k11 + h*a12*k21);x[i,1]+h*a11*k12 + h*a12*k22]

        # k2 = f(x[1,1]+h*a21*k1 + h*a22*k2,x[1,2]+h*a21*k1 + h*a22*k2)
        # f(x1,x2) = [sin(x2);x1]
        # k2 =[sin(x[1,2]+h*a21*k1 + h*a22*k2);x[1,1]+h*a21*k1 + h*a22*k2]

        eqns3_4 = k2 - [sin(x[i,2]+h*a21*k11 + h*a22*k21);x[i,1]+h*a21*k12 + h*a22*k22]

        eq1 = eqns1_2[2]
        eq2 = eqns3_4[2]
        diction = solve([eq1, eq2], [k12,k22])
        k12 = diction[k12]
        k22 = diction[k22]

        eq3 = eqns1_2[1]
        eq4 = eqns3_4[1]

        k11_0 = 0
        k21_0 = 0
        k11, k21 = getk11_k21(k11_0,k21_0)



        k1 = [k11,k12]
        k2 = [k21,k22]


        xnext = x[i,:] + h*(b1*k1 + b2*k2)
        x = vcat(x,xnext')

    end

    return x, time

end




function getk11_k21(k11_0,k21_0)
    for i = 1:500
        IsItZero1 = eq3 |> subs(k11=>k11_0,k21=>k21_0)
        IsItZero2 = eq4 |> subs(k11=>k11_0,k21=>k21_0)


        if IsItZero1 >= .001
            break
        end
        if IsItZero2 >= .001
            break
        end

        k11_0 = k11_0 + .001
        k21_0 = k21_0 + .001


    end
    return k11_0, k21_0

end



end


















i=2
@vars k11 k12 k21 k22
k1 = [k11,k12]
k2 = [k21,k22]
eqns1_2 = k1 - [sin(x[1,2]+h*a11*k11 + h*a12*k21);x[1,1]+h*a11*k12 + h*a12*k22]

    # k2 = f(x[1,1]+h*a21*k1 + h*a22*k2,x[1,2]+h*a21*k1 + h*a22*k2)
    # f(x1,x2) = [sin(x2);x1]
    # k2 =[sin(x[1,2]+h*a21*k1 + h*a22*k2);x[1,1]+h*a21*k1 + h*a22*k2]

eqns3_4 = k2 - [sin(x[i,2]+h*a21*k11 + h*a22*k21);x[i,1]+h*a21*k12 + h*a22*k22]

eq1 = eqns1_2[2]
eq2 = eqns3_4[2]
diction = solve([eq1, eq2], [k12,k22])
k12 = diction[k12]
k22 = diction[k22]

eq3 = eqns1_2[1]
eq4 = eqns3_4[1]

k11_0 = 0
k21_0 = 0
k11, k21 = getk11_k21(k11_0,k21_0)



k1 = [k11,k12]
k2 = [k21,k22]


xnext = x[i,:] + h*(b1*k1 + b2*k2)
x = vcat(x,xnext')


x_float = x0

x_float = vcat(x_float,[x[i,:][1] )

[x[i,:][1],x[i,:][2]]
