using LinearAlgebra
using SymPy
include("getk11_k21.jl")

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
        @vars k11 k12 k21 k22
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
        k12_SymVal = diction[k12]
        k22_SymVal = diction[k22]

        eq3 = eqns1_2[1]
        eq4 = eqns3_4[1]


        k11_SymVal, k21_SymVal = getk11_k21(eq3,eq4)
        k11 = convert(Float64,k11_SymVal)
        k21 = convert(Float64,k21_SymVal)
        k12 = convert(Float64,k12_SymVal)
        k22 = convert(Float64,k22_SymVal)

        k1 = [k11,k12]
        k2 = [k21,k22]


        xnext = x[i,:] + h*(b1*k1 + b2*k2)
        x = vcat(x,xnext')

    end

    return x, time

end
