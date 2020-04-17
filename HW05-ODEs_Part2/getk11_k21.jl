
function getk11_k21(eq3,eq4)
    k11_0 = 0
    k21_0 = 0
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
