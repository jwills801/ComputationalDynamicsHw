module Adams_Bashforth2

using LinearAlgebra

export ab2

function ab2(f, tf, h, x0)
    time = 0:h:tf
    n = length(time)


    # Initial Condition
    x = x0

    # 1 Euler Step (this is cheap and possibly bad)
    fn = f(x)
    x = vcat(x0, x0 + h*fn)

    # AB2 the rest
    for i = 2:n-1
        fn_m1 = fn # f(x[i-1,:][1], x[i-1,:][2])
        fn = f(x[i,:])
        #x[i+1,:] = x[i,:] + h/2*( 3*f(x[i,:], time[i]) - f(x[i-1,:], time[i-1])  )
        xnext = x[i,:]' + h/2*( 3*fn - fn_m1 )
        x = vcat(x,xnext)
    end

    return x, time
end



end
