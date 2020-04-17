module MyForwardEuler

export feuler

function feuler(f,tf,h,x0)
    x = x0
    time = 0:h:tf
    n = length(time)
    for i=1:n-1
        xnext = x[i,:]' + h*f(x[i,:])
        x = vcat(x,xnext)
    end
    return x,time


end # function feuler


end  # modul MyForwardEuler
