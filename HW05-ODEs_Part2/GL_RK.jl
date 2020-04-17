module GL_RK

using LinearAlgebra

export gl_rk

function gl_rk(f,tf,h,x0)
  time = 0:h:tf
  n = length(time)

  x = x0


  c1 = 1/2 - sqrt(3)/6
  c2 = 1/2 + sqrt(3)/6
  a11 = 1/4
  a12 = 1/4 - sqrt(3)/6
  a21 = 1/4 + sqrt(3)/6
  a22 = 1/4
  b1 = 1/2
  b2 = 1/2

  k1_guess = [0 0.]
  k2_guess = [0 0.]
  for i=1:n-1
    for j=1:500
      k1 = f(x[i,:]' + h*a11*k1_guess + h*a12*k2_guess)
      k2 = f(x[i,:]' + h*a21*k1 + h*a22*k2_guess)
      if norm(k1_guess - k1) <= .0001
        if norm(k2_guess - k2) <= .0001
          break
        end
      end
      k1_guess = k1
      k2_guess = k2

    end


    xnext = x[i,:]' + h*(b1*k1_guess + b2*k2_guess)
    x = vcat(x,xnext)
  end



 return x,time

end # function gl_rk

end  # module GL_RK
