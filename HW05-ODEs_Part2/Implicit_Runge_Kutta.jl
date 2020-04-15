
using SymPy
@vars p q theta

Ham = p^2/2 +cos(q)

# Assuming the p^2/2 is the Kinetic Energy Term
# And that cos(q) is the potential energy term...

T = p^2/2
V = cos(q)
L = T-V


# now I'll find the generalized momentum of the angle
pdot = -diff(Ham,q) |> subs(q=>theta)
thetadot = diff(Ham,p)

# therefore thetaddot = sin(theta)
