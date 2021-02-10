


using LinearAlgebra

function fun(x::Array{Float64,1})
    f = x[1]^2-x[1]-2*x[2]-x[1]*x[2]+x[2]^2
end
delta = 1E-2
alpha = delta
f1 = fun(x .+ (alpha - delta)*d)
f2 = fun(x .+ alpha*d)
f3 = fun(x .+ (alpha + delta)*d)
d1 = (f3 - f2) /(2*delta)
d2 = (f3 - 2 * f2 + f1)/(delta^2)
alpha = alpha - ((d1)/d2)
