
import Pkg

using LinearAlgebra
using Calculus

function fun(x::Array{Float64,1})
    f = x[1]^2-x[1]-2*x[2]-x[1]*x[2]+x[2]^2
end

function grad(x::Array{Float64,1}) #Testar DF no lugar do gradien
    g = Calculus.gradient(x -> x[1]^2-x[1]-2*x[2]-x[1]*x[2]+x[2]^2,[x[1],x[2]])
end

function NewtonMod(x::Array{Float64,1},d::Array{Float64,1},tol::Float64,n::Int64)
    alpha = 0
       for i=1:n
           delta = 1E-2
           alpha = delta
           f1 = fun(x .+ (alpha - delta)*d)
           f2 = fun(x .+ alpha*d)
           f3 = fun(x .+ (alpha + delta)*d)
           d1 = (f3 - f2) /(2*delta)
           d2 = (f3 - 2 * f2 + f1)/(delta^2)
           alpha = alpha - ((d1)/d2)
       end

       # temos o minimo
       println(" Minimo igual a ",alpha)
       return alpha

end

function Steepest_NewtonMod(x::Array{Float64}, n::Int64, tol::Float64, delta::Float64)

        for i= 1:n
            gradien = grad(x)
            norma = norm(gradien)

                if norma <= tol
                    println("X otimo é: ", x)
                    break

                else

                d = -gradien/norma
                alpha = NewtonMod(x,d,tol,n)
                x = x.+alpha*d

                end
        end
    println("Xotimo é", x)

    return x
    end

Steepest_NewtonMod([0.0,0.0],10, 1E-6, 1E-2)
