# Rotina para a determinação do intervalo inicial. Dado o ponto atual x0, uma direção de busca d,
# uma estimativa inicial de passo delta, uma função objetivo fun e uma tolerância, devolve o
# intervalo que deve conter o alpha ótimo.
#
# Eduardo Lenz - setembro de 2018
#
# fun -> função com entrada x::Array{Float64,1}
# x0 -> ponto atual
# d -> direção de minimização
# delta -> intervalo inicial
clearconsole()

using LinearAlgebra
using Calculus

delta=0.1

function fun(x::Array{Float64,1})
    f = (x[1]-3)^2+(x[2]-2)^2
end

function grad(x::Array{Float64,1})

    g =Calculus.gradient(x -> (x[1]-3)^2+(x[2]-2)^2,[x[1],x[2]])
end

function InterPolinom(x::Array{Float64,1}, d::Array{Float64,1}, tol::Float64, n::Int64)
        alpha = 0
        delta = 1E-2
        alpha = delta
        f1 = fun(x .+ (alpha - delta)*d)
        f2 = fun(x .+ alpha*d)
        f3 = fun(x .+ (alpha + delta)*d)
        a0 = ((f3 - 2*f2 + f1)*alpha^2 + (delta*f1 - delta*f3)*alpha + 2*delta^2*f2)/(2*delta^2)
        a1 = -((2*f3 - 4*f2 + 2*f1)*alpha - delta*f3 + delta*f1)/(2*delta^2)
        a2 = (f3 - 2*f2 + f1)/(2*delta^2)
        alpha =  ((2*f3-4*f2+2*f1)*alpha-delta*f3+delta*f1)/(2*(f3-2*f2+f1))
        d = 2*a2*alpha+a1
            while d>tol
                f1 = fun(x .+ (alpha - delta)*d)
                f2 = fun(x .+ alpha*d)
                f3 = fun(x .+ (alpha + delta)*d)
                a0 = ((f3 - 2*f2 + f1)*alpha^2 + (delta*f1 - delta*f3)*alpha + 2*delta^2*f2)/(2*delta^2)
                a1 = -((2*f3 - 4*f2 + 2*f1)*alpha - delta*f3 + delta*f1)/(2*delta^2)
                a2 = (f3 - 2*f2 + f1)/(2*delta^2)
                alpha =  ((2*f3-4*f2+2*f1)*alpha-delta*f3+delta*f1)/(2*(f3-2*f2+f1))
                d = 2*a2*alpha+a1
            end
       # temos o minimo
       println(" Minimo igual a ",alpha)
       return alpha
   end

function Steepest_InterPolinom(x::Array{Float64}, n::Int64, tol::Float64, delta::Float64)
        tol = 1E-4

        for i= 1:n
            gradien = grad(x)
            norma = norm(gradien)

                if norma <= tol
                    println("X otimo é: ", x)
                    break

                else

                d = -gradien/norma
                alpha = InterPolinom(x,d,tol,n)
                x = x.+alpha*d

                end
        end
    println("Xotimo é", x)

    return x
    end

Steepest_InterPolinom([0.0,0.0],10, 1E-6, 1E-2)
