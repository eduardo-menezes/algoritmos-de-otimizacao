# Rotina para a determinação do intervalo inicial. Dado o ponto atual x0, uma direção de busca d,
# uma estimativa inicial de passo delta, uma função objetivo fun e uma tolerância, devolve o
# intervalo que deve conter o alpha ótimo.

#
# fun -> função com entrada x::Array{Float64,1}
# x0 -> ponto atual
# d -> direção de minimização
# delta -> intervalo inicial
using LinearAlgebra
using Calculus

delta=0.1

function fun(x::Array{Float64,1})
    f = (x[1]-3)^2+(x[2]-2)^2
end

function grad(x::Array{Float64,1})
    g = Calculus.gradient(x -> (x[1]-3)^2+(x[2]-2)^2,[x[1],x[2]])
end


function Armijo(x::Array{Float64,1}, d::Array{Float64,1}, c1::Float64, reduc::Float64,n::Int64)
    alpha = 1
    df = grad(x)
    norma = norm(df)
    dir = df/norma
    f0 = fun(x) + c1*alpha*dot(dir,d)
    f1 = fun(x .+ alpha*dir) + c1*alpha*dot(dir,d)
        while f1 > f0
            alpha = alpha/reduc
            f0 = fun(x) + c1*alpha*dot(dir,d)
            f1 = fun(x .+ alpha*d) + c1*alpha*dot(dir,d)
        end
        alphaotm = alpha
        println("alpha otimo é: ", alphaotm)
        return alphaotm
end

function steepest_Armijo(x::Array{Float64},c1::Float64, reduc::Float64, n::Int64)
        tol = 1E-4

        for i= 1:n
            gradien = grad(x)
            norma = norm(gradien)

                if norma <= tol
                    println("X otimo é: ", x)
                    break

                else

                d = -gradien/norma
                alpha = Armijo(x,d,c1, reduc,n)
                x = x.+alpha*d

                end
        end
    println("Xotimo é", x)

    return x
    end

steepest_Armijo([1.0,1.0],0.1,5.0,100)
