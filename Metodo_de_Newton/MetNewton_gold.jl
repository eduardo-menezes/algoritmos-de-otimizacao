# Rotina para a determinação do intervalo inicial. Dado o ponto atual x0, uma direção de busca d,
# uma estimativa inicial de passo delta, uma função objetivo fun e uma tolerância, devolve o
# intervalo que deve conter o alpha ótimo.
#

#
# fun -> função com entrada x::Array{Float64,1}
# x0 -> ponto atual
# d -> direção de minimização
# delta -> intervalo inicial


using LinearAlgebra
using Calculus

delta=0.1

function fun(x::Array{Float64,1})
    f = x[1]^2-x[1]-2*x[2]-x[1]*x[2]+x[2]^2
end

function grad(x::Array{Float64,1})
    g = Calculus.gradient(x -> x[1]^2-x[1]-2*x[2]-x[1]*x[2]+x[2]^2,[x[1],x[2]])
end

function Hess(x::Array{Float64,1})
    H = Calculus.hessian(x -> x[1]^2-x[1]-2*x[2]-x[1]*x[2]+x[2]^2,[x[1],x[2]])

end

function MetNewtonMod(x::Array{Float64}, n::Int64)
        tol = 1E-4
        iter = 0
        for i= 1:n
            gradien = grad(x)
            @show gradien
            norma = norm(gradien)
            @show norma
                if norma <= tol
                    println("X otimo é: ", x)
                    break

                else
                H = Hess(x)
                autovalor = eigvals(H)
                detH = det(H)
                inversaH = inv(H)
                @show inversaH
                @show det(H)
                @show H
                @show autovalor
                d =  -H\gradien
                @show d
                x = x + d
                @show x
                end
            iter = iter + 1
            @show iter
        end

    return x
    end

MetNewtonMod([0.0,0.0],300)
