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
import Pkg
using LinearAlgebra
using Calculus


function fun(x::Array{Float64,1})
    f = (x[1]-3)^2+(x[2]-2)^2
end

function grad(x::Array{Float64,1})
    g = Calculus.gradient(x -> (x[1]-3)^2+(x[2]-2)^2,[x[1],x[2]])
end

function Intervalo_Inicial(x::Array{Float64,1},d::Array{Float64})

       # Vamos definir a razao aurea
       aurea =  (sqrt(5.0)+1.0)/2.0
       delta = 0.1
       # Primeiro temos que achar o ]intervalo inicial
       alpha =  0.0

       # Valor de referencia
       xt =  x .+ alpha*d
       f0 = fun(xt)

       # Se f(alpha=delta) ja for maior do que f0, podemos
       #    utilizar este intervalo diretamente. Se f1<f0,
       #    temos que fazer achar o interalo inicial
       xt =  x .+ delta*d
       f1 =  fun(xt)

       # Padrao eh a=0 e b=delta. So muda se entrarmos no if abaixo
       a = 0.0
       b = delta

       if f1<f0
          n =  1
          while f1<f0
              n = n+1
              alpha =  (aurea^n)*delta
              xt = x .+ alpha*d
              f1 = fun(xt)
         end # while
           b = (aurea^(n-1))*delta
           a = (aurea^(n-3))*delta
       end # if f1<f2*/

       println("Bracketing::Limites do intervalo: ",a," ",b)
       return a,b

end

# Rotina de Line Search - Golden Search
# Acha o passo para o mínimo de fun(x) na direção d
#
# Eduardo Lenz - Setembro de 2018
function BuscaExaust(d::Array{Float64,1},a::Float64,b::Float64,x::Array{Float64,1},n::Int64)

    delta = 0.1
    FR = 10
    alpha0 = 0
    x = x .+ alpha0 * d
    f0 = fun(x)
    alpha1 = delta
    x1 = x .+ alpha1 * d
    f1 = fun(x1)
    while f1 < f0
        alpha0 = alpha1
        f0 = f1
        alpha1 = alpha0 + delta
        x1 = x .+ alpha1 * d
        f1 = fun(x1)
    end
    delta = delta/FR
    alpha0 = alpha1*delta
    x = x .+ alpha0 * d
    f0 = fun(x)
    alpha1 = alpha0 + delta
    x1 = x .+ alpha1 * d
    f1 = fun(x1)
        while f1 < f0
            alpha0 = alpha1
            f0 = f1
            alpha1 = alpha0 + delta
            x1 = x .+ alpha1 * d
            f1 = fun(x1)
            end
        println("Ponto Minimo é: ", (alpha1+alpha0)/2)
        return (alpha1+alpha0)/2
    end

    function Steepest_buscaExaust(x::Array{Float64}, n::Int64)
            tol = 1E-4

            for i= 1:n
                gradien = grad(x)
                norma = norm(gradien)

                    if norma <= tol
                        println("X otimo é: ", x)
                        break

                    else

                    d = -gradien/norma
                    a,b = Intervalo_Inicial(x,d)
                    alpha = BuscaExaust(d,a,b,x,n)
                    x = x.+alpha*d

                    end
            end
        println("Xotimo é", x)

        return x
        end

Steepest_buscaExaust([0.0,0.0],100)
