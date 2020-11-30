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

using LinearAlgebra
using Calculus

delta=0.1

function fun(x::Array{Float64,1})
    f = x[1]^2-x[1]-2*x[2]-x[1]*x[2]+x[2]^2
end

function grad(x::Array{Float64,1})
    g = Calculus.gradient(x -> x[1]^2-x[1]-2*x[2]-x[1]*x[2]+x[2]^2,[x[1],x[2]])
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
function LS_GS(d::Array{Float64,1},a::Float64,b::Float64,x::Array{Float64,1})

       # Vamos garantir que o usuário não seja uma anta
       @assert a<b "LS_GS:: o intervalo inicial deve ser [a,b]"

       # Vamos definir a razao aurea
       aurea =  (sqrt(5.0)+1.0)/2.0
       tol=1E-5
       # Agora vamos para o LS mesmo - Razao Aurea
       I =  b-a

       alpha1 =  a + I/aurea
       alpha2 =  b - I/aurea

       while I>tol

            f1 =  fun(x .+ alpha1*d)
            f2 =  fun(x .+ alpha2*d)

           if f1<f2
                a =  alpha2
                b = b
           elseif f1>f2
             a =  a
             b =  alpha1
           else
              a =  alpha2
              b =  alpha1
           end
           I =  b-a
           alpha1 = a + I/aurea
           alpha2 = a + I/(aurea^2)
       end

       # temos o minimo
       println(" Minimo com tolerancia ",I)
       return (a+b)/2.0

end

function Stpst_dsct(x::Array{Float64}, n::Int64)
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
                alpha = LS_GS(d,a,b,x)
                x = x.+alpha*d

                end
        end
    println("Xotimo é", x)

    return x
    end

Stpst_dsct([0.0,0.0],300)
