clearconsole()
# Soluciona Min c'x + 0.5 x'Hx
#      T.q  Ax <= b
# com c nx1
#     H nxn
#     A mxn
#     b mx1
include("golden.jl")
#
# retorna as variáveis primais e as duais separadas e nesta ordem
#
# Esta rotina tem função puramente educacional, pois não tem nenhuma
# preocupação com eficiência (muito pelo contrário!!!)
#
# uso
#
# x,mu = QP(A,b,c,H,tol)
#
# Eduardo Lenz - outubro de 2018.
#

using LinearAlgebra


# Rotinas necessárias para a solução do problema de maximização dual

#
# Mapeamento Primal-Dual
# Retorna x para um dado mu
#
function Map(mu::Array{Float64,1},c::Array{Float64,1},H::Array{Float64,2},
             A::Array{Float64,2})

      -H\(transpose(A)*mu  + c)


end

#
# Derivada total de x em relação a mu_m
#
function dMap(m::Int64,H::Array{Float64,2},A::Array{Float64,2})

        # linha m de A [1 x n]
        lma = A[m,:]

        # derivada de x (map) em relação a mu_m [n x 1]
        return -H\(lma)

end

#
# Valor da função dual
#
function dual(mu::Array{Float64,1},c::Array{Float64,1},H::Array{Float64,2},
             A::Array{Float64,2},b::Array{Float64,1})

     # Mapping
     xm = Map(mu,c,H,A)

     # valor do dual
     dot(c,xm) + 0.5*dot(xm,H*xm) + dot(mu, A*xm-b)

end


# Derivada total do dual em relação a um mu
function dadual( mu::Array{Float64,1},c::Array{Float64,1},H::Array{Float64,2}, A::Array{Float64,2},b::Array{Float64,1})

   # Mapeamento PD
   xmap = Map(mu,c,H,A)

   # Aloca o vetor de saída
   m = size(A,1)

   # Primeiro termo é a derivada parcial de l em relação a mu
   dd = A*xmap - b

   # Derivada de parcial de l em relação a x
   dx = c + H*xmap + transpose(A)*mu

   # Loop pelas componentes da derivada total
   for i=1:m

       # Derivada total de xmap em relação a uma posição i de mu
       dx_mu_i = dMap(i,H,A)
       dd[i] += dot(dx,dx_mu_i)

   end

   return dd

end


# Derivada total do dual em relação a um mu
function dadual2( mu::Array{Float64,1},c::Array{Float64,1},H::Array{Float64,2}, A::Array{Float64,2},b::Array{Float64,1})

   # Mapeamento PD
   xmap = Map(mu,c,H,A)

   # Aloca o vetor de saída
   m = size(A,1)

   # Primeiro termo é a derivada parcial de l em relação a mu
   dldmu = A*xmap - b

   # Derivada de parcial de l em relação a x
   dldx = c + H*xmap + transpose(A)*mu

   # Derivada total de x em relação a mu
   dxdmu = -H\transpose(A)

   return dldmu + transpose(dxdmu)*dldx

end



#
# Estabelece o intervalo inicial para o LS (mas para um problema de maximização)
# e considera eventual restrição de passo (devido a restrição de que mu>=0)
#
function Intervalo_Inicial_blq(A,b,c,H,mu,d::Array{Float64,1},delta::Float64=0.1)


       # Vamos calcular o valor do passo limite
       passo_limite = 255E255
       for j in LinearIndices(mu)

           dd = d[j]

           if dd<0.0
               passo_limite = min(passo_limite,(0.0-mu[j])/dd)
           end
       end #j

       #println("Passo limite ",passo_limite)

       # Vamos definir a razao aurea
       aurea =  (sqrt(5.0)+1.0)/2.0

       # Primeiro temos que achar o intervalo inicial
       alpha =  0.0

       # Valor de referencia
       mut =  mu .+ alpha*d
       f0 = -dual(mut,c,H,A,b)

       # Se f(alpha=delta) ja for maior do que f0, podemos
       #    utilizar este intervalo diretamente. Se f1<f0,
       #    temos que fazer achar o interalo inicial
       mut =  mu .+ delta*d
       f1 =  -dual(mut,c,H,A,b)

       # Padrao eh a=0 e b=delta. So muda se entrarmos no if abaixo
       aa = 0.0
       bb = delta

       if f1<f0
          n =  1
          while f1<f0
              n = n+1
              alpha =  (aurea^n)*delta
              mut = mu .+ alpha*d
              f1 = -dual(mut,c,H,A,b)
              #println(f0," ",f1)
          end # while
           bb = (aurea^(n-1))*delta
           aa = (aurea^(n-3))*delta

           if bb>=passo_limite
              println("Intervalo inicial...não podemos violar o passo limite ", passo_limite)
              return 0.0,passo_limite
          end

       end # if f1<f2*/

       #println("Bracketing::Limites do intervalo: ",aa," ",bb)
       return aa,bb
end

#
# Rotina de Line Search - Golden Search
# Acha o passo para o máximo de fun(x) na direção d
#
function LS_MAX_GS_blq(A,b,c,H,mu,d::Array{Float64,1},
                       aa::Float64,bb::Float64, tol::Float64=1E-8)

       # Vamos garantir que o usuário não seja uma anta
       @assert aa<bb "LS_GS:: o intervalo inicial deve ser [a,b]"

              # Vamos calcular o valor do passo limite
              passo_limite = 255E255
              for j in LinearIndices(mu)

                  dd = d[j]

                  if dd<0.0
                      passo_limite = min(passo_limite,(0.0-mu[j])/dd)
                  end
              end #j

       #println("GS::passo limite é ", passo_limite)


       # Vamos definir a razao aurea
       aurea =  (sqrt(5)+1)/2

       # Agora vamos para o LS mesmo - Razao Aurea
       I =  bb-aa
       alpha1 =  aa + I/aurea
       alpha2 =  bb - I/aurea

       while (I>tol)

            f1 =  -dual(mu + alpha1*d ,c,H,A,b)
            f2 =  -dual(mu + alpha2*d ,c,H,A,b)

           if f1<f2
                aa =  alpha2
                bb = bb
           elseif f1>f2
             aa =  aa
             bb =  alpha1
           else
              aa =  alpha2
              bb =  alpha1
           end
           I =  bb-aa
           alpha1 = aa + I/aurea
           alpha2 = aa + I/(aurea^2)

           # Vamos verificar se não passamos do passo limite
           if (aa+bb)/2.0 > passo_limite
               return passo_limite
           end

       end

       # temos o minimo
       #println(" Maximo com tolerancia ",I)
       return (aa+bb)/2.0

end


#
# Soluciona o problema de programação quadrática
#
function QP(A,b,c,H,tol::Float64=1E-5, niter::Int64=100)


    # Vamos descobrir a dimensão do problema
    m = size(A,1)
    @show m

    # Começamos com mu=0
    mu = zeros(m)
    alpha = 0
    # Loop do Steepest Ascent
    for iter=1:niter

         # Gradiente
         dd = dadual(mu,c,H,A,b)

         # Bloqueio direto no gradiente. Aqui só temos que nos preocupar
         # com a situação onde o dd for negativo e mu for nulo
         for j=1:m
             if dd[j]<0.0 && mu[j]==0.0
                dd[j] = 0.0
             end
         end #j

         # Norma do gradiente
         norma = norm(dd)

         # Verifica pela convergência
         if norma <= tol
             println("Solução dual obtida::",mu)
             return  Map(mu,c,H,A),  mu
         end

         # Direção de MAXIMIZAÇÂO normalizada
         d = (dd) ./ norma

         println("Direção de maximização do dual ", d)

         # Line search
         aa,bb = Intervalo_Inicial_blq(A,b,c,H,mu,d)
         alpha = LS_MAX_GS_blq(A,b,c,H,mu,d,aa,bb)
         @show alpha

         # Evita que o passo seja muito pequeno
         if alpha<=1E-10
            break
         end

         # Novo ponto
         mu += alpha*d
         @show mu

  end

  # Retorna as variáveis primais e as duais
  xm = Map(mu,c,H,A)
  println("xm",xm)
  return xm, mu

end


#
# Problemas de verificação
#

#
# Exercício visto em sala: Solução em x = [2 ; -4] e mu = [0 ; 5]
#
function Test(x)

    function fun(x)
        ff = x[1]^2+x[2]^2-3*x[1]*x[2]
    end

    function drestricao(x)
        g11 = (1/3)*x[1]
        g12 = (1/3)*x[2]
        g21 = -1
        g22 = 0
        g31 = 0
        g32 = -1

        A = [g11 g12; g21 g22; g31 g32]
        return A
    end

    function Vetorb(x)
        g1 = (1/6)*x[1]^2+(1/6)*x[2]^2-1
        g2 = -x[1]
        g3= -x[2]

        g11 = (1/3)*x[1]
        g12 = (1/3)*x[2]
        g21 = -1
        g22 = 0
        g31 = 0
        g32 = -1


        b=[-g1+ g11*x[1]+g12*x[2];-g2+ g21*x[1]+g22*x[2];-g3+ g31*x[1]+g32*x[2]]
        return b
    end

    function Vetorc(x)
        c = [2*x[1]-3*x[2];2*x[2]-3x[1]]
        return c
    end

    function restricao(x)
        g1 = (1/6)*x[1]^2+(1/6)*x[2]^2-1
        g2 = -x[1]
        g3= -x[2]

        v = max(0,g1,g2,g3)

    end

        # Exercício visto em sala: Solução em x = [2 ; -4] e mu = [0 ; 5]
    tol = 1E-5
    H = [1.0 0.0; 0.0 1.0]
    @show H


    for i=1:5

    A = drestricao(x)
    @show A
    b0=Vetorb(x)
    @show b0
    c0 = Vetorc(x)
    @show c0

   x0, mu = QP(A,b0,c0,H,1E-5)
   @show x0 mu

   d = x0
   norma = norm(d)
   if norma <= tol
       println("Acabou")
   end

   s =  d
   z = H*s
   y0 = c0

   x =d
   @show x
   c1 = Vetorc(x)
   y1 = c1
   y = y1 - y0
   qsi1 = dot(s,y)
   qsi2 = dot(s,z)
   if qsi1>= 0.2*qsi2
       teta = 1
   else
       teta = (0.8*qsi2)/(qsi2-qsi1)
   end
   w = teta*y .+ (1 - teta)*z
   qsi3 = dot(s,w)
   D = (1/qsi3)*w*transpose(w)
   E = (1/qsi2)*z*transpose(z)
   H = H + D - E
   @show H


end

end
Test([1.0;1.0])
