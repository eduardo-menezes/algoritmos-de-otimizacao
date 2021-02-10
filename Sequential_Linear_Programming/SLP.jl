clearconsole()
#
# Pacotes necessários:
# LinearAlgebra
# Combinatorics
import Pkg

Pkg.add("CSV")
Pkg.update("Atom")
Pkg.add("Combinatorics")
Pkg.add("DataFrames")
using LinearAlgebra
using CSV
using DataFrames
using Combinatorics
#
# Soluciona Min c'x
#      T.q  Ax <= b
# com c nx1
#     A mxn
#     b mx1
#
# retorna as variáveis primais e as duais separadas e nesta ordem e um flag
# de sucesso/insucesso. flag = 1 (OK) / -1 (devolve nos limites móveis) / -2
# problema mal posto.
#
# Esta rotina tem função puramente educacional, pois não tem nenhuma
# preocupação com eficiência (muito pelo contrário!!!)
#
# uso
#
# x,mu = LP_SLP(c,A,b,ci,cs)
#
# onde ci e cs são as restrições laterais (obrigatórias aqui!!)
# Eduardo Lenz - setembro de 2018.
# Atualização: 5/11/2018 -> situação limite onde não existe solução viável,
#                           mas para utilização com um SLP a solução é
#                           devolvida como um ponto nos limites móveis. Para isto,
#                           introduzi as entradas ci e cs.
#
#
function LP_SLP(c,Ao,bo,ci,cs,x0)

     # Vamos fazer alguns testes básicos
     ma1,na1 = size(Ao)
     mb1 = size(bo,1)
     nc1 = size(c,1)
     nci1 = size(ci,1)
     ncs1 = size(cs,1)

     @assert ma1==mb1 "Número de linhas em A e b deve ser igual ao número de restrições (sem os limites móveis)"
     @assert nc1==na1 "Número de colunas em A deve ser igual ao número de linhas em c (número de váriáveis)"
     @assert nci1==nc1 "Número de linhas em ci deve ser igual ao de váriáveis"
     @assert ncs1==nc1 "Número de linhas em cs deve ser igual ao de váriáveis"

     # Já podemos garantir o número de variáveis.
     n = nc1
     m = ma1+2*n

     # Agora vou aumentar a matriz A e o vetor b, para conter os limites móveis
     A = zeros(m,n)
     b = zeros(m)

     # Vamos lá
     A[1:ma1,:] .= Ao
     b[1:mb1] .= bo
     cont = ma1
     for i=1:n
         cont +=1
         A[cont,i]=1.0
         b[cont]= cs[i]
         cont +=1
         A[cont,i]=-1.0
         b[cont] = -ci[i]

     end
     @show A b
     # Se tivermos m>n, temos que impor multiplicadores adicionais
     sobra = m-n
     @assert sobra>=0 "Não é permitido termos mais variáveis do que restrições (ao menos as laterais né!)"

     # A dimensão do sistema aumentado será
     dimen = m+n+sobra

     # Monta o sistema aumentado, mas ainda sem impor as restrições sobre os
     # multiplicadores de Lagrange
     M = zeros(dimen,dimen)
     F = zeros(dimen)

     # Vamos copiar os blocos
     M[1:n,n+1:n+m] .= transpose(A)
     M[n+1:n+m,1:n] .= A
     F[1:n] .= -c
     F[n+1:n+m] .= b
    

     # Agora vem a parte da solução. O caso mais simples é quando m==n => sobra=0
     # neste caso, devemos ter uma solução diretamente.
     if sobra==0.0
         println("Solução direta (m==n) ")
         try
             fact = lu(M)
         catch
              println("LP::Problema mal posto")
              return -2, zeros(n),zeros(m)
         end
         sol =  vec(fact\F)

         return sol[1:n]


    else

     # Caso sobra seja maior do que zero temos que impor as condições de
     # restrições inativas e testar pelo sinal dos multiplicadores (novos e
     # originais)

     # O número de testes é dado pelo número de combinações com dimensão sobra
     # que podemos fazer nas nossas restrições. Isto dá
     #
     # nc = m! / ( sobra! n!)
     #
     #
     nc = convert(Int64,factorial(m) / (factorial(sobra)*factorial(n)))
     #println("LP::número de combinações a testar é de ",nc," zerando ", sobra, " multiplicadores por vez")

     # O julia tem o pacote Combinatorials que também permite calcular este valor diretamente
     nc2 = binomial(m,sobra)

     #println("LP::número de combinações a testar é de ",nc2," zerando ", sobra, " multiplicadores por vez")

     # A lista de todas as combinações é dada por
     combinacoes = collect(combinations(1:m,sobra))

     #println("LP::as combinações são ", combinacoes)

     # Loop por todas as combinações. A alteração para usar esta rotina
     # com um SLP e no caso de não termos uma solução para um dado conjunto
     # de limites móveis é que se não tivermos sucesso após testar todas as
     # combinações, devolvemos a melhor solução dentro dos limites móveis,
     # mesmo que isto viole alguma restrição.
     for pos_comb=1:length(combinacoes) #comb in combinacoes

          # Recupera a combinação
          comb = combinacoes[pos_comb]

          # Mostra a combinação atual
          #println("Testando a combinação ", comb)

          # Vamos montar uma matriz local
          Ml = copy(M)

          # Vamos colocar os uns nos lugares adequados
          cont = 1
          for pos in comb
                       # Parte de baixo
              Ml[n+m+cont,n+pos] = 1.0

              # Parte de cima
              Ml[n+pos,m+n+cont] = 1.0

              cont += 1


         end # pos

         # Soluciona o problema atual
         sol_ok = true
         sol = zeros(n+m+sobra)
         try

            fact = lu(Ml)
            sol = vec(fact\F)

         catch
            sol_ok = false
        end

        if sol_ok

            # Procura por multiplicadores negativos.
            flag = true
            for i in sol[n+1:end]
                if i<0
                    flag=false
                    break
                end
            end

            if flag
                dual = sol[n+1:n+m]


                println("Solução encontrada ", sol[1:n],sol[n+1:n+m])

                # Retornamos a solução primal e a dual separadas
                return sol[1:n]


            end


        end #sol_ok

     end # comb

 end
     # Se chegamos aqui, então não tivemos sucesso em nenhuma combinação.
     # Assim, devemos usar o gradiente (c) e os limites móveis
     println("Não obtivemos uma solução viável. Devolvendo uma solução na
              fronteira de minimzação (limites móveis)")
    xs = zeros(n)
    for i=1:n
        if c[i]<0.0
            xs[i]=ci[i]
        elseif c[i]>0.0
            xs[i]=cs[i]
        else
            xs[i]=x0[i]
        end
   end
   println(xs)
   return xs



end


function Principal_SLP(x)
tol = 1E-5
count = 0
daux = zeros(0)
xaux = zeros(0)
delta = 1
deltad = 0
    while delta>=0.001

        ci = [(delta-x[1]);(delta-x[2])]
        cs = [(x[1]+delta);(x[2]+delta)]
        @show ci cs
        g1 = 100-(x[1]-5)^2 - (x[2]-5)^2
        g2 = -82.81 + (x[1]-6)^2 + (x[2]-5)^2
        g3= -x[1] +13.0
        g4=  x[1] - 100.0
        g5= -x[2]
        g6 = x[2] - 100.0
        g11 = -2*(x[1]-5)
        g12 = -2*(x[2]-5)
        g21 = 2*(x[1]-6)
        g22 = 2*(x[2]-5)
        g31 = -1.0
        g32 = 0.0
        g41 = 1.0
        g42 = 0.0
        g51 = 0.0
        g52 = -1.0
        g61 = 0.0
        g62 = 1.0

        A = [g11 g12; g21 g22;g31 g32; g41 g42;g51 g52; g61 g62]

        g1 = 100.0-(x[1]-5)^2 - (x[2]-5)^2
        g2 = -82.81 + (x[1]-6)^2 + (x[2]-5)^2
        g3= -x[1] +13
        g4=  x[1] - 100
        g5= -x[2]
        g6 = x[2] - 100
        g11 = -2*(x[1]-5)
        g12 = -2*(x[2]-5)
        g21 = 2*(x[1]-6)
        g22 = 2*(x[2]-5)
        g31 = -1.0
        g32 = 0.0
        g41 = 1.0
        g42 = 0.0
        g51 = 0.0
        g52 = -1.0
        g61 = 0.0
        g62 = 1.0

        b=[-g1+ g11*x[1]+g12*x[2];-g2+ g21*x[1]+g22*x[2];-g3+ g31*x[1]+g32*x[2]; -g4+ g41*x[1]+g42*x[2]; -g5+g51*x[1]+g52*x[2]; -g6+g61*x[1]+g62*x[2]]
        c = [3*(x[1]-10)^2;3*(x[2]-20)^2]
        x = LP_SLP(c,A,b,ci,cs,x)
        xaux = append!(xaux,x)
        daux = append!(daux,delta)

        count = count +1
        @show xaux
        @show count
       if count>=3
           x1 = xaux[2*(count-1)+1]-xaux[2*(count-1)-1]
           y1 = xaux[2*(count-1)+2]-xaux[2*(count-1)]
           @show x1 y1

           x2 = xaux[2*(count-1)-1] - xaux[2*(count-1)-3]
           y2 = xaux[2*(count-1)] - xaux[2*(count-1)-2]
           @show x2 y2

           deltax = x1 - x2
           deltay = y1 - y2
           @show deltax deltay

           if (x1 < tol) & (x2 < tol)

                delta = 0.1*delta
            end


        end

        @show delta
        @show x
    end

end

Principal_SLP([15.0;5.0])
