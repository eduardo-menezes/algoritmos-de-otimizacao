clearconsole()
#
# Steepest descent
#

include("golden_LAumentado.jl")



function Objetivo(x::Array{Float64,1},L::Array{Float64,1})

    @assert length(x)==length(L) "Usuario maneta!!!"

    soma = 0.0
    for ele in length(x)
        soma = soma + x[ele]*L[ele]
    end

    return soma

end

function dObjetivo(x::Array{Float64,1},L::Array{Float64,1})

    #=
    # Numero de elementos
    ne = length(x)

    # Inicializa o vetor de derivadas
    d = zeros(ne)

    for m =1:ne

       soma = 0.0
       for ele=1:ne
    	   soma = soma + δ(m,ele)*L[ele]
       end
    end #m
    =#

	return L

end




function Steepest(fun::Function,grad::Function,x::Array{Float64,1},)

    # A dimensão do problema deve vir de x0 e de grad.
    n = size(x,1)
    tol=1E-5
    # Loop do programa
    #while true
    cont = 0
    for i=1:6
         # Gradiente
         df = grad(x)
         @show df

         # Norma do gradiente
         norma = norm(df)
         @show norma

         # Verifica pela convergência
         if norma <= tol
             println("Solução obtida::",x)
             break
         end

         # Direção de minimização normalizada
         d = (-1.0 .* df) ./ norma
         @show d

         # Intervalo inicial para o LS
         a,b =  Intervalo_Inicial(fun,x,d)
         @show a,b

         # LS
         alpha = LS_GS(fun,x,d,a,b,tol)
         @show alpha

         # Novo ponto
         x = x .+ alpha*d
         @show x
         cont = cont +1
         @show cont
     end

    return x

 end


#
# Teste
#





#
# Solução de um problema utilizando LA
#
function Exemplo_LA()


    function objetivo(x)
        a = (x[1]-10)^3 + (x[2]-20)^3
        return a
    end

    function dobjetivo(x)
        b = [3*(x[1]-10)^2 ; 3*(x[2]-20)^2]
        return b
    end

    function restricoes(x)

        g1 = 100.0 - (x[1] - 5.0)^2 - (x[2] - 5.0)^2
        g2 = - 82.81 + (x[1]-6.0)^2 + (x[2]-5.0)^2
        g3 = 13.0 - x[1]
        g4 = x[1] - 100.0
        g5 = - x[2]
        g6 = x[2] - 100.0

        gg = [g1 ; g2 ; g3; g4 ; g5 ; g6]

        return gg

    end


    # Cada derivada em uma coluna
    function drestricoes(x)

     dg = [ -2*(x[1]-5) 2*(x[1]-6) -1.0 1.0 0.0 0.0;
            -2*(x[2]-5) 2*(x[2]-5) 0.0 0.0 -1.0 1.0 ]

            return dg

    end

    #
    # Operador <>
    #
    function Oper(a)
        op = max(0.0,a)

        return op
    end

    #
    # Algoritomo do LA
    #

    # Nova função Objetivo
    function LA(x,c,μ)

        # Calcula o objetivo
        f = objetivo(x)

        # Calcula as restrições
        g = restricoes(x)

        # Monta o LA
        L = f + (c/2)*sum((Oper.(μ/c+g)).^2)

        #=
        L = f + (c/2)*( Oper(mu[1]/c + g[1])^2 +
                        Oper(mu[2]/c + g[2])^2 +
                        Oper(mu[3]/c + g[3])^2  )

        L = f
        for j=1:3
            L =  L + (c/2)*Oper(mu[j]/c + g[j])^2
        end
        =#

        return L

    end

    #
    # Novo vetor gradiente
    #
    function dLA(x,c,μ)

         # Calculamos as restrições neste ponto
         g = restricoes(x)

         # Calcula o gradiente do objetivo
         df = dobjetivo(x)

         # Calcula os gradientes das restrições
         dg = drestricoes(x)

         # Monta o vetor gradiente "composto"
         dL = copy(df)

         # Monta por Loop
         for j=1:size(dg,2)

             dL .= dL + Oper(μ[j] + c*g[j]).*vec(dg[:,j])

         end

         return dL

    end


    ###################################################################
    ### JÁ TEMOS TODAS AS ROTINAS PARA SOLUCIONAR O PROBLEMA POR LA
    ###################################################################

    # Soluciona por LA

    x = [15.0 ; 5.0]

    c = 10.0

    μ = [0.0 ; 0.0 ; 0.0; 0.0 ; 0.0 ; 0.0] #-> primeira iter é penalização pura
    g = [0.0 ; 0.0 ; 0.0; 0.0 ; 0.0 ; 0.0]

    # Criar um alias para as funções
    ALA(x) = LA(x,c,μ)

    DLA(x) = dLA(x,c,μ)


    # Loop externo do LA
    τ = 1.1
    nexternas = 100

    for iter=1:nexternas

        # Chama o Steepest
        sol = Steepest(ALA,DLA,x)

        # Atualiza o ponto inicial
        x = copy(sol)
        @show x

        # Calcula as restrições neste novo ponto
        g = restricoes(x)
        @show g

        # Atualiza a penalização e os multiplicadores
        c = τ*c
        @show c
        μ .= Oper.( μ + c*g )
        @show μ
        @show nexternas

    end
    @show x μ g
    # Retorna as variáveis do problema
    return x, μ, g

end # Exemplo_LA

Exemplo_LA()
