#Eduardo O. Menezes
using LinearAlgebra
using Random

#Função objetivo
function fun(x::Array{Float64,1})
    f = (x[1]-3)^2+(x[2]-2)^2
end

#Função que converte bináio para real

function BinParaReal(Nbits::Int64,VetorBits::Array{Int64,1},xmin::Float64,xmax::Float64)
        N = 0
        for i = Nbits:1
                N = N + VetorBits(i) * 2^(Nbits - i)
        end
        delta = (xmax - xmin) / ((2^N) - 1)

        return delta
end

function IndividuoParaReal(BinIndividuo::Array{Int64,1},Xmin::Array{Float64,1},xmax::Array{Float64,1},n::Int64,Nbits::Int64)
        xs = zeros(n,1)
        inicio = 1
        fim = Nbits
        for i = 1:n
                bin = zeros(n,1)
                for j = inicio : fim
                        bin = append!(bin,BinIndividuo(j))
                end
                inicio = fim + 1
                fim = fim + Nbits
                valor = BinParaReal(Nbits,bin,xmin(i),xmax(i))
                xs = append!(xs,valor)
        end
        return xs
end

function Populacao(Np,Nbits,n)
        pop = zeros(2*Np,n*Nbits)
        for i = 1:Np
                individuo = zeros(1,n*Nbits)
                for j = 1: n*Nbits
                        individuo = append!(individuo, rand(1))
                end
                pop = append!(pop,individuo)
        end
        return pop
end

function ZerosParaUm(k)
        val = 0
        if k == 0
                val = 1
        end
        return val
end

function Mutacao(pop, Np,n,Nbits,TaxMut)
        nbitsMod = round(n*Nbits*TaxMut*Np)
        for i = 1:nbitsMod
                PosIndiv = rand(1:Np)
                Posbit = rand(1:n*Nbits)
                val = pop[PosIndiv,Posbit]
                pop[PosIndiv,Posbit] = ZerosParaUm(val)
        end
        return pop
end

function Fitness(pop,Np,n,Nbits,xmin,xmax)
        v = zeros(Np,1)
        for i = 1:Np
                xi = IndividuoParaReal(pop[i],xmin,xmax,n,Nbits)
                v = append!(v,fun(xi))
        end
        return v
end

function Selecao(P,v,Np)
        Npa = Np / 2
        pai = zeros(Npa,1)
        mae = zeros(Npa,1)
        for i = 1:Npa
                p1 = rand(1:Np)
                p2 = rand(1:Np)
                if v[p1] < v[p2]
                        pai = append!(pai,p1)
                else
                        pai = append!(pai,p2)
                end
                m1 = rand(1:Np)
                m2 = rand(1:Np)
                if v[m1] < v[m2]
                        mae = append!(mae,m1)
                else
                        mae = append!(mae,m2)
                end
        end
        par = [pai,mae]
        return par
end

function Cruzamento(P,pais,maes,Np,n,Nbits)
        P2 = zeros(Np,1)
        for i = 1 : Np/2
                p = P[pais[i]]
                m = P[maes[i]]
                crossover = rand(1:n*Nbits)
                if crossover == n*Nbits
                        crossover = crossover - 1
                end
                f1 = zeros(Np,1)
                f2 = zeros(Np,1)
                for j = 1:crossover
                        f1 = append!(f1,p[j])
                        f2 = append!(f2,m[j])
                end
                for j = crossover+1 : n*Nbits
                        f1 = append!(f1,m[j])
                        f2 = append!(f2,p[j])
                end
        P2 = append!(P2,f1)
        P2 = append!(P2,f2)
        end
        return P2
end

function Informacoes(geracao,v,P,xmin,xmax,Np,n,Nbits)
        melhor = 100E100
        pior = -100E100
        posicao = 0
        for i =1:Np
                if v[i] < melhor
                        melhor = v[i]
                        posicao = i
                end
                if v[i] > pior
                        pior = v[i]
                end
        end
        xm = IndividuoParaReal(P[posicao],xmin,xmax,n,Nbits)
        media = (melhor+pior) / 2
        println("Geração: ", geracao, "objetivo: ", melhor, "media: ", media, "individuo", xm)
end

function Geneticos(n,Nbits,Np,xmin,xmax,Ng)
        P = Populacao(Np,n,Nbits)
        for gera = 1:Ng
                v = Fitness(P,Np,n,Nbits,xmin,xmax)
                display(Informacoes(gera,v,P,xmin,xmax,Np,n,Nbits))
                pares = Selecao(P,v,Np)
                pais = pares[1]
                maes = pares[2]
                P2 = Cruzamento(P,pais,maes,Np,n,Nbits)
                PM = Mutacao(P2,Np,n,Nbits,1/100)
                P = PM
        end
        return P
end

Geneticos(2,10,130,[0,0],[10,10],100)
