clearconsole()

using LinearAlgebra
using Random

## Definução da função objetivo

function fun(x)
    f = (x[1]-3)^2 + (x[2]-2)^2
end

## Rotina para PSO

function PSO(Np::Int64,n::Int64,Niter::Int64,w::Float64,c1::Float64,c2::Float64)

    X = zeros(Np,n) ## Inicia o exame com zeros
     ## Inicia a velocidade do enxame com zeros
    ## Inicia o valor de F com zeros
    OMelhorGlob = 100E100




    for i = 1:Np ## Np é o número de partículas
        for j = 1:n ## n é o número de variáveis de projeto
            X[i,j] = 5*rand()
        end
    end
    veloc = zeros(Np,n)
    for i = 1:Np ## Np é o número de partículas
        for j = 1:n ## n é o número de variáveis de projeto
            veloc[i,j] =  -5 + 5*rand()
        end
    end

    @show X
    @show veloc
    F = zeros(Np)
    XMelhorPers = zeros(Np,n)
    OMelhorPers = zeros(Np)
    XMelhorGlob = zeros(n)
    for i = 1:Np
        F[i] =fun(vec(X[i,:]))
        OMelhorPers = F[i]
        XMelhorPers[i,:] .=vec(X[i,:])

        if F[i] < OMelhorGlob
            OMelhorGlob = F[i]
            XMelhorGlob .= vec(X[i,:])
        end
    end
    @show F
    for k = 1:Niter
        for i = 1:Np
            for j = 1:n
                veloc[i,j] = w*veloc[i,j]
                             + c1*rand()*(XMelhorPers[i,j] - X[i,j])
                             + c2*rand()*(XMelhorGlob[j] - X[i,j])
                veloc[i,j] = max(-1.0,min(1.0,veloc[i,j]))

            end
        end
        @show veloc
        for i = 1:Np
            for j = 1:n
                X[i,j] = X[i,j] + veloc[i,j]
                X[i,j] = max(-10.0,min(10.0,X[i,j]))
            end
            F[i] = fun(vec(X[i,:]))

            if F[i] < OMelhorPers[i]
                OMelhorPers[i] = F[i]
                XMelhorPers[i,:] .= vec(X[i,:])
            end

            if F[i] < OMelhorGlob
            OMelhorGlob = F[i]
            XMelhorGlob .= vec(X[i,:])

            end
            @show XMelhorGlob
        end
    end
    
    return OMelhorGlob
end

PSO(50,2,100,0.5,1.5,1.5)
