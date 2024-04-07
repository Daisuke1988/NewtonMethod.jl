module NewtonMethod
    # julia REPL で add Revise PkgTemplates を忘れないこと！

    # Newton 1
    using NLsolve, LinearAlgebra

    function newtonroot(f::Function, f′::Function; x_0::Any, tol = 1.0E-7, maxiter=1000)
        # error handling of x_0
        if isa(x_0, Real)
            x = x_0
        elseif isa(x_0, Nothing)
            x = rand()
        else
            println("x_0 is NOT real number. I replace it with random real number.")
            x = rand()
        end

        normdiff = Inf
        iter = 0

        #f′ = Polynomials.derivative(f)
        while normdiff > tol && iter < maxiter
            x1 = x - f(x) / f′(x)
            normdiff = norm(x1 - x)
            x = x1
            iter += 1
        end
        
        # error handling of nonconvergence
        if iter == maxiter || isnan(x)
            x = nothing
            return (; x, iter)
            println("Number of iterations reaches maxiter!! NOT CONVERGE!!")
        else
            return (; x, iter)  # note: multiple returns are not defined as tuple but array 
        end

    end

    export newtonroot

    # for checking 
    @show f(x) = 3x^3 + 2x^2 - 4x + 3
    #@show f′(x) = 9x^2 + 4x -4
    #@show xsol = newtonroot(f, f′; x_0 = 0.0)
 



    # Newton 2===============================================================
    using ForwardDiff
    function newtonroot(f::Function; x_0::Any, tol = 1.0E-7, maxiter=1000)
        # error-handling of x_0
        if isa(x_0, Real)
            x = x_0
        elseif is(x_0, Nothing)
            x = rand()
        else
            println("x_0 is NOT real number. I replace it with random real number.")
            x = rand()
        end
        
        normdiff = Inf
        iter = 0
        D(f) = x -> ForwardDiff.derivative(f, x) # define taking-derivative operator
        f′ = D(f)

        while normdiff > tol && iter < maxiter
            x1 = x - f(x) / f′(x)  # 
            normdiff = norm(x1 - x)     # do NOT forget updating normdiff
            x = x1
            iter += 1
        end

        if iter == maxiter || isnan(x)
            x = nothing
            return (; x, iter)
            println("Number of iterations reaches maxiter!! NOT CONVERGE!!")
        end

        return (; x, iter)  # for making tuple as struct-like from (; ...) 
    end

    export newtonroot

    # for checking
    #f(x) = x^7 - 6x -1 #   f(x) = x^6 -6x + 1
    #@show xsol = newtonroot(f; x_0 = 100)



end