using NewtonMethod
using Test, LinearAlgebra

@testset "NewtonMethod.jl" begin
    x_0 = BigFloat(rand())

    # NewtonMethod 1-1, nonconvergence  
    f(x) = x^2 + 2
    f′(x)= 2x - 2
    xsol = newtonroot(f, f′; x_0 = x_0)
    @test isa(xsol[1],Nothing) 

    # NewtonMethod 1-2-1, convergence---------------------
    f(x) = x^2 - 2x + 1
    f′(x) = 2x - 2
    xsol = newtonroot(f, f′; x_0 = x_0)
    @test norm(xsol[1] - 1.0) < 1.0E-7

    # NewtonMethod 1-2-2, convergence
    f(x) = x^2 - 2x + 1
    f′(x) = 2x - 2
    xsol = newtonroot(f, f′; x_0 = nothing)
    @test norm(xsol[1] - 1.0) < 1.0E-7

    # NewtonMethod 1-2-3, convergence
    f(x) = x^2 - 2x + 1
    f′(x) = 2x - 2
    xsol = newtonroot(f, f′; x_0 = "test")
    @test norm(xsol[1] - 1.0) < 1.0E-7
    #-----------------------------------------------------

    # NewtonMethod 1-3, reduce maxiter
    f(x) = x^2 - 2x + 1
    f′(x) = 2x - 2
    xsol = newtonroot(f, f′; x_0 = 100.0, maxiter = 5)
    @test isa(xsol[1], Nothing)

    # NewtonMethod 1-4, loosen tolerance
    f(x) = x^2 - 2x + 1
    xsol = newtonroot(f, f′; x_0 = 100.0, tol = 1E-3)
    @test norm(xsol[1] - 1.0) >1E-4
    
    #--------------------------------------------------------
    # NewtonMethod 2-1, nonconvergence
    f(x) = x^2 + 2
    xsol = newtonroot(f; x_0 = x_0)
    @test isa(xsol[1],Nothing)  

    # NewtonMethod 2-2-1, convergence-----------
    f(x) = x^2 - 2x + 1
    xsol = newtonroot(f; x_0 = x_0)
    @test norm(xsol[1] - 1.0) < 1.0E-7

    # NewtonMethod 2-2-2, convergence
    f(x) = x^2 - 2x + 1
    xsol = newtonroot(f; x_0 = nothing)
    @test norm(xsol[1] - 1.0) < 1.0E-7

    # NewtonMethod 2-2-3, convergence
    f(x) = x^2 - 2x + 1
    xsol = newtonroot(f; x_0 = "test")
    @test norm(xsol[1] - 1.0) < 1.0E-7
    #----------------------------------------

    # NewtonMethod 2-3, reduce maxiter
    f(x) = x^2 - 2x + 1
    xsol = newtonroot(f; x_0 = 100.0, maxiter = 5)
    @test isa(xsol[1], Nothing)

    # NewtonMethod 2-4, loosen tolerance
    f(x) = x^2 - 2x + 1
    xsol = newtonroot(f; x_0 = 100.0, tol = 1E-3)
    @test norm(xsol[1] - 1.0) >1E-4
end
