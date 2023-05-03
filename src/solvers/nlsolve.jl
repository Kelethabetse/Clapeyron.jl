#nlsolve functionality

function _identity(v::Vector{T},mstyle::Union{NLSolvers.OutOfPlace,NLSolvers.InPlace}) where T
    out = Matrix{T}(undef, length(v), length(v))
    out .= false
    for i in 1:length(v)
        out[i,i] = true
    end
    return out
end

function _identity(x,mstyle)
     return I + x .* x' .* false
end

function _identity(x::StaticArrays.MArray{S, T, N, L},mstyle::Union{NLSolvers.OutOfPlace,NLSolvers.InPlace}) where {S, T, N, L}
    res = StaticArraysCore.MMatrix{L, L, T, L*L}(undef)
    res .= false
    for i in LinearAlgebra.diagind(res)
        res[i] = true
    end
    return res
end


#this is to allow our problems to use the optimization interface instead of the NEqProblem interfas
struct MeritObjective{TP,T1,T2,T3,T4,T5}
    prob::TP
    F::T1
    FJ::T2
    Fx::T3
    Jx::T4
    d::T5
end

function _opt_options(neq::NEqOptions)
    return OptimizationOptions(;
    x_abstol = 2*neq.f_abstol,
    x_reltol = 2*neq.f_reltol,
    x_norm = x -> norm(x, Inf),
    g_abstol = 2*neq.f_abstol,
    g_reltol = 2*neq.f_reltol,
    g_norm = x -> norm(x, Inf),
    f_limit = neq.f_limit,
    f_abstol = 0.0, #TODO, maybe neq.f_abstol
    f_reltol = 0.0,
    nm_tol = 1e-8,
    maxiter = neq.maxiter,
    show_trace = false,
)
end

function NLSolvers.value(mo::MeritObjective, x)
    Fx = mo.F(mo.Fx, x)
    dot(Fx,Fx)/2
end

#not correct but it is only used for Newton
function NLSolvers.upto_hessian(mo::MeritObjective, ∇f, ∇²f, x)
    Fx_ne, Jx_ne = mo.FJ(∇f, ∇²f, x)
    f = dot(Fx_ne,Fx_ne) / 2
    #mul!(∇f,transpose(Jx_ne),Fx_ne)
    return f, ∇f, ∇²f
end

function NLSolvers.upto_gradient(mo::MeritObjective, ∇f, x)  #Fx is the gradient and Jx is the Hessian
    Fx_ne, Jx_ne = mo.FJ(mo.Fx, mo.Jx, x)
    f = dot(Fx_ne,Fx_ne)/2
    # this is the gradient
    mul!(∇f,transpose(Jx_ne),Fx_ne)
    # As you may notice, this can be expensive... Because the gradient
    # is going to be very simple. May want to create a
    # special type or way to hook into trust regions here. We can exploit
    # that we only need the cauchy and the newton steps, not any shifted
    # systems. There is no need to get the full hessian. because these two
    # steps are don't need these multiplies
    # This is the Hessian
    #Jx .= Jx_sq' * Jx_sq
    return f, ∇f
end


function solve(problem::NEqProblem,
    x0,
    method::LineSearch = LineSearch(Newton()),
    options = NEqOptions())
    
    mstyle = NLSolvers.mstyle(problem)
    cache = NLSolvers.preallocate_qn_caches(mstyle, x0)

    Fx = cache.d #we use the cache here.
    F = problem.R.F
    FJ = problem.R.FJ
    
    B0 = _identity(x0,mstyle)
    #@show B0
    merit = MeritObjective(problem, F, FJ, Fx, B0, nothing)
    s0 = (x0,B0) 
    optprob = OptimizationProblem(merit, problem.bounds, problem.manifold, nothing, mstyle, nothing)
    opt_options = _opt_options(options)
    return NLSolvers._solve(mstyle,optprob,s0,method,opt_options,cache)
end

solve(problem,x,method,options) = NLSolvers.solve(problem,x,method,options)

"""
    function nlsolve(f!,x0,method=TrustRegion(Newton(), NWI()), options=NEqOptions(),chunk = ForwardDiff.Chunk{2}())


Given a function `f!(result,x)` that returns a system of equations,
`nlsolve(f!,x0)` returns a `NLSolvers.ConvergenceInfo` struct that contains the results of the non-linear solving procedure.

Uses `NLSolvers.jl` as backend, the jacobian is calculated with `ForwardDiff.jl`, with the specified `chunk` size

To obtain the underlying solution vector, use [`x_sol`](@ref)

To see available solvers and options, check `NLSolvers.jl`
"""
function nlsolve(f!,x0,method=TrustRegion(Newton(), NWI()),options=NEqOptions(),chunk = ForwardDiff.Chunk{2}())
    vector_objective = autoVectorObjective(f!,x0,chunk)
    nl_problem = NEqProblem(vector_objective)
    return nlsolve(nl_problem, x0,method, options)
end

function nlsolve(nl_problem::NEqProblem,x0,method =TrustRegion(Newton(), NWI()),options=NEqOptions())
    #hook our MeritObjective here.
    return NLSolvers.solve(nl_problem, x0, method, options)
end

function autoVectorObjective(f!,x0,chunk)
    Fcache = x0 .* false
    Tag = typeof(ForwardDiff.Tag(f!, typeof(Fcache)))
    jconfig = ForwardDiff.JacobianConfig(f!,x0,x0,chunk)
    function j!(J,x)
        ForwardDiff.jacobian!(J,f!,Fcache,x,jconfig)
        J
    end
    function fj!(F,J,x)
        ForwardDiff.jacobian!(J,f!,F,x,jconfig)
        F,J
    end
    #=
    jv_dual = ForwardDiff.Dual{Tag, eltype(Fcache), 1}
    jv_cache = similar(Fcache,jv_dual)
    #y = ForwardDiff.Dual{Tag, eltype(Fcache), 1}
    function jv!(x)
        #set cache
        for i in eachindex(x)
            xi = x[i]
            jv_cache[i] = jv_dual(xi,ForwardDiff.Partials((xi,)))
        end
        f!(jv_cache,x)
        for i in eachindex(x)
            x[i] = ForwardDiff.partials(jv_cache[i],1)
        end
        return x
    end =#

    return NLSolvers.VectorObjective(f!,j!,fj!,nothing)
end

#=
function auto_jacvec(f, x, v)
    T = typeof(ForwardDiff.Tag(f, typeof(x)))
    vv = reshape(v, axes(x))
    y = ForwardDiff.Dual{T, eltype(x), 1}.(x, ForwardDiff.Partials.(tuple.(vv)))
    vec(partials.(vec(f(y)), 1))
end =#


#= only_fj!: NLsolve.jl legacy form:

function only_fj!(F, J, x)
    # shared calculations begin
    # ...
    # shared calculation end
    if !(F == nothing)
        # mutating calculations specific to f! goes here
    end
    if !(J == nothing)
        # mutating calculations specific to j! goes
    end
end
=#
function only_fj!(fj!::T) where T
    function _f!(F,x)
        fj!(F,nothing,x)
        F
    end

    function _fj!(F,J,x)
        fj!(F,J,x)
        F,J
    end

    function _j!(J,x)
        fj!(nothing,J,x)
        J
    end
    _jv!(x) = nothing
    return NLSolvers.VectorObjective(_f!,_j!,_fj!,_jv!) |> NEqProblem
    # return NLSolvers.VectorObjective(f!,j!,fj!,jv!) |> NEqProblem
end


