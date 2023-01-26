

abstract type AbstractFixPoint end

"""
  Solvers.fixpoint(f,x0::Real,method=SSFixPoint())
does a fixpoint iteration convergence on a series of real numbers.
f is a function that should 
the following strategies:
 - `SSFixPoint(dampingfactor = 1.0)`: performs succesive substitutions until convergence is met. 
α = `dampingfactor` is determines a buffer for each iteration, defined as `x1 = α*f(x1) + (1-α)*x1`

- `Aitken()`: uses Aitken's delta-squared process to accelerate convergence of the series. recommended for harmonic iterates.

"""
function fixpoint end

struct SSFixPoint{T<:Real} <: AbstractFixPoint 
    dampingfactor::T
end

SSFixPoint(;dampingfactor=1.0) = SSFixPoint(dampingfactor)

function promote_method(method::SSFixPoint,T)
    return SSFixPoint(T(method.dampingfactor))
end

struct AitkenFixPoint <: AbstractFixPoint end

function promote_method(method::AitkenFixPoint,T)
    return method
end

function convergence(xold,xi,atol,rtol)
    not_finite = false
    for xii in xi
        if !isfinite(xii)
            not_finite = true
            break
        end
    end
    not_finite && return (true,false) #terminate, with nan
    xi == xold && return (true,true) #terminate, with current number
    if xi isa Number
        Δx = abs(xi-xold)
    else
        Δx = norm((xi[i] - xold[i] for i in eachindex(xold,xi)))
    end
    normxi = norm(xi)
    if abs(Δx) <= max(atol,normxi*rtol)
        return (true,true) #terminate, with current number
    end
    return (false,false) #keep iterating
end

function fixpoint(f,x0,
    method::AbstractFixPoint = SSFixPoint();
    atol=zero(eltype(x0)),
    rtol=8eps(one(eltype(x0))), 
    max_iters=100,
    return_last = false,
    inplace = false)
    _,atol,rtol = promote(one(eltype(x0)),atol,rtol)
    method = promote_method(method,eltype(x0))
    return _fixpoint(f,x0,method,atol,rtol,max_iters,return_last,inplace)
end

function _fixpoint(f::F,
    x0::T,
    method::SSFixPoint,
    atol::T = zero(T),
    rtol::T =8*eps(T),
    max_iters=100,
    return_last = false,
    inplace = false) where {F,T<:Real}
    
    nan = (0*atol)/(0*atol)
    xi = f(x0)
    converged,finite = convergence(x0,xi,atol,rtol)
    converged && return ifelse(finite,xi,nan)
    itercount = 1
    xold = x0
    α = method.dampingfactor
    while itercount < max_iters
        xi = α*f(xi) + (1-α)*xi  
        converged,finite = convergence(xold,xi,atol,rtol)
        converged && return ifelse(finite,xi,nan)    
        itercount +=1
        xold = xi
    end
    return ifelse(return_last,xi,nan)
end

function _fixpoint(f::F,
    x00::T,
    method::AitkenFixPoint,
    atol::T = zero(T),
    rtol::T =8*eps(T),
    max_iters=100,
    return_last = false,
    inplace = false) where {F,T<:Real}

    nan = (0*atol)/(0*atol)
    itercount = 2
    x3 = x00
    x1 = f(x3)
    converged,finite = convergence(x3,x1,atol,rtol)
    converged && return ifelse(finite,x1,nan)  
    x2 = f(x1)
    converged,finite = convergence(x1,x2,atol,rtol)
    converged && return ifelse(finite,x2,nan)

    while itercount < max_iters
        itercount += 1
        λ2 = (x2 - x1)/(x1-x3)
        dx = -(λ2/(1-λ2))*(x2-x1)
        x3 = x2 + dx
        converged,finite = convergence(x2,x3,atol,rtol)
        converged && return ifelse(finite,x3,nan)
        
        itercount += 1
        x1 = f(x3)
        converged,finite = convergence(x3,x1,atol,rtol)
        converged && return ifelse(finite,x1,nan)
        
        itercount += 1
        x2 = f(x1)
        converged,finite = convergence(x1,x2,atol,rtol)
        converged && return ifelse(finite,x2,nan)
    end
    return ifelse(return_last,x2,nan)
end

function swap!(x,y)
    for i in eachindex(x)
        xi,yi = x[i],y[i]
        y[i] = xi
        x[i] = yi
    end
    return x,y
end

function _fixpoint(f!::F,
    x0::X where {X <:AbstractVector{T}},
    method::SSFixPoint,
    atol::T = zero(T),
    rtol::T =8*eps(T),
    max_iters=100,
    return_last = false,
    inplace = false) where {F,T<:Real}
    
    nan = (0*atol)/(0*atol)
    if inplace isa Bool #hacks!. we save an argument here.
        xi = copy(x0)
    else
        xi = inplace
        inplace = true
    end
    xi = f!(xi,x0)
    converged,finite = convergence(x0,xi,atol,rtol)
    if converged
        if finite
            return xi
        else
            xi .= nan
            return xi 
        end
    end
    itercount = 1
    α = method.dampingfactor
    
    if inplace
        xi,xold = swap!(x0,xi) 
        #our initial x0 now its our latest point
        #whereas the storage created for xi is used for xold
    else
        xold = copy(x0)
    end

    while itercount < max_iters
        xi = f!(xi,xold)
        xi .*= α
        xi .+= (1 .- α) .* xold
        converged,finite = convergence(xold,xi,atol,rtol)
        if converged
            if finite
                return xi
            else
                xi .= nan
                return xi 
            end
        end
        itercount +=1
        xold .= xi
    end
    !return_last && (xi .= nan)
    return xi
    
end
