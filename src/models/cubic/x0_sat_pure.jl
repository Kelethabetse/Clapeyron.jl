#this plugs the Superancillary equations as the initial point for the supported cubics

function x0_sat_pure(model::Union{vdWModel,PRModel,RKModel},T)
    Tc = model.params.Tc.values[1]
    if Tc < T
        nan = zero(T)/zero(T)
        return (nan,nan)
    end
    a,b,_ = cubic_ab(model,1e-3,T)
    T̃ = T*R̄*b/a
    Vv = chebyshev_vapour_volume(model,T̃,b)
    Vl = chebyshev_liquid_volume(model,T̃,b)
    return vl,vv
end