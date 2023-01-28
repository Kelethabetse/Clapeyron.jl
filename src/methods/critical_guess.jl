#ExtrapolatedCritical, to evaluate properties really close to the critical point, as fast as possible

struct ExtrapolatedCritical{𝕄,𝕋} <: EoSModel
    model::𝕄
    Tc::𝕋
    Pc::𝕋
    Vc::𝕋
    dp_dT::𝕋
    d2p_dVdT::𝕋
    d3p_dV3::𝕋
end

function Base.show(io::IO,::MIME"text/plain",model::ExtrapolatedCritical)
    print(io,"ExtrapolatedCritical(Tc = $(model.Tc),Pc = $(model.Pc), Vc = $(model.Vc))")
end

function Base.show(io::IO,model::ExtrapolatedCritical)
    print(io,"ExtrapolatedCritical(Tc = $(model.Tc),Pc = $(model.Pc), Vc = $(model.Vc))")
end

function crit_pure(model::ExtrapolatedCritical)
    return model.Tc,model.Pc,model.Vc
end

"""
    ExtrapolatedCritical(model,crit = crit_pure(model))

returns an extrapolation of a pure model's V-T surface, based on a taylor expansion at the critical point:
```
Δp = (∂p/∂T)*ΔT + (∂²p/∂T∂V)*(ΔV*ΔT) + (1/6)*(∂³p/∂V³)*ΔV^3
ΔT = T - Tc
Δp = p - pc
ΔV = V - Vc
```
"""
function ExtrapolatedCritical(model::EoSModel,crit = crit_pure(model))
    Tc,Pc,Vc = crit
    #we call this because it should be already compiled
    x0 = vec2(one(Tc),log10(Vc))
    obj = ObjCritPure(model,Tc,x0)
    function f(x)
        V,T = x
        ∂²A∂V², ∂³A∂V³ = ∂²³f(model, V, T, SA[1.0])
        return SA[∂²A∂V²,∂³A∂V³]
    end

    J = ForwardDiff.jacobian(f,vec2(Vc,Tc))
    #f(x) = ∂²A∂V², ∂³A∂V³ = ∂p∂V, ∂²p∂V²
    #J(x) = [∂2p∂V∂T, ∂2p∂V2
    #          ∂3p∂V2∂T, ∂3p∂V3]

    dp_dT = Solvers.derivative(_T -> pressure(model,Vc,_T),Tc)
    d2p_dVdT,d3p_dV3 = diagvalues(J)
    return ExtrapolatedCritical(model,Tc,Pc,Vc,dp_dT,d2p_dVdT,d3p_dV3)
end

lb_volume(model::ExtrapolatedCritical,z = SA[1.0]) = zero(model.Vc)*only(z)

function pressure(model::ExtrapolatedCritical,V,T,z = SA[1.0])
    Tc,Pc,Vc = model.Tc,model.Pc,model.Vc
    ΔT = T - Tc
    ΔV = V/only(z) - Vc
    ∂p_∂T = model.dp_dT
    ∂²p_∂T∂V = model.d2p_dVdT
    ∂³p_∂V³ = model.d3p_dV3
    Δp = (∂p_∂T)*ΔT + (∂²p_∂T∂V)*(ΔV*ΔT) + (1/6)*(∂³p_∂V³)*ΔV^3
    return Δp + Pc
end

function x0_sat_pure(model::ExtrapolatedCritical,T)
    #we do this in V-T base
    Tc,Pc,Vc = model.Tc,model.Pc,model.Vc
    ΔT = T - Tc
    #Δp = p - Pc
    #∂p_∂T = model.dp_dT
    ∂²p_∂T∂V = model.d2p_dVdT
    ∂³p_∂V³ = model.d3p_dV3

    ΔV = sqrt(-6*ΔT*∂²p_∂T∂V/∂³p_∂V³)
    Vx = ΔV + Vc
    ρc = 1/Vc
    Δρ = 1/Vx - ρc
    Bρ = Δρ/sqrt(-ΔT/Tc)
    ρl,ρv = ρc + abs(Bρ)*sqrt(-ΔT/Tc), ρc - abs(Bρ)*sqrt(-ΔT/Tc)
    Vl0,Vv0 = 1/ρl, 1/ρv
    psat = pressure(model,Vl0,T)
    
    Vv = volume(model.model,psat,T,phase = :v)
    ΔV_corrected = abs(Vv - Vc)
    lb_v = lb_volume(model.model,T)
    Vl = Vc - ΔV_corrected
    if Vl < 2*lb_v
        Vl = volume(model.model,psat,T,phase = :l)
    end
    return Vl,Vv
end

function x0_saturation_temperature(model::ExtrapolatedCritical,p,crit::Tuple)
    Tc,Pc,Vc = model.Tc,model.Pc,model.Vc
    #ΔT = T - Tc
    Δp = p - Pc
    #∂p_∂T = model.dp_dT
    ∂²p_∂T∂V = model.d2p_dVdT
    ∂³p_∂V³ = model.d3p_dV3
    fp(ΔT) = pressure(model,Vc + sqrt(-6*ΔT*∂²p_∂T∂V/∂³p_∂V³),ΔT + Tc) - p
    prob = Roots.ZeroProblem(fp,-0.01*Tc)
    Tsat  = Tc + Roots.solve(prob)
    Vv = volume(model.model,p,Tsat,phase = :v)
    ΔV_corrected = abs(Vv - Vc)
    lb_v = lb_volume(model.model)
    Vl = Vc - ΔV_corrected
    if Vl < 2*lb_v
        Vl = volume(model.model,p,Tsat,phase = :l)
    end
    return Tsat,Vl,Vv
end

export ExtrapolatedCritical