#ExtrapolatedCritical, to evaluate properties really close to the critical point, as fast as possible

struct ExtrapolatedCritical{𝕄,𝕋} <: EoSModel
    model::𝕄
    Tc::𝕋
    Pc::𝕋
    Vc::𝕋
    dp_dT::𝕋
    d2p_dVdT::𝕋
    d3p_dV3::𝕋
    failed::Bool
end

function Base.show(io::IO,::MIME"text/plain",model::ExtrapolatedCritical)
    print(io,"ExtrapolatedCritical(Tc = $(model.Tc),Pc = $(model.Pc), Vc = $(model.Vc))")
    if model.failed
        print(io,", failed")
    end
    print(io,')')
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
    #x0 = vec2(one(Tc),log10(Vc))
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

    failed = false

    if !isfinite(d2p_dVdT) || !isfinite(d3p_dV3)
        failed = true
    end

    if d2p_dVdT <= 0 || d3p_dV3 <= 0
        failed = true
    end

    return ExtrapolatedCritical(model,Tc,Pc,Vc,dp_dT,d2p_dVdT,d3p_dV3,failed)
end

lb_volume(model::ExtrapolatedCritical,z = SA[1.0]) = zero(model.Vc)*only(z)

function ∂f∂V(model::ExtrapolatedCritical,V,T,z)
    Tc,Pc,Vc = model.Tc,model.Pc,model.Vc
    ΔT = T - Tc
    ΔV = V/only(z) - Vc
    ∂p_∂T = model.dp_dT
    ∂²p_∂T∂V = model.d2p_dVdT
    ∂³p_∂V³ = model.d3p_dV3
    Δp = (∂p_∂T)*ΔT + (∂²p_∂T∂V)*(ΔV*ΔT) + (1/6)*(∂³p_∂V³)*ΔV^3
    return -(Δp + Pc)
end

function x0_sat_pure(model::ExtrapolatedCritical,T)
    #we do this in V-T base
    Tc,Pc,Vc = model.Tc,model.Pc,model.Vc
    ΔT = T - Tc
    if model.failed
        T0 = 369.89*T/Tc
        psat0 = _propaneref_psat(T0)
        psat = Pc*psat0/4.2512e6
    else
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
    end
    Vv = volume(model.model,psat,T,phase = :v)
    ΔV_corrected = abs(Vv - Vc)
    lb_v = lb_volume(model.model,T)
    Vl = Vc - ΔV_corrected
    if Vl < 2*lb_v
        Vl = volume(model.model,psat,T,phase = :l)
    end
    return Vl,Vv
end


#in case that there isn't any antoine coefficients:
#We aproximate to RK, use the cubic antoine, and perform refinement with one Clapeyron Saturation iteration 

function x0_saturation_temperature(model::EoSModel,p,crit::Tuple)
    
    Tc,Pc,Vc = crit
    Pr = p/Pc

    if Pr > 1
        nan = zero(p)/zero(p)
        return (nan,nan,nan)
    end

    if Pr < 0.95
        A,B,C = (6.668322465137264,6.098791871032391,-0.08318016317721941)
        lnp̄ = log(p / Pc)
        T0 = Tc*(B/(A-lnp̄)-C)
        
        pii,vli,vvi = saturation_pressure(model,T0,ChemPotVSaturation(;crit))

        if isnan(pii)
            nan = zero(p)/zero(p)
            return (nan,nan,nan)
        end
        Δp = (p-pii)
        S_v = VT_entropy(model,vvi,T0)
        S_l = VT_entropy(model,vli,T0)
        ΔS = S_v - S_l
        ΔV = vvi - vli
        dpdt = ΔS/ΔV #≈ (p - pii)/(T-Tnew)
        T = T0 + Δp/dpdt
        vv = volume_virial(model,p,T)
        vl = 0.3*lb_volume(model) + 0.7*vli
        return (T,vl,vv)
    else
        model_crit = ExtrapolatedCritical(model,crit)
        return x0_saturation_temperature(model_crit,p,crit)
    end
end

function x0_saturation_temperature(model::ExtrapolatedCritical,p,crit::Tuple)
    Tc,Pc,Vc = model.Tc,model.Pc,model.Vc
    #ΔT = T - Tc
    Δp = p - Pc
    #∂p_∂T = model.dp_dT
    if model.failed
        T0 = 369.89*T/Tc
        p0 = Pc*p/4.2512e6
        T0 = _propaneref_tsat(p0)
        Tsat = 369.89*T0/Tc
    else
        ∂²p_∂T∂V = model.d2p_dVdT
        ∂³p_∂V³ = model.d3p_dV3
        fp(ΔT) = pressure(model,Vc + sqrt(-6*ΔT*∂²p_∂T∂V/∂³p_∂V³),ΔT + Tc) - p
        prob = Roots.ZeroProblem(fp,-0.01*Tc)
        Tsat  = Tc + Roots.solve(prob)
    end
    Vl = volume(model.model,p,Tsat,phase = :l)
    
    ΔV_corrected = abs(Vl - Vc)
    lb_v = lb_volume(model.model)
    #Vmax = -2*second_virial_coefficient(model.model,Tsat)
    Vv0 = Vc + 2*ΔV_corrected #overestimate vapor volume to use it as an initial guess
    Vv = volume(model.model,p,Tsat,vol0 = Vv0)
    #if Vl < 2*lb_v
    #    Vl = volume(model.model,p,Tsat,phase = :l)
    #end
    return Tsat,Vl,Vv
end

export ExtrapolatedCritical