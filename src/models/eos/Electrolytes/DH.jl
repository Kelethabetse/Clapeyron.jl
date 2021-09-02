struct DHParam <: EoSParam
    sigma::SingleParam{Float64}
    charge::SingleParam{Float64}
end

abstract type IonModel <: EoSModel end
abstract type DHModel <: IonModel end

struct DH{ϵ<:RSPModel} <: DHModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::DHParam
    RSPmodel::ϵ
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel DH
export DH
function DH(components; RSPmodel=ConstW, SAFTlocations=String[], userlocations=String[], ideal_userlocations=String[], verbose=false)
    params,sites = getparams(components, append!(["Electrolytes/charges.csv","properties/molarmass.csv"],SAFTlocations); userlocations=userlocations, verbose=verbose)
    icomponents = 1:length(components)
    println(params)
    params["sigma"].values .*= 1E-10
    sigma = params["sigma"]
    charge = params["charge"]

    packagedparams = DHParam(sigma,charge)

    references = [""]

    init_RSPmodel = RSPmodel(components)

    model = DH(components,icomponents, packagedparams, init_RSPmodel, 1e-12,references)
    return model
end

function a_ion(model::DHModel, V, T, z)
    σ = model.params.sigma.values
    Z = model.params.charge.values
    ϵ_r = RSP(model.RSPmodel,V,T,z)

    x = z ./ sum(z)
    ρ = N_A*sum(z)/V

    s = e_c^2/(4π*ϵ_0*ϵ_r*k_B*T)
    κ = (4π*s*ρ*sum(x[i]*Z[i]^2 for i ∈ @comps))^(1/2)
    y = σ*κ
    χ = @. 3/y^3*(3/2+log(1+y)-2*(1+y)+1/2*(1+y)^2)
    return -1/3*s*κ*sum(x[i]*Z[i]^2*χ[i] for i ∈ @comps)
end