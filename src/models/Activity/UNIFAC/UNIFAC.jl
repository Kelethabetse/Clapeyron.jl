struct UNIFACParam <: EoSParam
    A::PairParam{Float64}
    B::PairParam{Float64}
    C::PairParam{Float64}
    R::SingleParam{Float64}
    Q::SingleParam{Float64}
end

abstract type UNIFACModel <: ActivityModel end

struct UNIFAC{c<:EoSModel} <: UNIFACModel
    components::Array{String,1}
    groups::GroupParam
    params::UNIFACParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
    unifac_cache::UNIFACCache
end

@registermodel UNIFAC
const modUNIFAC = UNIFAC
export UNIFAC

"""
    UNIFACModel <: ActivityModel
    UNIFAC(components::Vector{String};
    puremodel = PR,
    userlocations = String[], 
    pure_userlocations = String[],
    verbose = false)
## Input parameters
- `R`: Single Parameter (`Float64`)  - Normalized group Van der Vals volume
- `Q`: Single Parameter (`Float64`) - Normalized group Surface Area
- `A`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter
- `B`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter
- `C`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter
## Input models
- `puremodel`: model to calculate pure pressure-dependent properties
## Description
UNIFAC (UNIQUAC Functional-group Activity Coefficients) activity model.
Modified UNIFAC (Dortmund) implementation.
The Combinatorial part corresponds to an GC-averaged modified [`UNIQUAC`](@ref) model. The residual part iterates over groups instead of components.
```
Gᴱ = nRT(gᴱ(comb) + gᴱ(res))
```
Combinatorial part:
```
gᴱ(comb) = ∑[xᵢlog(Φ'ᵢ) + 5qᵢxᵢlog(θᵢ/Φᵢ)]
θᵢ = qᵢxᵢ/∑qᵢxᵢ
Φᵢ = rᵢxᵢ/∑rᵢxᵢ
Φ'ᵢ = rᵢ^(0.75)/∑xᵢrᵢ^(0.75)
rᵢ = ∑Rₖνᵢₖ for k ∈ groups
qᵢ = ∑Qₖνᵢₖ for k ∈ groups
```
Residual Part:
```
gᴱ(residual) = -v̄∑XₖQₖlog(∑ΘₘΨₘₖ)
v̄ = ∑∑xᵢνᵢₖ for k ∈ groups,  for i ∈ components
Xₖ = (∑xᵢνᵢₖ)/v̄ for i ∈ components 
Θₖ = QₖXₖ/∑QₖXₖ
Ψₖₘ = exp(-(Aₖₘ + BₖₘT + CₖₘT²)/T)
```
## References
1. Fredenslund, A., Gmehling, J., Michelsen, M. L., Rasmussen, P., & Prausnitz, J. M. (1977). Computerized design of multicomponent distillation columns using the UNIFAC group contribution method for calculation of activity coefficients. Industrial & Engineering Chemistry Process Design and Development, 16(4), 450–462. [doi:10.1021/i260064a004](https://doi.org/10.1021/i260064a004)
2. Weidlich, U.; Gmehling, J. A modified UNIFAC model. 1. Prediction of VLE, hE, and.gamma..infin. Ind. Eng. Chem. Res. 1987, 26, 1372–1381.
"""
UNIFAC

function UNIFAC(components::Vector{String};
    puremodel = PR,
    userlocations = String[],
    group_userlocations = String[], 
    pure_userlocations = String[],
    verbose = false)
    
    groups = GroupParam(components, ["Activity/UNIFAC/UNIFAC_groups.csv"]; group_userlocations = group_userlocations, verbose = verbose)

    params = getparams(groups, ["Activity/UNIFAC/UNIFAC_like.csv", "Activity/UNIFAC/UNIFAC_unlike.csv"]; userlocations=userlocations, asymmetricparams=["A","B","C"], ignore_missing_singleparams=["A","B","C"], verbose=verbose)
    A  = params["A"]
    B  = params["B"]
    C  = params["C"]
    R  = params["R"]
    Q  = params["Q"]
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = UNIFACParam(A,B,C,R,Q)
    references = String["10.1021/i260064a004"]
    cache = UNIFACCache(groups,packagedparams)
    model = UNIFAC(components,groups,packagedparams,_puremodel,references,cache)
    return model
end

function recombine_impl!(model::UNIFACModel)
    recombine_unifac_cache!(model.unifac_cache,model.groups,model.params)
    recombine!(model.puremodel)
end

function activity_coefficient(model::UNIFACModel,V,T,z)
    return exp.(@f(lnγ_comb)+ @f(lnγ_res))
end

function lnγ_comb(model::UNIFACModel,V,T,z)
    x = z ./ sum(z)
    r =model.unifac_cache.r
    q =model.unifac_cache.q
    q_p = model.unifac_cache.q_p
    Φ = r/dot(x,r)
    Φ_p = q_p/dot(x,q_p)
    θ = q/dot(x,q)
    lnγ_comb = @. log(Φ_p)+(1-Φ_p)-5*q*(log(Φ/θ)+(1-Φ/θ))
    return lnγ_comb
end

function lnγ_SG(model::UNIFACModel,V,T,z)

    x = z ./ sum(z)

    r =model.unifac_cache.r
    q =model.unifac_cache.q

    Φ = r/dot(x,r)
    θ = q/dot(x,q)
    lnγ_SG = @. -5*q*(log(Φ/θ)+(1-Φ/θ))
    return lnγ_SG
end

function lnγ_res(model::UNIFACModel,V,T,z)
    v  = model.groups.n_flattenedgroups
    _ψ = @f(Ψ)
    lnΓ_ = @f(lnΓ,_ψ)
    lnΓi_ = @f(lnΓi,_ψ)
    lnγ_res_ =  [sum(v[i][k].*(lnΓ_[k].-lnΓi_[i][k]) for k ∈ @groups) for i ∈ @comps]
    return lnγ_res_
end

function lnΓ(model::UNIFACModel,V,T,z,ψ = @f(ψ))
    Q = model.params.Q.values
    v  = model.groups.n_flattenedgroups
    x = z ./ sum(z)
    X = sum(v[i][:]*x[i] for i ∈ @comps) ./ sum(sum(v[i][k]*x[i] for k ∈ @groups) for i ∈ @comps)
    θ = X.*Q / dot(X,Q)
    lnΓ_ = Q.*(1 .-log.(sum(θ[m]*ψ[m,:] for m ∈ @groups)) .- sum(θ[m]*ψ[:,m]./sum(θ[n]*ψ[n,m] for n ∈ @groups) for m ∈ @groups))
    return lnΓ_
end

function lnΓi(model::UNIFACModel,V,T,z,ψ = @f(ψ))
    Q = model.params.Q.values
    v  = model.groups.n_flattenedgroups
    ψ = @f(Ψ)
    X = [v[i][:] ./ sum(v[i][k] for k ∈ @groups) for i ∈ @comps]
    θ = [X[i][:].*Q ./ sum(X[i][n]*Q[n] for n ∈ @groups) for i ∈ @comps]
    lnΓi_ = [Q.*(1 .-log.(sum(θ[i][m]*ψ[m,:] for m ∈ @groups)) .- sum(θ[i][m]*ψ[:,m]./sum(θ[i][n]*ψ[n,m] for n ∈ @groups) for m ∈ @groups)) for i ∈ @comps]
    return lnΓi_
end

function Ψ(model::UNIFACModel,V,T,z)
    A = model.params.A.values
    B = model.params.B.values
    C = model.params.C.values
    return @. exp(-(A+B*T+C*T^2)/T)
end

function excess_g_SG(model::UNIFACModel,p,T,z)
    lnγ = lnγ_SG(model,p,T,z)
    return sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
end

function excess_g_res(model::UNIFACModel,p,T,z)
    lnγ = lnγ_res(model,p,T,z)
    return sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
end

UNIFACGroups = [raw"[CX4;H3;!R]" "CH3";
raw"[CX4;H2;!R]" "CH2";
raw"[CX4;H1;!R]" "CH";
raw"[CX4;H0;!R]" "C";
raw"[CX3;H2]=[CX3;H1]" "CH2=CH";
raw"[CX3;H1]=[CX3;H1]" "CH=CH";
raw"[CX3;H2]=[CX3;H0]" "CH2=C";
raw"[CX3;H1]=[CX3;H0]" "CH=C";
raw"[cX3;H1]" "ACH";
raw"[cX3;H0]" "AC";
raw"[cX3;H0][CX4;H3]" "ACCH3";
raw"[cX3;H0][CX4;H2]" "ACCH2";
raw"[cX3;H0][CX4;H1]" "ACCH";
raw"[OH1;$([OH1][CX4H2])]" "OH(P)";
raw"[CX4;H3][OX2;H1]" "CH3OH";
raw"[OH2]" "H2O";
raw"[cX3;H0;R][OX2;H1]" "ACOH";
raw"[CX4;H3][CX3](=O)" "CH3CO";
raw"[CX4;H2][CX3](=O)" "CH2CO";
raw"[CX3H1](=O)" "CHO";
raw"[CH3][CX3;H0](=[O])[O]" "CH3COO";
raw"[CX4;H2][CX3](=[OX1])[OX2]" "CH2COO";
raw"[CX3;H1](=[OX1])[OX2]" "HCOO";
raw"[CH3;!R][OH0;!R]" "CH3O";
raw"[CH2;!R][OH0;!R]" "CH2O";
raw"[C;H1;!R][OH0;!R]" "CHO";
raw"[CX4;H2;R][OX2;R][CX4;H2;R]" "THF";
raw"[CX4;H3][NX3;H2]" "CH3NH2";
raw"[CX4;H2][NX3;H2]" "CH2NH2";
raw"[CX4;H1][NX3;H2]" "CHNH2";
raw"[CX4;H3][NX3;H1]" "CH3NH";
raw"[CX4;H2][NX3;H1]" "CH2NH";
raw"[CX4;H1][NX3;H1]" "CHNH";
raw"[CX4;H3][NX3;H0]" "CH3N";
raw"[CX4;H2][NX3;H0]" "CH2N";
raw"[c][NX3;H2]" "ACNH2";
raw"[cX3H1][n][cX3H1]" "AC2H2N";
raw"[cX3H0][n][cX3H1]" "AC2HN";
raw"[cX3H0][n][cX3H0]" "AC2N";
raw"[CX4;H3][CX2]#[NX1]" "CH3CN";
raw"[CX4;H2][CX2]#[NX1]" "CH2CN";
raw"[CX3,cX3](=[OX1])[OX2H0,oX2H0]" "COO";
raw"[CX3](=[OX1])[O;H1]" "COOH";
raw"[CX3;H1](=[OX1])[OX2;H1]" "HCOOH";
raw"[CX4;H2;!$(C(Cl)(Cl))](Cl)" "CH2CL";
raw"[CX4;H1;!$(C(Cl)(Cl))](Cl)" "CHCL";
raw"[CX4;H0](Cl)" "CCL";
raw"[CX4;H2;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)" "CH2CL2";
raw"[CX4;H1;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)" "CHCL2";
raw"[CX4;H0;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)" "CCL2";
raw"[CX4;H1;!$([CX4;H0](Cl)(Cl)(Cl)(Cl))](Cl)(Cl)(Cl)" "CHCL3";
raw"[CX4;H0;!$([CX4;H0](Cl)(Cl)(Cl)(Cl))](Cl)(Cl)(Cl)" "CCL3";
raw"[CX4;H0]([Cl])([Cl])([Cl])([Cl])" "CCL4";
raw"[c][Cl]" "ACCL";
raw"[CX4;H3][NX3](=[OX1])([OX1])" "CH3NO2";
raw"[CX4;H2][NX3](=[OX1])([OX1])" "CH2NO2";
raw"[CX4;H1][NX3](=[OX1])([OX1])" "CHNO2";
raw"[cX3][NX3](=[OX1])([OX1])" "ACNO2";
raw"C(=S)=S" "CS2";
raw"[SX2H][CX4;H3]" "CH3SH";
raw"[SX2H][CX4;H2]" "CH2SH";
raw"c1cc(oc1)C=O" "FURFURAL";
raw"[OX2;H1][CX4;H2][CX4;H2][OX2;H1]" "DOH";
raw"[I]" "I";
raw"[Br]" "BR";
raw"[CX2;H1]#[CX2;H0]" "CH=-C";
raw"[CX2;H0]#[CX2;H0]" "C=-C";
raw"[SX3H0](=[OX1])([CX4;H3])[CX4;H3]" "DMSO";
raw"[CX3;H2]=[CX3;H1][CX2;H0]#[NX1;H0]" "ACRY";
raw"[$([Cl;H0]([C]=[C]))]" "CL-(C=C)";
raw"[CX3;H0]=[CX3;H0]" "C=C";
raw"[cX3][F]" "ACF";
raw"[CX4;H3][N]([CX4;H3])[CX3;H1]=[O]" "DMF";
raw"[NX3]([CX4;H2])([CX4;H2])[CX3;H1](=[OX1])" "HCON(CH2)2";
raw"C(F)(F)F" "CF3";
raw"C(F)F" "CF2";
raw"C(F)" "CF";
raw"[CH2;R]" "CY-CH2";
raw"[CH1;R]" "CY-CH";
raw"[CH0;R]" "CY-C";
raw"[OH1;$([OH1][CX4H1])]" "OH(S)";
raw"[OH1;$([OH1][CX4H0])]" "OH(T)";
raw"[CX4H2;R][OX2;R]" "CY-CH2O";
raw"[CX4H2;R][OX2;R]" "TRIOXAN";
raw"[CX4H0][NH2]" "CNH2";
raw"[OX1H0]=[C;R][NX3H0;R][CH3]" "NMP";
raw"[OX1H0]=[CH0X3;R][H0;R][CH2]" "NEP";
raw"[OX1H0;!R]=[CX3H0;R][NX3H0;R][C;!R]" "NIPP";
raw"[OX1H0;!R]=[CH0X3;R][NX3H0;R][CH0;!R]" "NTBP";
raw"[CX3H0](=[OX1H0])[NX3H2]" "CONH2";
raw"[OX1H0;!R]=[CX3H0;!R][NH1X3;!R][CH3;!R]" "CONHCH3";
raw"[CH2X4;!R][NH1X3;!R][CX3H0;!R]=[OX1H0;!R]" "CONHCH2"]