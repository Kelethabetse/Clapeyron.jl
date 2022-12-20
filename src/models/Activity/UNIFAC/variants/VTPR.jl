abstract type VTPRUNIFACModel <: UNIFACModel end


struct VTPRUNIFACCache <: EoSModel
    components::Vector{String}
    m::Vector{Float64}
end

struct VTPRUNIFAC{c<:EoSModel} <: VTPRUNIFACModel
    components::Array{String,1}
    groups::GroupParam
    params::UNIFACParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
    unifac_cache::VTPRUNIFACCache
end

@registermodel VTPRUNIFAC
export VTPRUNIFAC

"""
    VTPRUNIFACModel <: UNIFACModel

    VTPRUNIFAC(components::Vector{String};
    puremodel = BasicIdeal,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

## Input parameters
- `Q`: Single Parameter (`Float64`) - Normalized group Surface Area
- `A`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter
- `B`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter
- `C`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

## Description
UNIFAC (UNIQUAC Functional-group Activity Coefficients) activity model.

Modified UNIFAC (Dortmund) implementation, only residual part, activity model used for the Volume-Translated Peng-Robinson (VTPR) EoS.

The residual part iterates over groups instead of components.

```
Gᴱ = nRT(gᴱ(res))
gᴱ(res) = -v̄∑XₖQₖlog(∑ΘₘΨₘₖ)
v̄ = ∑∑xᵢνᵢₖ for k ∈ groups,  for i ∈ components
Xₖ = (∑xᵢνᵢₖ)/v̄ for i ∈ components
Θₖ = QₖXₖ/∑QₖXₖ
Ψₖₘ = exp(-(Aₖₘ + BₖₘT + CₖₘT²)/T)
```

## References
1. Fredenslund, A., Gmehling, J., Michelsen, M. L., Rasmussen, P., & Prausnitz, J. M. (1977). Computerized design of multicomponent distillation columns using the UNIFAC group contribution method for calculation of activity coefficients. Industrial & Engineering Chemistry Process Design and Development, 16(4), 450–462. "doi:[10.1021/i260064a004](https://doi.org/10.1021/i260064a004)"
2. Weidlich, U.; Gmehling, J. A modified UNIFAC model. 1. Prediction of VLE, hE, and.gamma..infin. Ind. Eng. Chem. Res. 1987, 26, 1372–1381.
3. Ahlers, J., & Gmehling, J. (2001). Development of an universal group contribution equation of state. Fluid Phase Equilibria, 191(1–2), 177–188. [doi:10.1016/s0378-3812(01)00626-4](https://doi.org/10.1016/s0378-3812(01)00626-4)
"""
VTPRUNIFAC

function VTPRUNIFAC(components;
    puremodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false, kwargs...)

    groups = GroupParam(components, ["Activity/UNIFAC/VTPR/VTPR_groups.csv"];group_userlocations = group_userlocations, verbose=verbose)

    params = getparams(groups, ["Activity/UNIFAC/VTPR/VTPR_like.csv", "Activity/UNIFAC/VTPR/VTPR_unlike.csv"]; userlocations=userlocations,  asymmetricparams=["A","B","C"], ignore_missing_singleparams=["A","B","C"], verbose=verbose)
    A  = params["A"]
    B  = params["B"]
    C  = params["C"]
    Q  = params["Q"]
    R = deepcopy(Q)
    R.values .= 0
    cache = VTPRUNIFACCache(groups)
    _puremodel = init_puremodel(puremodel,groups.components,pure_userlocations,verbose)
    packagedparams = UNIFACParam(A,B,C,R,Q)
    references = String["10.1016/S0378-3812(01)00626-4"]
    model = VTPRUNIFAC(groups.components,groups,packagedparams,_puremodel,references,cache)
    return model
end

function recombine_unifac_cache!(cache::VTPRUNIFACCache,groups,params)
    group_sum!(cache.m,groups,nothing)
    return cache
end

function VTPRUNIFACCache(groups::GroupParam)
    m = group_sum(groups,nothing)
    return VTPRUNIFACCache(groups.components,m)
end

function excess_gibbs_free_energy(model::VTPRUNIFACModel,V,T,z)
    lnγ = lnγ_res(model,V,T,z)
    return sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
end

VTPRGroups = [raw"[CX4;H3]" "CH3";
raw"[CX4;H2]" "CH2";
raw"[CX4;H1]" "CH";
raw"[CX4;H0]" "C";
raw"[CX3;H2]=[CX3;H1]" "CH2=CH";
raw"[CX3;H1]=[CX3;H1]" "CH=CH";
raw"[CX3;H2]=[CX3;H0]" "CH2=C";
raw"[CX3;H1]=[CX3;H0]" "CH=C";
raw"[CX3;H0]=[CX3;H0]" "C=C";
raw"[cX3;H1]" "ACH";
raw"[cX3;H0]" "AC";
raw"[cX3;H0][CX4;H3]" "ACCH3";
raw"[cX3;H0][CX4;H2]" "ACCH2";
raw"[cX3;H0][CX4;H1]" "ACCH";
raw"[OH1;$([OH1][CX4H2])]" "OH(P)";
raw"[OH1;$([OH1][CX4H1])]" "OH(S)";
raw"[OH1;$([OH1][CX4H0])]" "OH(T)";
raw"[CX4;H3][OX2;H1]" "CH3OH";
raw"[OH2]" "H2O";
raw"[cX3;H0;R][OX2;H1]" "ACOH";
raw"[CX4;H3][CX3](=O)" "CH3CO";
raw"[CX4;H2][CX3](=O)" "CH2CO";
raw"[CX3H1](=O)" "CHO";
raw"[CH3][CX3;H0](=[O])[O]" "CH3COO";
raw"[CX4;H2][CX3](=[OX1])[OX2]" "CH2COO";
raw"[CX3;H1](=[OX1])[OX2]" "HCOO";
raw"[CH3][O]" "CH3O";
raw"[CH2][O]" "CH2O";
raw"[C;H1][O]" "CHO";
raw"[CX4;H3][NX3;H2]" "CH3NH2";
raw"[CX4;H2][NX3;H2]" "CH2NH2";
raw"[CX4;H1][NX3;H2]" "CHNH2";
raw"[CX4H0][NH2]" "CNH2";
raw"[CX4;H3][NX3;H1]" "CH3NH";
raw"[CX4;H2][NX3;H1]" "CH2NH";
raw"[CX4;H1][NX3;H1]" "CHNH";
raw"[CX4;H3][NX3;H0]" "CH3N";
raw"[CX4;H2][NX3;H0]" "CH2N";
raw"[c][NX3;H2]" "ACNH2";
raw"[CX4;H3][CX2]#[NX1]" "CH3CN";
raw"[CX4;H2][CX2]#[NX1]" "CH2CN";
raw"[CX4;H2;!$(C(Cl)(Cl))](Cl)" "CH2CL";
raw"[CX4;H1;!$(C(Cl)(Cl))](Cl)" "CHCL";
raw"[CX4;H0](Cl)" "CCL";
raw"[CX4;H2;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)" "CH2CL2";
raw"[CX4;H1;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)" "CHCL2";
raw"[CX4;H0;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)" "CCL2";
raw"[CX4;H0;!$([CX4;H0](Cl)(Cl)(Cl)(Cl))](Cl)(Cl)(Cl)" "CCL3";
raw"[CX4;H0]([Cl])([Cl])([Cl])([Cl])" "CCL4";
raw"[c][Cl]" "ACCL";
raw"[SX2H][CX4;H3]" "CH3SH";
raw"[SX2H][CX4;H2]" "CH2SH";
raw"C(=S)=S" "CS2";
raw"c1cc(oc1)C=O" "FURFURAL";
raw"[OX2;H1][CX4;H2][CX4;H2][OX2;H1]" "DOH";
raw"[I]" "I";
raw"[Br]" "BR";
raw"[SX3H0](=[OX1])([CX4;H3])[CX4;H3]" "DMSO";
raw"[CX4;H3][N]([CX4;H3])[CX3;H1]=[O]" "DMF";
raw"[NX3]([CX4;H2])([CX4;H2])[CX3;H1](=[OX1])" "HCON(..";
raw"[CH2;R]" "CY-CH2";
raw"[CH1;R]" "CY-CH";
raw"[CH0;R]" "CY-C";
raw"[CX4;H2;R][OX2;R]" "THF";
raw"[CX3;H1;R][OX2;R]" "THF";
raw"[CX4H2;R][OX2;R]" "CY-CH2O";
raw"[CX4H2;R][OX2;R]" "TRIOXAN";
raw"[CX4;H1;!$([CX4;H0](Cl)(Cl)(Cl)(Cl))](Cl)(Cl)(Cl)" "CHCL3";
raw"[OX1H0]=[C;R][NX3H0;R][CH3]" "NMP";
raw"[OX1H0]=[CH0X3;R][H0;R][CH2]" "NEP";
raw"[OX1H0;!R]=[CX3H0;R][NX3H0;R][C;!R]" "NIPP";
raw"[OX1H0;!R]=[CH0X3;R][NX3H0;R][CH0;!R]" "NTBP";
raw"[CX4;H3][NX3](=[OX1])([OX1])" "CH3NO2";
raw"[CX4;H2][NX3](=[OX1])([OX1])" "CH2NO2";
raw"[CX4;H1][NX3](=[OX1])([OX1])" "CHNO2"]