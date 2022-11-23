struct ogUNIFACParam <: EoSParam
    A::PairParam{Float64}
    R::SingleParam{Float64}
    Q::SingleParam{Float64}
end

abstract type ogUNIFACModel <: UNIFACModel end

struct ogUNIFAC{c<:EoSModel} <: ogUNIFACModel
    components::Array{String,1}
    groups::GroupParam
    params::ogUNIFACParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
    unifac_cache::UNIFACCache
end

@registermodel ogUNIFAC
export ogUNIFAC

"""
    ogUNIFACModel <: UNIFACModel

    ogUNIFAC(components::Vector{String};
    puremodel = PR, 
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

## Input parameters
- `R`: Single Parameter (`Float64`)  - Normalized group Van der Vals volume
- `Q`: Single Parameter (`Float64`) - Normalized group Surface Area
- `A`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

UNIFAC (UNIQUAC Functional-group Activity Coefficients) activity model.

Original formulation.

The Combinatorial part corresponds to an GC-averaged modified UNIQUAC model. The residual part iterates over groups instead of components.

```
Gᴱ = nRT(gᴱ(comb) + gᴱ(res))
```

Combinatorial part:
```
gᴱ(comb) = ∑[xᵢlog(Φᵢ/xᵢ) + 5qᵢxᵢlog(θᵢ/Φᵢ)]
θᵢ = qᵢxᵢ/∑qᵢxᵢ
Φᵢ = rᵢxᵢ/∑rᵢxᵢ
rᵢ = ∑Rₖνᵢₖ for k ∈ groups
qᵢ = ∑Qₖνᵢₖ for k ∈ groups
```
Residual Part:
```
gᴱ(residual) = -v̄∑XₖQₖlog(∑ΘₘΨₘₖ)
v̄ = ∑∑xᵢνᵢₖ for k ∈ groups,  for i ∈ components
Xₖ = (∑xᵢνᵢₖ)/v̄ for i ∈ components 
Θₖ = QₖXₖ/∑QₖXₖ
Ψₖₘ = exp(-(Aₖₘ/T)
```

## References
1. Fredenslund, A., Gmehling, J., Michelsen, M. L., Rasmussen, P., & Prausnitz, J. M. (1977). Computerized design of multicomponent distillation columns using the UNIFAC group contribution method for calculation of activity coefficients. Industrial & Engineering Chemistry Process Design and Development, 16(4), 450–462. [doi:10.1021/i260064a004](https://doi.org/10.1021/i260064a004)

"""
ogUNIFAC

function ogUNIFAC(components;
    puremodel = PR,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false, kwargs...)

    groups = GroupParam(components, ["Activity/UNIFAC/ogUNIFAC/ogUNIFAC_groups.csv"];group_userlocations = group_userlocations, verbose=verbose)

    params = getparams(groups, ["Activity/UNIFAC/ogUNIFAC/ogUNIFAC_like.csv", "Activity/UNIFAC/ogUNIFAC/ogUNIFAC_unlike.csv"]; userlocations=userlocations, asymmetricparams=["A"], ignore_missing_singleparams=["A"], verbose=verbose)
    A  = params["A"]
    R  = params["R"]
    Q  = params["Q"]
    
    _puremodel = init_puremodel(puremodel,groups.components,pure_userlocations,verbose)
    packagedparams = ogUNIFACParam(A,R,Q)
    references = String[]
    cache = UNIFACCache(groups,packagedparams)
    model = ogUNIFAC(groups.components,groups,packagedparams,_puremodel,references,cache)
    return model
end

function lnγ_comb(model::ogUNIFACModel,V,T,z)
    Q = model.params.Q.values
    R = model.params.R.values

    v  = model.groups.n_flattenedgroups

    x = z ./ sum(z)

    r =[sum(v[i][k]*R[k] for k in @groups) for i in @comps]
    q =[sum(v[i][k]*Q[k] for k in @groups) for i in @comps]

    Φ = r/sum(x[i]*r[i] for i ∈ @comps)
    θ = q/sum(x[i]*q[i] for i ∈ @comps)
    lnγ_comb = @. log(Φ)+(1-Φ)-5*q*(log(Φ/θ)+(1-Φ/θ))
    return lnγ_comb
end

function Ψ(model::ogUNIFACModel,V,T,z)
    A = model.params.A.values
    return @. exp(-A/T)
end

ogUNIFACGroups = [raw"[CX4;H3]" "CH3";
raw"[CX4;H2]" "CH2";
raw"[CX4;H1]" "CH";
raw"[CX4;H0]" "C";
raw"[CX3;H2]=[CX3;H1]" "CH2=CH";
raw"[CX3;H1]=[CX3;H1]" "CH=CH";
raw"[CX3;H2]=[CX3;H0]" "CH2=C";
raw"[CX3;H1]=[CX3;H0]" "CH=C";
raw"[cX3;H1]" "ACH";
raw"[cX3;H0]" "AC";
raw"[cX3;H0][CX4;H3]" "ACCH3";
raw"[cX3;H0][CX4;H2]" "ACCH2";
raw"[cX3;H0][CX4;H1]" "ACCH";
raw"[OX2;H1]" "OH";
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
raw"[CX4;H2;R][OX2;R]" "THF";
raw"[CX3;H1;R][OX2;R]" "THF";
raw"[CX4;H3][NX3;H2]" "CH3NH2";
raw"[CX4;H2][NX3;H2]" "CH2NH2";
raw"[CX4;H1][NX3;H2]" "CHNH2";
raw"[CX4;H3][NX3;H1]" "CH3NH";
raw"[CX4;H2][NX3;H1]" "CH2NH";
raw"[CX4;H1][NX3;H1]" "CHNH";
raw"[CX4;H3][NX3;H0]" "CH3N";
raw"[CX4;H2][NX3;H0]" "CH2N";
raw"[c][NX3;H2]" "ACNH2";
raw"[cX3;H1]1:[cX3;H1]:[cX3;H1]:[nX2;H0]:[cX3;H1]:[cX3;H1]:1" "C5H5N";
raw"[cX3;H0]1:[cX3;H1]:[cX3;H1]:[nX2;H0]:[cX3;H1]:[cX3;H1]:1" "C5H4N";
raw"[cX3;H1]1:[cX3;H0]:[cX3;H1]:[nX2;H0]:[cX3;H1]:[cX3;H1]:1" "C5H4N";
raw"[cX3;H1]1:[cX3;H1]:[cX3;H0]:[nX2;H0]:[cX3;H1]:[cX3;H1]:1" "C5H4N";
raw"[cX3;H1]1:[cX3;H1]:[cX3;H1]:[nX2;H0]:[cX3;H0]:[cX3;H1]:1" "C5H4N";
raw"[cX3;H1]1:[cX3;H1]:[cX3;H1]:[nX2;H0]:[cX3;H1]:[cX3;H0]:1" "C5H4N";
raw"[cX3;H0]1:[cX3;H0]:[cX3;H1]:[nX2;H0]:[cX3;H1]:[cX3;H1]:1" "C5H3N";
raw"[cX3;H0]1:[cX3;H1]:[cX3;H0]:[nX2;H0]:[cX3;H1]:[cX3;H1]:1" "C5H3N";
raw"[cX3;H0]1:[cX3;H1]:[cX3;H1]:[nX2;H0]:[cX3;H0]:[cX3;H1]:1" "C5H3N";
raw"[cX3;H0]1:[cX3;H1]:[cX3;H1]:[nX2;H0]:[cX3;H1]:[cX3;H0]:1" "C5H3N";
raw"[cX3;H1]1:[cX3;H0]:[cX3;H0]:[nX2;H0]:[cX3;H1]:[cX3;H1]:1" "C5H3N";
raw"[cX3;H1]1:[cX3;H0]:[cX3;H1]:[nX2;H0]:[cX3;H0]:[cX3;H1]:1" "C5H3N";
raw"[cX3;H1]1:[cX3;H0]:[cX3;H1]:[nX2;H0]:[cX3;H1]:[cX3;H0]:1" "C5H3N";
raw"[cX3;H1]1:[cX3;H1]:[cX3;H0]:[nX2;H0]:[cX3;H0]:[cX3;H1]:1" "C5H3N";
raw"[cX3;H1]1:[cX3;H1]:[cX3;H0]:[nX2;H0]:[cX3;H1]:[cX3;H0]:1" "C5H3N";
raw"[cX3;H1]1:[cX3;H1]:[cX3;H1]:[nX2;H0]:[cX3;H0]:[cX3;H0]:1" "C5H3N";
raw"[CX4;H3][CX2]#[NX1]" "CH3CN";
raw"[CX4;H2][CX2]#[NX1]" "CH2CN";
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
raw"[CX3,cX3](=[OX1])[OX2,oX2]" "COO";
raw"[SiX4,SiX3,SiX5;H3]" "SIH3";
raw"[SiX4,SiX3,SiX5,SiX2;H2]" "SIH2";
raw"[SiX4,SiX3,SiX5,SiX2,SiX1;H1]" "SIH";
raw"[Si]" "SI";
raw"[SiH2][O]" "SIH2O";
raw"[SiH1][O]" "SIHO";
raw"[SiH0][O]" "SIO";
raw"[CX4;H3][NX3;H0]1[CX4;H2][CX4;H2][CX4;H2][CX3;H0]1=[OX1;H0]" "NMP";
raw"[CX4;H0]([F])([Cl])([Cl])[Cl]" "CCL3F";
raw"C(F)(Cl)Cl" "CCL2F";
raw"ClC(Cl)F" "HCCL2F";
raw"C(Cl)F" "HCCLF";
raw"Cl[CX4;H0](F)(F)" "CCLF2";
raw"Cl[CX4;H1](F)F" "HCCLF2";
raw"ClC(F)(F)F" "CCLF3";
raw"ClC(Cl)(F)F" "CCL2F2";
raw"[CX3;H0](=[OX1])[NX3;H2]" "AMH2";
raw"[CX3;H0](=[OX1])[NX3;H1][CX4;H3]" "AMHCH3";
raw"[CX3;H0](=[OX1])[NX3;H1][CX4;H2]" "AMHCH2";
raw"[CX3;H0](=[OX1])[NX3;H0]([CX4;H3])[CX4;H3]" "AM(CH3)2";
raw"[CX3;H0](=[OX1])[NX3;H0]([CX4;H3])[CX4;H2]" "AMCH3CH2";
raw"[CX3;H0](=[OX1])[NX3;H0]([CX4;H2])[CX4;H2]" "AM(CH2)2";
raw"[CX4;H2]([OX2;H1])[CX4;H2][OX2;H0]" "C2H5O2";
raw"[CX4;H1]([OX2;H1])[CX4;H2][OX2;H0]" "C2H4O2";
raw"[CX4;H2]([OX2;H1])[CX4;H1][OX2;H0]" "C2H4O2";
raw"[CX4;H3][SX2]" "CH3S";
raw"[CX4;H2][SX2]" "CH2S";
raw"[CX4,CX3,CX2;H1][S]" "CHS";
raw"[CX4;H2]1[CX4;H2][OX2;H0][CX4;H2][CX4;H2][NX3;H1]1" "MORPH";
raw"[cX3;H1]1[cX3;H1][cX3;H1][sX2;H0][cX3;H1]1" "C4H4S";
raw"[cX3;H1]1[cX3;H1][cX3;H1][sX2;H0][cX3;H0]1" "C4H3S";
raw"[cX3;H1]1[cX3;H0][cX3;H1][sX2;H0][cX3;H1]1" "C4H3S";
raw"[cX3;H0]1[cX3;H0][cX3;H1][sX2;H0][cX3;H1]1" "C4H2S";
raw"[cX3;H0]1[cX3;H1][cX3;H0][sX2;H0][cX3;H1]1" "C4H2S";
raw"[cX3;H0]1[cX3;H1][cX3;H1][sX2;H0][cX3;H0]1" "C4H2S";
raw"[cX3;H1]1[cX3;H0][cX3;H0][sX2;H0][cX3;H1]1" "C4H2S";
raw"[cX3;H1]1[cX3;H0][cX3;H1][sX2;H0][cX3;H0]1" "C4H2S";
raw"[cX3;H1]1[cX3;H1][cX3;H0][sX2;H0][cX3;H0]1" "C4H2S";
raw"[NX2H0]=[CX2H0]=[OX1H0]" "NCO";
raw"[CX4;H2][SX4](=O)(=O)[CX4;H2]" "(CH2)2SU";
raw"[CX4;H2][SX4](=O)(=O)[CX4;H1]" "CH2CHSU";
raw"[c]1:[c]:[n]:[c]:[n]:1" "IMIDAZOL";
raw"C(F)(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F" "BTI"]