abstract type GCsPCSAFTModel <: sPCSAFTModel end

struct GCsPCSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    msigma3::SingleParam{Float64}
    mepsilon::SingleParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

struct GCsPCSAFT{I,PC} <: GCsPCSAFTModel
    components::Vector{String}
    groups::GroupParam
    sites::SiteParam
    params::GCsPCSAFTParam
    idealmodel::I
    pcsaftmodel::PC
    assoc_options::AssocOptions
    references::Array{String,1}
end

export GCsPCSAFT

"""
    GCsPCSAFT <: PCSAFTModel
    GCsPCSAFT(components; 
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions())
## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `m`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `k`: Pair Parameter (`Float64`) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`
## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume
## Input models
- `idealmodel`: Ideal Model
## Description
Simplified Perturbed-Chain SAFT (sPC-SAFT)
## References
1. von Solms, N., Michelsen, M. L., & Kontogeorgis, G. M. (2003). Computational and physical performance of a modified PC-SAFT equation of state for highly asymmetric and associating mixtures. Industrial & Engineering Chemistry Research, 42(5), 1098–1105. [doi:10.1021/ie020753p](https://doi.org/10.1021/ie020753p)
"""
GCsPCSAFT

function GCsPCSAFT(components;
    idealmodel=BasicIdeal,
    userlocations=String[],
    group_userlocations = String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions())
    
    groups = GroupParam(components,["SAFT/PCSAFT/GCsPCSAFT/GCsPCSAFT_groups.csv"]; group_userlocations = group_userlocations,verbose=verbose)
    gc_params,sites = getparams(groups, ["SAFT/PCSAFT/GCsPCSAFT/","properties/molarmass_groups.csv"]; userlocations=userlocations, verbose=verbose)
    components = groups.components
    params = getparams(components, ["SAFT/PCSAFT/sPCSAFT/sPCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    
    gc_mw = gc_params["Mw"]
    mw = group_sum(groups,gc_mw)

    gc_segment = gc_params["segment"]
    segment = group_sum(groups,gc_segment)

    gc_msigma3 = gc_params["msigma3"]
    gc_msigma3.values .*= 1E-30
    sigma = group_sum(groups,gc_msigma3)
    sigma.values ./= segment.values
    sigma.values .= cbrt.(sigma.values)
    sigma = sigma_LorentzBerthelot(sigma)

    gc_mepsilon = gc_params["mepsilon"]
    epsilon = group_sum(groups,gc_mepsilon)
    epsilon.values ./= segment.values

    k = params["k"]
    epsilon = epsilon_LorentzBerthelot(epsilon,k)

    gc_epsilon_assoc = gc_params["epsilon_assoc"]
    gc_bondvol = gc_params["bondvol"]
    gc_bondvol,gc_epsilon_assoc = assoc_mix(gc_bondvol,gc_epsilon_assoc,[],assoc_options) #combining rules for association

    comp_sites,idx_dict = gc_to_comp_sites(sites,groups)
    assoc_idx = gc_to_comp_assoc_idx(gc_bondvol,comp_sites,idx_dict)
    assoc_idxs,outer,inner,outer_size,inner_size = assoc_idx.values,assoc_idx.outer_indices,assoc_idx.inner_indices,assoc_idx.outer_size,assoc_idx.inner_size
    _comp_bondvol = Float64[gc_bondvol.values.values[i] for i ∈ assoc_idxs]
    _comp_epsilon_assoc = Float64[gc_epsilon_assoc.values.values[i] for i ∈ assoc_idxs]
    compval_bondvol = Compressed4DMatrix(_comp_bondvol,outer,inner,outer_size,inner_size)
    compval_epsilon_assoc = Compressed4DMatrix(_comp_epsilon_assoc,outer,inner,outer_size,inner_size)
    bondvol = AssocParam{Float64}("epsilon assoc",components,compval_bondvol,comp_sites.sites,String[],String[])
    epsilon_assoc = AssocParam{Float64}("bondvol",components,compval_epsilon_assoc,comp_sites.sites,String[],String[])

    gcparams = GCsPCSAFTParam(gc_mw, gc_segment, gc_msigma3, gc_mepsilon, gc_epsilon_assoc,gc_bondvol)
    params = PCSAFTParam(mw, segment, sigma, epsilon, epsilon_assoc, bondvol)
    
    idmodel = init_model(idealmodel,components,ideal_userlocations,verbose)

    references = ["10.1021/ie020753p"]
    pc = sPCSAFT(params,comp_sites,idmodel;ideal_userlocations, references, verbose, assoc_options)
    model = GCsPCSAFT(components, groups, sites, gcparams,idmodel,pc, assoc_options, references)
    return model
end

function a_res(model::GCsPCSAFTModel,V,T,z)
    return a_res(model.pcsaftmodel,V,T,z)
end

JobackGroups = [raw"[CX4H3]" "-CH3";
raw"[!R;CX4H2]" "-CH2-";
raw"[!R;CX4H]" ">CH-";
raw"[!R;CX4H0]" ">C<";
raw"[CX3H2][CX3H1]" "CH2=CH-";
raw"[CX3H1][CX3H1]" "-CH=CH-";
raw"[$([!R;#6X3H0]);!$([!R;#6X3H0]=[#8])]" "=C<";
raw"[$([CX2H0](=*)=*)]" "=C=";
raw"[$([CX2H1]#[!#7])]" "CH";
raw"[$([CX2H0]#[!#7])]" "C";
raw"[R;CX4H2]" "ring-CH2-";
raw"[R;CX4H]" "ring>CH-";
raw"[R;CX4H0]" "ring>C<";
raw"[R;CX3H1,cX3H1]" "ring=CH-";
raw"[$([R;#6X3H0]);!$([R;#6X3H0]=[#8])]" "ring=C<";
raw"[F]" "-F";
raw"[Cl]" "-Cl";
raw"[Br]" "-Br";
raw"[I]" "-I";
raw"[OX2H;!$([OX2H]-[#6]=[O]);!$([OX2H]-a)]" "-OH (alcohol)";
raw"[O;H1;$(O-!@c)]" "-OH (phenol)";
raw"[OX2H0;!R;!$([OX2H0]-[#6]=[#8])]" "-O- (non-ring)";
raw"[#8X2H0;R;!$([#8X2H0]~[#6]=[#8])]" "-O- (ring)";
raw"[$([CX3H0](=[OX1]));!$([CX3](=[OX1])-[OX2]);!R]=O" ">C=O (non-ring)";
raw"[$([#6X3H0](=[OX1]));!$([#6X3](=[#8X1])~[#8X2]);R]=O" ">C=O (ring)";
raw"[CH;D2](=O)" "O=CH- (aldehyde)";
raw"[OX2H]-[C]=O" "-COOH (acid)";
raw"[#6X3H0;!$([#6X3H0](~O)(~O)(~O))](=[#8X1])[#8X2H0]" "-COO- (ester)";
raw"[OX1H0;!$([OX1H0]~[#6X3]);!$([OX1H0]~[#7X3]~[#8])]" "=O (other than above)";
raw"[NX3H2]" "-NH2";
raw"[NX3H1;!R]" ">NH (non-ring)";
raw"[#7X3H1;R]" ">NH (ring)";
raw"[#7X3H0;!$([#7](~O)~O)]" ">N- (non-ring)";
raw"[#7X2H0;!R]" "-N= (non-ring)";
raw"[#7X2H0;R]" "-N= (ring)";
raw"[#7X2H1]" "=NH";
raw"[#6X2]#[#7X1H0]" "-CN";
raw"[$([#7X3,#7X3+][!#8])](=[O])~[O-]" "-NO2";
raw"[SX2H]" "-SH";
raw"[#16X2H0;!R]" "-S- (non-ring)";
raw"[#16X2H0;R]" "-S- (ring)"]