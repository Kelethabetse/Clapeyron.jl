"""
    PSRKUNIFAC(components::Vector{String};
    puremodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
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

Modified [UNIFAC](@ref) (Dortmund) implementation, with parameters tuned to the Predictive Soave-Redlich-Kwong (PSRK) EoS.

## References
1. Fredenslund, A., Gmehling, J., Michelsen, M. L., Rasmussen, P., & Prausnitz, J. M. (1977). Computerized design of multicomponent distillation columns using the UNIFAC group contribution method for calculation of activity coefficients. Industrial & Engineering Chemistry Process Design and Development, 16(4), 450–462. [doi:10.1021/i260064a004](https://doi.org/10.1021/i260064a004)
2. Weidlich, U.; Gmehling, J. A modified UNIFAC model. 1. Prediction of VLE, hE, and.gamma..infin. Ind. Eng. Chem. Res. 1987, 26, 1372–1381.
3. Horstmann, S., Jabłoniec, A., Krafczyk, J., Fischer, K., & Gmehling, J. (2005). PSRK group contribution equation of state: comprehensive revision and extension IV, including critical constants and α-function parameters for 1000 components. Fluid Phase Equilibria, 227(2), 157–164. [doi:10.1016/j.fluid.2004.11.002](https://doi.org/10.1016/j.fluid.2004.11.002)"
"""
function PSRKUNIFAC(components::Vector{String};
    puremodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false, kwargs...)

    PRSK_userlocations = vcat("@REMOVEDEFAULTS","@DB/Activity/UNIFAC/PSRK",userlocations)
    PRSK_group_userlocations = vcat("@REMOVEDEFAULTS","@DB/Activity/UNIFAC/PSRK/PSRK_groups.csv",group_userlocations)
    model =  UNIFAC(components,
    puremodel = puremodel,
    userlocations = PRSK_userlocations,
    group_userlocations = PRSK_group_userlocations,
    pure_userlocations = pure_userlocations,
    verbose = verbose
    )
    setreferences!(model,String["10.1021/i260064a004","10.1016/j.fluid.2004.11.002"])
    return model
end

export PSRKUNIFAC

PSRKGroups = [raw"[CX4;H3]" "CH3";
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
raw"[CX3H2]=[CX3H2]" "H2C=CH2";
raw"[CX2;H1]#[CX2;H1]" "CH=-CH";
raw"[NX3H3]" "NH3";
raw"[C-]#[O+]" "CO";
raw"[HH]" "H2";
raw"[SH2]" "H2S";
raw"N#N" "N2";
raw"[ArX0]" "AR";
raw"[CX2H0](=[OX1H0])=[OX1H0]" "CO2";
raw"[CX4H4]" "CH4";
raw"[OX1H0]=[OX1H0]" "O2";
raw"[2H][2H]" "D2";
raw"[OX1H0]=[SX2H0]=[OX1H0]" "SO2";
raw"[NX1H0]=[OX1H0]" "NO";
raw"[NX1H0]#[N+X2H0][O-X1H0]" "N2O";
raw"[FX1H0][SX6H0]([FX1H0])([FX1H0])([FX1H0])([FX1H0])[FX1H0]" "SF6";
raw"[HeX0H0]" "HE";
raw"[NeX0]" "NE";
raw"[KrX0]" "KR";
raw"[XeX0]" "XE";
raw"[FX1H1]" "HF";
raw"[ClX1H1]" "HCL";
raw"[BrX1H1]" "HBR";
raw"[IX1H1]" "HI";
raw"[CX2H0](=[OX1H0])=[SX1H0]" "COS";
raw"[CX4H1][SX2H1]" "CHSH";
raw"[CX4H0][SX2H1]" "CSH";
raw"[CX4H2]1[CX4H1][OX2H0]1" "H2COCH";
raw"[CX4H1]1[CX4H1][OX2H0]1" "HCOCH";
raw"[CX4H1]1[CX4H0][OX2H0]1" "HCOC";
raw"[CX4H2]1[CX4H2][OX2H0]1" "H2COCH2";
raw"[CX4H2]1[CX4H0][OX2H0]1" "H2COC";
raw"[CX4H0]1[CX4H0][OX2H0]1" "COC";
raw"[F][F]" "F2";
raw"[Cl][Cl]" "CL2";
raw"[Br][Br]" "BR2";
raw"[CX2H1]#[NX1H0]" "HCN";
raw"[OX1H0][NX2H0]=[OX1H0]" "NO2";
raw"[CX4H0]([FX1])([FX1])([FX1])[FX1]" "CF4";
raw"[O-X1H0][O+X2H0]=[OX1H0]" "O3";
raw"[NX2H0](=[OX1H0])[ClX1H0]" "CLNO"]