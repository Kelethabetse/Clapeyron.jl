
@newmodelsingleton LinearMixing MixingRule

"""
    LinearMixing <: MultiFluidDepartureModel
    LinearMixing(components;
    userlocations=String[],
    verbose=false)

## Input parameters
none

## Description
Linear mixing rule for MultiParameter EoS models:

```
τ = T̄/T
δ = V̄/V
V̄ = ∑xᵢVcⱼ
T̄ = ∑xᵢTcᵢ
```
"""
LinearMixing

function v_scale(model::MultiFluid,z,mixing::LinearMixing,∑z)
    Vc = model.params.Vc.values
    dot(Vc,z)/∑z
end

function T_scale(model::MultiFluid,z,mixing::LinearMixing,∑z)
    Tc = model.params.Tc.values
    return dot(Tc,z)/∑z
end

export LinearMixing