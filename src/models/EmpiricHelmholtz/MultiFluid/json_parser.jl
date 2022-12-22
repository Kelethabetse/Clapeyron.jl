struct IdealHelmholtzLead
    active::Bool
    a1::Float64
    a2::Float64
end

@inline function eval_f0(data::IdealHelmholtzLead,δ,τ)
    if data.active
    return log(δ) + data.a1 + data.a2*τ
    else
        return zero(δ+τ)
    end
end

struct IdealHelmholtzPlanckEinstein
    active::Bool
    n::Vector{Float64}
    t::Vector{Float64}
end

@inline function eval_f0(data::IdealHelmholtzPlanckEinstein,δ,τ)
    res = zero(δ+τ)
    if data.active
        n,t = data.n,data.t
        @inbounds for i in 1:length(n)
            res += n[i]*log1p(-exp(-t[i]*τ))
        end
    end
    return res
end

struct IdealHelmholtzLogTau
    active::Bool
    a::Float64
end

@inline function eval_f0(data::IdealHelmholtzLogTau,δ,τ)
    res = zero(τ)
    if data.active
        res += data.a*log(τ)
    end
    return res
end

@inline function eval_f0(data::IdealHelmholtzPlanckEinstein,δ,τ)
    res = zero(δ+τ)
    if data.active
        n,t = data.n,data.t
        @inbounds for i in 1:length(n)
            res += n[i]*log1p(-exp(-t[i]*τ))
        end
    end
    return res
end

function build_ideal(ideal_dict)
    _lead = get(ideal_dict,"IdealGasHelmholtzLead",nothing)
    if isnothing(_lead)
        lead = IdealHelmholtzLead(false,0.0,0.0)
    else
        lead = IdealHelmholtzLead(true,_lead[:a1],_lead[:a2])
    end

    _logtau = get(ideal_dict,"IdealGasHelmholtzLogTau",nothing)
    if isnothing(_lead)
        logtau = IdealHelmholtzLogTau(false,0.0)
    else
        logtau = IdealHelmholtzLogTau(true,_logtau[:a])
    end

    _pe1 = get(ideal_dict,"IdealGasHelmholtzPlanckEinstein",nothing)
    if isnothing(_pe1)
        pe1 = IdealHelmholtzPlanckEinstein(false,Float64[],Float64[])
    else
        pe1 = IdealHelmholtzPlanckEinstein(true,copy(_pe1[:n]),copy(_pe1[:t]))
    end
    #todo:flattening everything
    return lead,logtau,pe1
end


function parse_multifluid(path,type = :single)
    fullpath = only(flattenfilepaths(String[],String[path]))
    #@show fullpath
    json_string = read(fullpath, String) #TODO: see why JSON3.read(filepath) fails
    json = JSON3.read(json_string)

    eos_data = first(json[:EOS])
    comp_states = json[:STATES]
    comp_info = json[:INFO] # TODO: parse aliases to get a working list to reference
    
    
    #ideal part
    eos_ideal = eos_data[:alpha0]
    ideal_dict = Dict{String,Any}()
    for idealpart in eos_ideal
        ideal_dict[idealpart[:type]] = idealpart
    end
    id = build_ideal(ideal_dict)
    
    return json
end