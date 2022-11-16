function get_groups(component::String,groups::Array{String})
    res = search_chemical(component)
    mol = get_mol(res.smiles)
    mol_list = get_substruct_matches(mol,mol)

    queries = get_qmol.(groups[:,1])

    atoms = mol_list[1]["atoms"]
    group_list = []
    group_id = []
    group_occ_list = []
    atoms_list = []
    for i in 1:length(groups[:,1])    
        if !isempty(get_substruct_match(mol,queries[i]))
            smatch = get_substruct_matches(mol,queries[i])
            for j in 1:length(smatch)
                if isempty(atoms_list)
                    append!(group_list,[groups[i,1]])
                    append!(group_id,i)
                    append!(group_occ_list,1)
                    append!(atoms_list,smatch[j]["atoms"])
                else
                    if sum(smatch[j]["atoms"] .∈ [atoms_list])==0
                        append!(atoms_list,smatch[j]["atoms"])
                        if !(groups[i,1] in group_list)
                            append!(group_list,[groups[i,1]])
                            append!(group_id,i)
                            append!(group_occ_list,1)
                        else
                            group_occ_list[end] += 1
                        end
                    end
                end
            end
        end
    end

    if !(sum(atoms_list .∈ [atoms])==length(atoms))
        warning("Could not find all groups")
    end

    return (component,[groups[group_id[i],2] => group_occ_list[i] for i in 1:length(group_id)])
end

export get_groups