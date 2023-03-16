# julia 1.6
#
# interface for parsing PGLib ACOPF data with linear costs to POP 
# (no generation power variables) 
#
# adapted from the TSSOS package https://github.com/wangjie212/TSSOS

using PowerModels
using MomentHierarchy
MH = MomentHierarchy


function normalize_coeff(coe)
    mc = maximum(abs.(coe))
    return coe./mc
end

function move_zero!(supp,coe)
    ind=[abs(coe[i])>=1e-8 for i=1:length(coe)]
    return supp[ind],coe[ind]
end

function fl_sum(vector)
    return mapreduce(x->x, +, vector, init = 0.0)
end

function bfind(A, l, a)
    if l==0
        return 0
    end
    low=1
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        if ndims(A)==2
            temp=A[:, mid]
        else
            temp=A[mid]
        end
        if temp==a
           return mid
        elseif temp<a
           low=mid+1
        else
           high=mid-1
        end
    end
    return 0
end

function parse_opf_linear_costs_to_pop(data::Dict{String, Any}; normal=true)

    silence()
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:it][pm_it_sym][:nw][0] 

    nbus=length(ref[:bus])
    nb=length(ref[:branch])
    ng=length(ref[:gen])
    
    n=2*nbus
    m=2*nbus+4*ng
    numeq = 2*(nbus-ng) + length(keys(ref[:ref_buses]))
    m = m + numeq
    
    startpoint=zeros(n)

    supp=Vector{Vector{Vector{UInt16}}}(undef, m+1)
    coe=Vector{Vector{Float64}}(undef, m+1)

    supp_pf = Vector{Vector{Vector{UInt16}}}(undef, n)
    coe_pf = Vector{Vector{Float64}}(undef, n)

    gens=collect(keys(ref[:gen]))
    sort!(gens)

    # minimal sparsity sets
    sets = Vector{Vector{UInt16}}()

    # objective function
    nc=ng+1
    coe[1]=Vector{Float64}(undef, 1)
    supp[1]=Vector{Vector{UInt16}}(undef, 1)
    coe[1][1]=sum(gen["cost"][3] for (i,gen) in ref[:gen])
    supp[1][1]=[]

    # move on
    
    k=2
    k_pf = 1

    bus=collect(keys(ref[:bus]))
    sort!(bus)

    # voltage magnitude constraints
    for i=1:nbus
        supp[k]=[[], [i;i], [i+nbus;i+nbus]]
        coe[k]=[-ref[:bus][bus[i]]["vmin"]^2;1;1]
        supp[k+1]=[[], [i;i], [i+nbus;i+nbus]]
        coe[k+1]=[ref[:bus][bus[i]]["vmax"]^2;-1;-1]
        startpoint[i]=sqrt(ref[:bus][bus[i]]["vmin"]*ref[:bus][bus[i]]["vmax"]/2)
        startpoint[i+nbus]=sqrt(ref[:bus][bus[i]]["vmin"]*ref[:bus][bus[i]]["vmax"]/2)
        k+=2
    end

    ## power flow constraints ##
    
    for (r, i) in enumerate(bus)
        
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        coe_pf[k_pf]=zeros(4*length(ref[:bus_arcs][i])+3)
        supp_pf[k_pf]=Vector{Vector{UInt16}}(undef, 4*length(ref[:bus_arcs][i])+3)
        
        supp_pf[k_pf][1:3]=[[], [r;r], [r+nbus;r+nbus]]
        coe_pf[k_pf][1]=fl_sum(load["pd"] for load in bus_loads) # L

        coe_pf[k_pf+1]=zeros(4*length(ref[:bus_arcs][i])+3)
        supp_pf[k_pf+1]=Vector{Vector{UInt16}}(undef, 4*length(ref[:bus_arcs][i])+3)
        
        supp_pf[k_pf+1][1:3]=[[], [r;r], [r+nbus;r+nbus]]
        coe_pf[k_pf+1][1]=fl_sum(load["qd"] for load in bus_loads)
        
        sgs=fl_sum(shunt["gs"] for shunt in bus_shunts)
        sbs=fl_sum(shunt["bs"] for shunt in bus_shunts)

        coe_pf[k_pf][2:3]=[sgs;sgs] # L + Shunt
        
        coe_pf[k_pf+1][2:3]=[-sbs;-sbs]

        j=1

        pf_set = [r, r+nbus]
        
        for flow in ref[:bus_arcs][i]
            branch=ref[:branch][flow[1]]
            vr = bfind(bus,nbus,branch["f_bus"])
            vt = bfind(bus,nbus,branch["t_bus"])
            srt=sort([vr;vt])
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]
            tm = branch["tap"]
            a1=(g+g_fr)/tm^2
            b1=-(b+b_fr)/tm^2
            c1=(-g*tr+b*ti)/tm^2
            d1=(b*tr+g*ti)/tm^2
            a2=g+g_to
            b2=-(b+b_to)
            c2=-(g*tr+b*ti)/tm^2
            d2=-(-b*tr+g*ti)/tm^2

            supp_pf[k_pf][j+3:j+6]=[srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus]
            supp_pf[k_pf+1][j+3:j+6]=[srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus]

            # check branch orientation 

            if vr==r
                coe_pf[k_pf][2:3].+=a1
                coe_pf[k_pf][j+3:j+6]=[c1;-d1;d1;c1]
                coe_pf[k_pf+1][2:3].+=b1
                coe_pf[k_pf+1][j+3:j+6]=[d1;c1;-c1;d1]

                append!(pf_set, [vt, vt+nbus])

            else
                coe_pf[k_pf][2:3].+=a2
                coe_pf[k_pf][j+3:j+6]=[c2;d2;-d2;c2]
                coe_pf[k_pf+1][2:3].+=b2
                coe_pf[k_pf+1][j+3:j+6]=[d2;-c2;c2;d2]

                append!(pf_set, [vr, vr+nbus])

            end

            j+=4

        end

        k_pf += 2

        # minimal sparsity sets

        push!(sets, pf_set)

    end # L + Shunt + sum S_ij

    # ready !

    no_gen = Int64[]

    for (r, i) in enumerate(bus)

        pf_id = 2*r-1

        if !isempty(ref[:bus_gens][i])

            gen=ref[:gen][ref[:bus_gens][i]...] # ???

            # add cost
            append!(coe[1], coe_pf[pf_id]*gen["cost"][2])
            append!(supp[1], supp_pf[pf_id])

            # bounds on Sg

            # min P/Q
            coe[k] = [coe_pf[pf_id]..., -gen["pmin"]]
            supp[k] = [supp_pf[pf_id]..., []]

            coe[k+1] = [coe_pf[pf_id+1]..., -gen["qmin"]]
            supp[k+1] = [supp_pf[pf_id+1]..., []]

            if normal==true
                coe[k]=normalize_coeff(coe[k])
                coe[k+1]=normalize_coeff(coe[k+1])
            end

            supp[k],coe[k]=move_zero!(supp[k],coe[k])
            supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])

            k += 2

            # max P/Q
            coe[k] = [-1*coe_pf[pf_id]..., gen["pmax"]] 
            supp[k] = [supp_pf[pf_id]..., []]

            coe[k+1] = [-1*coe_pf[pf_id+1]..., gen["qmax"]]
            supp[k+1] = [supp_pf[pf_id+1]..., []]

            if normal==true
                coe[k]=normalize_coeff(coe[k])
                coe[k+1]=normalize_coeff(coe[k+1])
            end

            supp[k],coe[k]=move_zero!(supp[k],coe[k])
            supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])

            k += 2

        else append!(no_gen, [pf_id, pf_id+1])

        end

    end

    supp[1],coe[1]=move_zero!(supp[1],coe[1])

    # power flow 
    for r in no_gen
        coe[k] = coe_pf[r]
        supp[k] = supp_pf[r]
        if normal==true
            coe[k]=normalize_coeff(coe[k])
        end
        supp[k], coe[k] = move_zero!(supp[k], coe[k])
        k += 1
    end

    # reference voltage
    for key in keys(ref[:ref_buses])
        i=bfind(bus,nbus,key)
        supp[k]=[[i+nbus;i+nbus]]
        coe[k]=[1]
        k+=1
    end

    ####### assembling the model #########

    if normal == true
        max_coefficient = maximum(abs.(coe[1]))
        coe[1] = coe[1] ./ max_coefficient
    else
        max_coefficient = nothing
    end

    objective = MH.SparsePolynomial(supp[1], coe[1]) 
    inequality_constraints = [MH.SparsePolynomial(supp[i], coe[i]) for i in 2:m-numeq+1]
    equality_constraints = [MH.SparsePolynomial(supp[i], coe[i]) for i in m-numeq+2:m+1]

    return MH.POP(objective, n, inequality_constraints, equality_constraints), startpoint, sets, max_coefficient

end