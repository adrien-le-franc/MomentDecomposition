# julia 1.6
#
# interface for parsing PGLib ACOPF data to POP 
# adapted from the TSSOS package https://github.com/wangjie212/TSSOS

using PowerModels
using MomentHierarchy
MH = MomentHierarchy


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

function parse_opf_to_pop(data::Dict{String, Any}; AngleCons=false, 
    LineLimit=false, decomposition=false, nlp=false, n_max::Int64=16, quartic_v=false)

    silence()
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:it][pm_it_sym][:nw][0] 
    nbus=length(ref[:bus])
    nb=length(ref[:branch])
    ng=length(ref[:gen])
    n=2*nbus+2*ng
    m=4*nbus+2*ng+length(keys(ref[:ref_buses]))
    if quartic_v == true
        m -= nbus
    end
    if AngleCons==true
        m+=2*nb
    end
    if LineLimit==true||LineLimit=="relax"
        m+=2*nb
    end
    numeq=2*nbus+length(keys(ref[:ref_buses]))
    startpoint=zeros(n)
    supp=Vector{Vector{Vector{UInt16}}}(undef, m+1)
    coe=Vector{Vector{Float64}}(undef, m+1)

    println("initially $(n) variables")

    # adaptative power flow equations
    new_vars = 0
    new_cons = 0
    new_supp = Vector{Vector{Vector{UInt16}}}()
    new_coe = Vector{Vector{Float64}}()

    # minimal sparsity sets
    sets = Vector{Vector{UInt16}}()

    ## objective function ##

    gens=collect(keys(ref[:gen]))
    sort!(gens)

    nc=2*ng+1
    coe[1]=Vector{Float64}(undef, nc)
    supp[1]=Vector{Vector{UInt16}}(undef, nc)
    coe[1][1]=sum(gen["cost"][3] for (i,gen) in ref[:gen])
    supp[1][1]=[]
    for i=1:ng
        gen=ref[:gen][gens[i]]
        coe[1][2*(i-1)+2:2*(i-1)+3]=[gen["cost"][2];gen["cost"][1]]
        supp[1][2*(i-1)+2:2*(i-1)+3]=[[2*nbus+i], [2*nbus+i;2*nbus+i]]
    end
    supp[1],coe[1]=move_zero!(supp[1],coe[1])
    k=2

    bounds = Dict{UInt16, Vector{Float64}}()

    ## voltage magnitude constraints ##

    bus=collect(keys(ref[:bus]))
    sort!(bus)
    
    for i=1:nbus

        if quartic_v

            prod_v = (ref[:bus][bus[i]]["vmin"]^2)*(ref[:bus][bus[i]]["vmax"]^2)
            sum_v = ref[:bus][bus[i]]["vmin"]^2 + ref[:bus][bus[i]]["vmax"]^2

            supp[k]=[[], [i;i], [i+nbus;i+nbus], [i;i; i+nbus;i+nbus], [i;i;i;i], [i+nbus;i+nbus;i+nbus;i+nbus]]
            coe[k]=[-1*prod_v; sum_v; sum_v; -2.; -1.; -1.]            
            k+=1

        else

            supp[k]=[[], [i;i], [i+nbus;i+nbus]]
            coe[k]=[-ref[:bus][bus[i]]["vmin"]^2;1;1]
            supp[k+1]=[[], [i;i], [i+nbus;i+nbus]]
            coe[k+1]=[ref[:bus][bus[i]]["vmax"]^2;-1;-1]
            k+=2

        end

        startpoint[i]=sqrt(ref[:bus][bus[i]]["vmin"]*ref[:bus][bus[i]]["vmax"]/2)
        startpoint[i+nbus]=sqrt(ref[:bus][bus[i]]["vmin"]*ref[:bus][bus[i]]["vmax"]/2)

        bounds[convert(UInt16, i)] = [ref[:bus][bus[i]]["vmax"], -ref[:bus][bus[i]]["vmax"]]
        bounds[convert(UInt16, i+nbus)] = [ref[:bus][bus[i]]["vmax"], -ref[:bus][bus[i]]["vmax"]]

    end

    ## angle and line limit constraints ##

    if AngleCons==true||LineLimit==true||LineLimit=="relax"
        for (i, branch) in ref[:branch]
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]
            tm = branch["tap"]
            vr = bfind(bus,nbus,branch["f_bus"])
            vt = bfind(bus,nbus,branch["t_bus"])
            srt=sort([vr;vt])
            dsrt=sort([vr;vr;vt;vt])
            ab1=(g+g_fr)^2+(b+b_fr)^2
            cd1=(-g*tr+b*ti)^2+(b*tr+g*ti)^2
            acbd1=(g+g_fr)*(-g*tr+b*ti)-(b+b_fr)*(b*tr+g*ti)
            bcad1=-(b+b_fr)*(-g*tr+b*ti)-(g+g_fr)*(b*tr+g*ti)
            ab2=(g+g_to)^2*tm^4+(b+b_to)^2*tm^4
            cd2=(g*tr+b*ti)^2+(-b*tr+g*ti)^2
            acbd2=-(g+g_to)*tm^2*(g*tr+b*ti)+(b+b_to)*tm^2*(-b*tr+g*ti)
            bcad2=(b+b_to)*tm^2*(g*tr+b*ti)+(g+g_to)*tm^2*(-b*tr+g*ti)

            # angle differences
            if AngleCons==true
                coe[k]=[tan(branch["angmax"]);tan(branch["angmax"]);-1;1]
                supp[k]=[srt, srt.+nbus, [vt;vr+nbus], [vr;vt+nbus]]
                coe[k+1]=[1;-1;-tan(branch["angmin"]);-tan(branch["angmin"])]
                supp[k+1]=[[vt;vr+nbus], [vr;vt+nbus], srt, srt.+nbus]
                k+=2
            end

            # thermal limits
            if LineLimit==true
                coe[k]=[branch["rate_a"]^2*tm^4;-ab1;-2*ab1;-ab1;-cd1;-cd1;-cd1;-cd1;-2*acbd1;-2*acbd1;-2*acbd1;-2*acbd1;-2*bcad1;2*bcad1;-2*bcad1;2*bcad1]
                supp[k]=[[], [vr;vr;vr;vr], [vr;vr;vr+nbus;vr+nbus], [vr+nbus;vr+nbus;vr+nbus;vr+nbus], [vt;vt;vr+nbus;vr+nbus], [vr;vr;vt+nbus;vt+nbus],
                dsrt, dsrt.+nbus, sort([vr;vr;vr;vt]), [vr;vr;sort([vr+nbus;vt+nbus])], [srt;vr+nbus;vr+nbus], sort([vr;vr;vr;vt]).+nbus, [sort([vr;vr;vt]);vr+nbus],
                [vr;vr;vr;vt+nbus], [vt;vr+nbus;vr+nbus;vr+nbus], [vr;sort([vr+nbus;vr+nbus;vt+nbus])]]
                supp[k],coe[k]=move_zero!(supp[k],coe[k])
                coe[k+1]=[branch["rate_a"]^2*tm^4;-ab2;-2*ab2;-ab2;-cd2;-cd2;-cd2;-cd2;-2*acbd2;-2*acbd2;-2*acbd2;-2*acbd2;-2*bcad2;2*bcad2;-2*bcad2;2*bcad2]
                supp[k+1]=[[], [vt;vt;vt;vt], [vt;vt;vt+nbus;vt+nbus], [vt+nbus;vt+nbus;vt+nbus;vt+nbus], [vr;vr;vt+nbus;vt+nbus], [vt;vt;vr+nbus;vr+nbus],
                dsrt, dsrt.+nbus, sort([vt;vt;vt;vr]), [vt;vt;sort([vt+nbus;vr+nbus])], [srt;vt+nbus;vt+nbus], sort([vt;vt;vt;vr]).+nbus, [sort([vt;vt;vr]);vt+nbus],
                [vt;vt;vt;vr+nbus], [vr;vt+nbus;vt+nbus;vt+nbus], [vt;sort([vt+nbus;vt+nbus;vr+nbus])]]
                supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
                k+=2
            elseif LineLimit=="relax"
                mvr=ref[:bus][bus[vr]]["vmin"]^2
                mvt=ref[:bus][bus[vt]]["vmin"]^2
                coe[k]=[branch["rate_a"]^2*tm^4/mvr;-ab1;-ab1;-cd1;-cd1;-2*acbd1;2*bcad1;-2*bcad1;-2*acbd1]
                supp[k]=[[], [vr;vr], [vr+nbus;vr+nbus], [vt;vt], [vt+nbus;vt+nbus], srt, [vr;vt+nbus], [vt;vr+nbus], srt.+nbus]
                supp[k],coe[k]=move_zero!(supp[k],coe[k])
                coe[k+1]=[branch["rate_a"]^2*tm^4/mvt;-ab2;-ab2;-cd2;-cd2;-2*acbd2;2*bcad2;-2*bcad2;-2*acbd2]
                supp[k+1]=[[], [vt;vt], [vt+nbus;vt+nbus], [vr;vr], [vr+nbus;vr+nbus], srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus]
                supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
                k+=2
            end
        end
    end

    ## power generation bound ##

    zero_pgen=UInt16[]
    zero_qgen=UInt16[]

    for i=1:ng

        gen=ref[:gen][gens[i]]

        if gen["pmax"] - gen["pmin"] >= 1e-6

            coe[k]=[-gen["pmin"]*gen["pmax"];gen["pmin"]+gen["pmax"];-1]
            supp[k]=[[], [i+2*nbus], [i+2*nbus;i+2*nbus]]
            supp[k],coe[k]=move_zero!(supp[k],coe[k])
            startpoint[i+2*nbus]=(gen["pmin"]+gen["pmax"])/2
            k+=1
        
        else
            
            push!(zero_pgen, i)
            startpoint[i+2*nbus]=0
            numeq+=1
        
        end

        if gen["qmax"] - gen["qmin"] >= 1e-6

            coe[k]=[-gen["qmin"]*gen["qmax"];gen["qmin"]+gen["qmax"];-1]
            supp[k]=[[], [i+2*nbus+ng], [i+2*nbus+ng;i+2*nbus+ng]]
            supp[k],coe[k]=move_zero!(supp[k],coe[k])
            startpoint[i+2*nbus+ng]=0
            k+=1
        
        else
            
            push!(zero_qgen, i)
            startpoint[i+2*nbus+ng]=0
            numeq+=1
        
        end

        bounds[convert(UInt16, i+2*nbus)] = [gen["pmax"], gen["pmin"]]
        bounds[convert(UInt16, i+2*nbus+ng)] = [gen["qmax"], gen["qmin"]]

    end

    ## power flow constraints ##

    for (r, i) in enumerate(bus)

        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
        sgs=fl_sum(shunt["gs"] for shunt in bus_shunts)
        sbs=fl_sum(shunt["bs"] for shunt in bus_shunts)

        supp[k] = [[], [r;r], [r+nbus;r+nbus]]
        coe[k]=[fl_sum(load["pd"] for load in bus_loads), sgs, sgs]
        
        supp[k+1] = [[], [r;r], [r+nbus;r+nbus]]
        coe[k+1] = [fl_sum(load["qd"] for load in bus_loads), -sbs, -sbs] 
        
        # number of variables to write the pair of Real & Imaginary PF equations:
        n_gens = length(ref[:bus_gens][i])
        n_lines = length(ref[:bus_arcs][i])
        n_vars = 2 + 2*n_gens + 2*n_lines

        pf_set = [r, r+nbus]

        n_line_clusters = 0

        if n_vars > n_max

            #println("PF at bus $(r) - $(i) has $(n_vars) variables")
            #println("$(n_gens) gens - $(n_lines) lines")

            n_line_clusters = ceil(Int, n_lines / (floor(Int, n_max / 2) - 2))
            max_lines_per_cluster = ceil(Int, n_lines / n_line_clusters)

            #println("reducing to $(n_line_clusters) clusters of max $(2*max_lines_per_cluster + 4) variables")

            new_sets = Vector{Vector{UInt16}}()

            var = 1

            for _ in 1:n_line_clusters

                # update old PF equations
                push!(supp[k], [n + new_vars + var])
                push!(coe[k], -1.)
                push!(supp[k+1], [n + new_vars + var+1])
                push!(coe[k+1], -1.)

                # initialize new PF equations
                push!(new_supp, [[n + new_vars + var], [r;r], [r+nbus;r+nbus]])
                push!(new_coe, [1., 0., 0.])
                push!(new_supp, [[n + new_vars + var+1], [r;r], [r+nbus;r+nbus]])
                push!(new_coe, [1., 0., 0.])

                # new variables in set
                append!(pf_set, [n + new_vars + var, n + new_vars + var+1])

                # initialize new sets
                push!(new_sets, [r, r+nbus, n + new_vars + var, n + new_vars + var+1])

                # initialize new variables bounds
                bounds[convert(UInt16, n + new_vars + var)] = [0., 0.]
                bounds[convert(UInt16, n + new_vars + var+1)] = [0., 0.]

                var += 2                                       

            end
       
            n_lines_in_cluster = 0
            n_cluster = 1
            cons = 1
            var = 1

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

                append!(new_supp[new_cons + cons], [srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus])
                append!(new_supp[new_cons + cons+1], [srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus])

                bounds[convert(UInt16, n + new_vars + var)] .+= [branch["rate_a"]*tm^2, -branch["rate_a"]*tm^2]
                bounds[convert(UInt16, n + new_vars + var+1)] .+= [branch["rate_a"]*tm^2, -branch["rate_a"]*tm^2]

                if vr==r

                    new_coe[new_cons + cons][2:3].+=a1
                    append!(new_coe[new_cons + cons], [c1;-d1;d1;c1])
                    new_coe[new_cons + cons+1][2:3].+=b1
                    append!(new_coe[new_cons + cons+1], [d1;c1;-c1;d1])

                    append!(new_sets[n_cluster], [vt, vt+nbus])

                else

                    new_coe[new_cons + cons][2:3].+=a2
                    append!(new_coe[new_cons + cons], [c2;d2;-d2;c2])
                    new_coe[new_cons + cons+1][2:3].+=b2
                    append!(new_coe[new_cons + cons+1], [d2;-c2;c2;d2])

                    append!(new_sets[n_cluster], [vr, vr+nbus])

                end

                n_lines_in_cluster += 1

                if n_lines_in_cluster == max_lines_per_cluster

                    cons += 2
                    var += 2
                    n_cluster += 1
                    n_lines_in_cluster = 0

                end

            end

        new_vars += 2*n_line_clusters
        new_cons += 2*n_line_clusters

        append!(sets, new_sets)

        else

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

                append!(supp[k], [srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus])
                append!(supp[k+1], [srt, [vt;vr+nbus], [vr;vt+nbus], srt.+nbus])

                if vr==r

                    coe[k][2:3].+=a1
                    append!(coe[k], [c1;-d1;d1;c1])
                    coe[k+1][2:3].+=b1
                    append!(coe[k+1], [d1;c1;-c1;d1])

                    append!(pf_set, [vt, vt+nbus])

                else
                    coe[k][2:3].+=a2
                    append!(coe[k], [c2;d2;-d2;c2])
                    coe[k+1][2:3].+=b2
                    append!(coe[k+1], [d2;-c2;c2;d2])

                    append!(pf_set, [vr, vr+nbus])

                end
            end

        end

        # add gens

        n_vars = 2 + 2*n_gens + 2*n_line_clusters

        if n_vars > n_max && n_gens > 1

            n_gen_clusters = ceil(Int, n_gens / (floor(Int, n_max / 2) - 2))
            max_gens_per_cluster = ceil(Int, n_gens / n_gen_clusters)

            new_sets = Vector{Vector{UInt16}}()

            var = 1

            for _ in 1:n_gen_clusters

                # update old PF equations
                push!(supp[k], [n + new_vars + var])
                push!(coe[k], -1.)
                push!(supp[k+1], [n + new_vars + var+1])
                push!(coe[k+1], -1.)

                # initialize new gen equations
                push!(new_supp, [[n + new_vars + var]])
                push!(new_coe, [1.])
                push!(new_supp, [[n + new_vars + var+1]])
                push!(new_coe, [1.])

                # new variables in set
                append!(pf_set, [n + new_vars + var, n + new_vars + var+1])

                # initialize new sets
                push!(new_sets, [r, r+nbus, n + new_vars + var, n + new_vars + var+1])

                # initialize new variables bounds
                bounds[convert(UInt16, n + new_vars + var)] = [0., 0.]
                bounds[convert(UInt16, n + new_vars + var+1)] = [0., 0.]

                var += 2
            
            end

            n_gens_in_cluster = 0
            n_cluster = 1
            cons = 1
            var = 1

            for gen_id in ref[:bus_gens][i]

                gen=bfind(gens, ng, gen_id)

                # update new gen equations
                push!(new_supp[new_cons + cons], [gen+2*nbus]) 
                push!(new_coe[new_cons + cons], -1)
                push!(new_supp[new_cons + cons+1], [gen+2*nbus+ng]) 
                push!(new_coe[new_cons + cons+1], -1)

                # update new sets
                append!(new_sets[n_cluster], [gen+2*nbus, gen+2*nbus+ng])

                # update new bounds
                gen_info = ref[:gen][gen_id]
                bounds[convert(UInt16, n + new_vars + var)] .+= [gen_info["pmax"], gen_info["pmin"]]
                bounds[convert(UInt16, n + new_vars + var+1)] .+= [gen_info["qmax"], gen_info["qmin"]]

                # flags

                n_gens_in_cluster += 1

                if n_gens_in_cluster == max_gens_per_cluster

                    cons += 2
                    var += 2
                    n_cluster += 1
                    n_gens_in_cluster = 0

                end

            end

            new_vars += 2*n_gen_clusters
            new_cons += 2*n_gen_clusters

            append!(sets, new_sets)

            println("added $(n_gen_clusters) gen variables")

        else
            
            if !isempty(ref[:bus_gens][i]) 
                for gen_id in ref[:bus_gens][i]
                    gen=bfind(gens, ng, gen_id)
                    push!(supp[k], [gen+2*nbus])
                    push!(coe[k], -1)
                    push!(supp[k+1], [gen+2*nbus+ng])
                    push!(coe[k+1], -1)

                    append!(pf_set, [gen+2*nbus, gen+2*nbus+ng])

                end
            end

        end

        # minimal sparsity sets

        push!(sets, pf_set)

        # zeros

        supp[k],coe[k]=move_zero!(supp[k],coe[k])
        supp[k+1],coe[k+1]=move_zero!(supp[k+1],coe[k+1])
        
        # next !

        k+=2

    end

    ## reference voltage constraints ##

    for key in keys(ref[:ref_buses])
        i=bfind(bus,nbus,key)
        supp[k]=[[i+nbus;i+nbus]]
        coe[k]=[1]
        k+=1
    end

    # zero power generation

    for i in zero_pgen
        supp[k]=[[i+2*nbus]]
        coe[k]=[1]
        
        bounds[convert(UInt16, i+2*nbus)] = [1., 0.]

        k+=1
    end

    for i in zero_qgen
        supp[k]=[[i+2*nbus+ng]]
        coe[k]=[1]
        
        bounds[convert(UInt16, i+2*nbus+ng)] = [1., 0.]

        k+=1
    end

    ####### assembling the model #########

    objective = MH.SparsePolynomial(supp[1], coe[1]) 
    inequality_constraints = [MH.SparsePolynomial(supp[i], coe[i]) for i in 2:m-numeq+1]
    equality_constraints = [MH.SparsePolynomial(supp[i], coe[i]) for i in m-numeq+2:m+1]

    # add new eq
    for i in 1:new_cons
        new_supp[i], new_coe[i] = move_zero!(new_supp[i], new_coe[i])
        push!(equality_constraints, MH.SparsePolynomial(new_supp[i], new_coe[i]))
    end

    # add new ineq
    """
    for i in 1:new_vars
        var = n + new_vars
        coefficients = [-bounds[var][1]*bounds[var][2];bounds[var][1]+bounds[var][2];-1] 
        support = [[], [var], [var;var]]
        push!(inequality_constraints, MH.SparsePolynomial(support, coefficients)) 
    end
    """
    
    # add new vars
    n += new_vars

    println("added $(new_vars) variables")

    # add new startpoints

    new_x = zeros(new_vars)

    for i in 1:new_vars

        x_i = 0.

        for (support, coeff) in zip(new_supp[i][2:end], new_coe[i][2:end])

            if length(support) == 0
                x_i += coeff
            else
                x_i += coeff*prod(startpoint[k] for k in support)
            end

        end

        new_x[i] = -x_i

    end

    startpoint = vcat(startpoint, new_x)
 
    return MH.POP(objective, n, inequality_constraints, equality_constraints), startpoint, sets, bounds

end