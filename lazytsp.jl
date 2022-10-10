using JuMP, Gurobi, GLPK
include("utils.jl")

function build_tsp_model(optimizer, d, n, lower = -Inf, upper = +Inf)

    @debug "build_tsp_model..."

    model = Model(optimizer)
    @variable(model, x[1:n, 1:n], Bin, Symmetric)
    @objective(model, Min, sum(d .* x)/2)

    @constraint(model, [i in 1:n], sum(x[i, :]) == 2)
    @constraint(model, [i in 1:n], x[i, i] == 0)
    if lower != -Inf
        @constraint(model, lower <= sum(d .* x)/2)
    end
    if upper != +Inf
        @constraint(model, sum(d .* x)/2 <= upper)
    end

    @debug "build_tsp_model...done"

    return model
end

function subtour(edges::Vector{Tuple{Int,Int}}, n)

    @debug "subtour..."

    shortest_subtour, unvisited = collect(1:n), Set(collect(1:n))
    found_subtour_cycle = false
    while !isempty(unvisited)
        this_cycle, neighbors = Int[], unvisited
        while !isempty(neighbors)
            current = pop!(neighbors)
            push!(this_cycle, current)
            if length(this_cycle) > 1
                pop!(unvisited, current)
            end
            neighbors =
                [j for (i, j) in edges if i == current && j in unvisited]
        end

        if length(this_cycle) < length(shortest_subtour)
            shortest_subtour = this_cycle
            found_subtour_cycle = true
        end

        if length(this_cycle) > 10 && !found_subtour_cycle
            @debug "subtour...done (1)"
            return shortest_subtour, this_cycle
        end

        found_subtour_cycle = false

    end

    @debug "subtour...done (2)"

    return shortest_subtour, shortest_subtour

end

function tasubtour(edges::Vector{Tuple{Int,Int}}, edges1::Vector{Tuple{Int,Int}}, n) # matching

    @debug "subtour..."
    new_graph = SimpleGraph(n)

	for (i, j) in edges
        add_edge!(new_graph, i, j)
	end

    deg = degree(new_graph)

    @info "degree: " deg
    index = zeros(Int64, nv(new_graph))
    k = 1
    for d in 1:nv(new_graph)
        if (deg[d] & 1) == 1
            index[k] = d
            k+=1
        end
    end
    @info "index odd degree: " index
    unvisited = Set([i for i in index if i != 0])
    while !isempty(unvisited)
        this_cycle, neighbors = Int[], unvisited
        while !isempty(neighbors)
            current = pop!(neighbors)
            push!(this_cycle, current)
            if length(this_cycle) > 1
                pop!(unvisited, current)
            end
            neighbors =
                [j for (i, j) in edges if i == current && j in unvisited]
        end
        
        comb = []
        if length(this_cycle) >= 3
            @info "ENCONTREI UM CICLO odd = ", this_cycle
            for atual in this_cycle
                neighbors =
                    [j for (i, j) in edges1 if i == atual && !(j in this_cycle)]
                if length(neighbors) > 0
                    push!(comb, (atual, neighbors))
                end
            end
            @info "ENCONTREI UM CICLO comb = ", comb
        end
        
        if length(comb) & 1 == 1 && length(comb) >= 3
            return true, this_cycle, comb
        end
    end

    @debug "subtour...done (2)"

    return false, [], []

end

function casubtour(edges::Vector{Tuple{Int,Int}}, edges1::Vector{Tuple{Int,Int}}, n) #comb

    @debug "subtour..."
    new_graph = SimpleGraph(n)

	for (i, j) in edges
        add_edge!(new_graph, i, j)
	end

    deg = degree(new_graph)

    @info "degree: " deg
    index = zeros(Int64, nv(new_graph))
    k = 1
    for d in 1:nv(new_graph)
        if (deg[d] & 1) == 1
            index[k] = d
            k+=1
        end
    end
    @info "index odd degree: " index
    unvisited = Set([i for i in index if i != 0])
    while !isempty(unvisited)
        this_cycle, neighbors = Int[], unvisited
        while !isempty(neighbors)
            current = pop!(neighbors)
            push!(this_cycle, current)
            if length(this_cycle) > 1
                pop!(unvisited, current)
            end
            neighbors =
                [j for (i, j) in edges if i == current && j in unvisited]
        end
        
        comb = []
        if length(this_cycle) >= 3
            @info "ENCONTREI UM CICLO odd = ", this_cycle
            for atual in this_cycle
                neighbors = [j for (i, j) in edges1 if i == atual && !(j in this_cycle)]
                neighborss = []
                for meudeus in neighbors
                    temp = [(meudeus, j) for (i, j) in edges if i == meudeus && j in unvisited]
                    if length(temp) > 0
                        push!(neighborss, temp)
                    end
                end
                if length(neighbors) > 0
                    push!(comb, (atual, neighbors, neighborss))
                end
            end
            @info "ENCONTREI UM CICLO comb = ", comb
        end
        
        # if length(this_cycle) > 10
        #     @debug "subtour...done (1)"
        #     return this_cycle
        # end
        if length(comb) & 1 == 1 || length(comb) >= 3
            return true, this_cycle, comb
        end
    end

    @debug "subtour...done (2)"

    return false, [], []

end


function selected_edges(x::Matrix{Float64}, n)
    return Tuple{Int,Int}[(i, j) for i in 1:n, j in 1:n if x[i, j] > 0.5]
end

function bselected_edges(x::Matrix{Float64}, n)
    return Tuple{Int,Int}[(i, j) for i in 1:n, j in 1:n if x[i, j] == 1]
end

function aselected_edges(x::Matrix{Float64}, n)
    return Tuple{Int,Int}[(i, j) for i in 1:n, j in 1:n if x[i, j] > 0 && x[i, j] < 1]
end

subtour(x::Matrix{Float64}) = subtour(selected_edges(x, size(x, 1)), size(x, 1))
subtour(x::AbstractMatrix{VariableRef}) = subtour(value.(x))

tasubtour(x::Matrix{Float64}, y::Matrix{Float64}) = tasubtour(aselected_edges(x, size(x, 1)), selected_edges(y, size(y, 1)), size(x, 1))
tasubtour(x::AbstractMatrix{VariableRef}, y::AbstractMatrix{VariableRef}) = tasubtour(value.(x), value.(y))

casubtour(x::Matrix{Float64}, y::Matrix{Float64}) = casubtour(aselected_edges(x, size(x, 1)), selected_edges(y, size(y, 1)), size(x, 1))
casubtour(x::AbstractMatrix{VariableRef}, y::AbstractMatrix{VariableRef}) = casubtour(value.(x), value.(y))

function plot(exp_id::String, edges::Vector{Tuple{Int,Int}}, n)

    new_graph = SimpleGraph(n)

	for (i, j) in edges
        add_edge!(new_graph, i, j)
	end

    draw(PNG(string("./work/",exp_id,"/new_graph_", get_next_id(),".png"), 45cm, 45cm), gplot(new_graph, nodelabel=1:nv(new_graph)))

    deg = degree(new_graph)

    @info "degree: " deg
    index = zeros(Int64, nv(new_graph))
    k = 1
    for d in 1:nv(new_graph)
        if (deg[d] & 1) == 1
            index[k] = d
            k+=1
        end
    end
    @info "index odd degree: " index


end

plot(exp_id::String, x::Matrix{Float64}) = plot(exp_id, selected_edges(x, size(x, 1)), size(x, 1))
plot(exp_id::String, x::AbstractMatrix{VariableRef}) = plot(exp_id, value.(x))

#
# Function lazy
# obtains the lagrangean_relaxation with christofides method 
#  
# Parameters
# exp_params -> experiment parameters
# Returns 
# void
function lazy(exp_params)

    exp_id = exp_params["exp_id"]

    save_step(exp_id,"lazy","start","algorithm")

    lower = parse(Float64, exp_params["lower_bound"])
    upper = parse(Float64, exp_params["upper_bound"])
    testdata = String(read(exp_params["testdatafile"]))
    optimizer = Gurobi.Optimizer
    if exp_params["strategy"] == "glpk"
        optimizer = GLPK.Optimizer 
    end

    xml_graph = parsexml(testdata)

    cost, n = graph_to_cost_matrix(xml_graph)

    total_iterations = 0
    experiment_start_time = get_time_in_ms()
    experiment_finish_time = get_time_in_ms()


    # lower, upper = lagrangean_relaxation(900, false, file)
    # println("lower = ", lower)
    # println("upper = ", upper)
    lazy_model = build_tsp_model(optimizer, cost, n, lower, upper)
    function subtour_elimination_callback(cb_data)

        @debug "subtour_elimination_callback..."

        total_iterations = total_iterations + 1

#         status = callback_node_status(cb_data, lazy_model)
#         if status != MOI.CALLBACK_NODE_STATUS_INTEGER
#             @debug "subtour_elimination_callback...done: only run at integer solution"
#             return  # Only run at integer solutions
#         end
        cycle, c = subtour(callback_value.(cb_data, lazy_model[:x]))
        shouldI, h, t = casubtour(callback_value.(cb_data, lazy_model[:x]), callback_value.(cb_data, lazy_model[:x]))

        @debug begin
            plot(exp_id,callback_value.(cb_data, lazy_model[:x]))            
        end

        if exp_params["rac_strategy"] in ["both combs"] && shouldI
            @info "generating combs..."
            # 2-matching inequalities
            # S = []
            # for i in 1:length(h)-1
            #     push!(S,(h[i], h[i+1]))
            # end
            # push!(S,(h[1], h[length(h)]))
            # @info "generating combs... S = ", S
            # T = []
            # @info "generating combs... t = ", t
            # ajuda = (length(t) & 1) == 1 ? 0 : 1
            # for i in 1:length(t)-ajuda
            #     @info "generating combs... t[i] = ", t[i]
            #     @info "generating combs... t[i][1] = ", t[i][1]
            #     @info "generating combs... t[i][2] = ", t[i][2]
            #     if length(t[i][2]) > 0
            #         push!(T,(t[i][1], t[i][2][1])) # t[i][1] = vemh
            #     end
            # end
            # @info "CRIAMOS A 2 MATCH - VITORIA"
            # @info "CRIAMOS A 2 MATCH - VITORIA len T = ", length(T)
            # @info "CRIAMOS A 2 MATCH - VITORIA T = ", T
            # @info "CRIAMOS A 2 MATCH - VITORIA len S = ", length(S)
            # @info "CRIAMOS A 2 MATCH - VITORIA S = ", S
            # con = @build_constraint(
            #     sum(lazy_model[:x][i, j] for (i, j) in S) + sum(lazy_model[:x][i, j] for (i, j) in T) <= length(h) + floor(length(T)/2),
            # )
            # MOI.submit(lazy_model, MOI.LazyConstraint(cb_data), con)
            
            # COMB
            S = []
            
            for i in 1:length(h)-1
                push!(S,(h[i], h[i+1]))
                # push!(nodeadded,h[i])
            end
            # push!(nodeadded,h[length(h)])
            push!(S,(h[1], h[length(h)]))
            @info "generating combs... S = ", S
            T = []
            @info "generating combs... t = ", t
            contador = 0
            HT = []
            AT = 0
            # ajuda = (length(t) & 1) == 1 ? 0 : 1
            nodeadded = []
            for i in 1:length(t)
                @info "generating combs... t[i] = ", t[i]
                @info "generating combs... t[i][1] = ", t[i][1]
                @info "generating combs... t[i][2] = ", t[i][2]
                @info "generating combs... t[i][3] = ", t[i][3]
                if length(t[i][2]) > 0
                    contador+=1
                    for aqui in t[i][2]
                        AT+=1
                        push!(HT,(t[i][1], aqui))
                        # push!(nodeadded,aqui)
                    end
                    if length(t[i][3]) > 0
                        for semcriatividade in t[i][3]
                            if length(semcriatividade) > 0
                                for aindasem in semcriatividade
                                    @info "CRIAMOS A COMB - NOT YET len nodeadded = ", length(nodeadded)
                                    @info "CRIAMOS A COMB - NOT YET nodeadded = ", nodeadded
                                    if !(any(x->x==aindasem[2], nodeadded))
                                        AT += 1
                                        push!(T,aindasem)
                                        push!(nodeadded,aindasem[2])
                                    end
                                end
                            end
                        end
                    end
                end
            end

            if contador > 2 && (contador & 1) == 1
                @info "CRIAMOS A COMB - VITORIA len S = ", length(S)
                @info "CRIAMOS A COMB - VITORIA S = ", S
                @info "CRIAMOS A COMB - VITORIA len HT = ", length(HT)
                @info "CRIAMOS A COMB - VITORIA HT = ", HT
                @info "CRIAMOS A COMB - VITORIA len T = ", length(T)
                @info "CRIAMOS A COMB - VITORIA T = ", T
                @info "CRIAMOS A COMB - VITORIA AT = ", AT
                @info "CRIAMOS A COMB - VITORIA ceil = ", ceil(3*contador/2)
                @info "CRIAMOS A COMB - VITORIA tudo = ", (length(h) + AT - ceil(3*contador/2))
                con = @build_constraint(
                    sum(lazy_model[:x][i, j] for (i, j) in S) + sum(lazy_model[:x][i, j] for (i, j) in T) <= length(h) + AT - ceil(3*contador/2),
                )
                MOI.submit(lazy_model, MOI.LazyConstraint(cb_data), con)
            end
        end


        if !(1 < length(cycle) < n)
            @debug "subtour_elimination_callback...done: only add a constraint if there is a cycle"
            return  # Only add a constraint if there is a cycle
        end

        # @debug "Found cycle of length $(length(cycle))"
        # S = []
        # for i in 1:length(cycle)-1
        #     push!(S, (cycle[i], cycle[i+1]))
        # end
        # push!(S, (cycle[length(cycle)], cycle[1]))
        S = [(i, j) for (i, j) in Iterators.product(cycle, cycle) if i < j]
        # @info "1 - cycle = ", cycle

        @debug "1 - S = ", S
        con = @build_constraint(
            sum(lazy_model[:x][i, j] for (i, j) in S) <= length(cycle) - 1,
        )
        MOI.submit(lazy_model, MOI.LazyConstraint(cb_data), con)
        return
    end

    MOI.set(lazy_model, MOI.LazyConstraintCallback(), subtour_elimination_callback)

    @debug "optimizing..."
    optimize!(lazy_model)
    @debug "optimizing...done"

    is_optimal = false
    stop_condition = ""
    status = MOI.get(lazy_model, MOI.TerminationStatus())
    if status == MOI.OPTIMAL
        # Uhuu! we solved the problem!
        @debug "Problem solved: MOI.OPTIMAL"
        is_optimal = true
        stop_condition = "Problem solved!"
        show_info("Problem solved!")
        show_info(string("Objective value found:",objective_value(lazy_model)))
        plot(exp_id,lazy_model[:x])
    else
        @debug "Problem remains: $status"
        stop_condition = "Aborted!"
    end

    min_gap = relative_gap(lazy_model)

    experiment_finish_time = get_time_in_ms()

    results = Dict(
        "exp_id" => exp_id,
        "method" => "lazy",
        "status" => status,
        "objective_value" => objective_value(lazy_model),
        "min_gap" => min_gap,
        "current_lower_bound" => "-",
        "best_lower_bound" => "-",
        "current_upper_bound" => "-",
        "best_upper_bound" => "-",
        "min_upper_bound" => "-",
        "max_lower_bound" => "-",
        "optimality_gap" => "-",
        "iterations_ran" => total_iterations,
        "is_optimal" => is_optimal,
        "stop_condition" => stop_condition,
        "start" => experiment_start_time,
        "finish" => experiment_finish_time,
        "duration" => experiment_finish_time - experiment_start_time,
    )

    for (key, value) in results
        show_result(exp_id, "lazy", string(key," "), string(value))
    end

    show_info("Solution Summary Report:")
    show_info(string(solution_summary(lazy_model,verbose=true)))

    save_step(exp_id,"lazy","finish","algorithm")

    @debug "status = ", status

    return results
    
end
