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

    return shortest_subtour, []

end

function selected_edges(x::Matrix{Float64}, n)
    return Tuple{Int,Int}[(i, j) for i in 1:n, j in 1:n if x[i, j] > 0.5]
end

subtour(x::Matrix{Float64}) = subtour(selected_edges(x, size(x, 1)), size(x, 1))
subtour(x::AbstractMatrix{VariableRef}) = subtour(value.(x))

function plot(exp_id::String, edges::Vector{Tuple{Int,Int}}, n)

    new_graph = SimpleGraph(n)

	for (i, j) in edges
        add_edge!(new_graph, i, j)
	end

    draw(PNG(string("./work/",exp_id,"/new_graph_", get_next_id(),".png"), 25cm, 25cm), gplot(new_graph, nodelabel=1:nv(new_graph)))

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

    # lower, upper = lagrangean_relaxation(900, false, file)
    # println("lower = ", lower)
    # println("upper = ", upper)
    lazy_model = build_tsp_model(optimizer, cost, n, lower, upper)
    function subtour_elimination_callback(cb_data)

        @debug "subtour_elimination_callback..."

        status = callback_node_status(cb_data, lazy_model)
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            @debug "subtour_elimination_callback...done: only run at integer solution"
            return  # Only run at integer solutions
        end
        cycle, c = subtour(callback_value.(cb_data, lazy_model[:x]))

        @debug begin
            plot(exp_id,callback_value.(cb_data, lazy_model[:x]))            
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
        # return
        if length(c) == 0
            @debug "subtour_elimination_callback...done: length(c) == 0"
            return
        end

        if exp_params["rac_strategy"] in "both|combs"
            @debug "generating combs..."
            @debug "1 - c = ", c
            T = [pop!(c)]
            H = [pop!(c)]
            push!(T, pop!(c))
            t = [(T[length(T)-1], T[length(T)])]
            ht = [(T[length(T)], H[length(H)])]
            push!(H, pop!(c))
            h = [(H[length(H)-1], H[length(H)])]
            push!(H, pop!(c))
            push!(h, (H[length(H)-1], H[length(H)]))
            push!(T, pop!(c))
            push!(ht, (T[length(T)], H[length(H)]))
            push!(H, pop!(c))
            push!(ht, (T[length(T)], H[length(H)]))
            push!(H, pop!(c))
            push!(h, (H[length(H)-1], H[length(H)]))
            push!(H, pop!(c))
            push!(h, (H[length(H)-1], H[length(H)]))
            push!(T, pop!(c))
            push!(ht, (T[length(T)], H[length(H)]))
            push!(H, pop!(c))
            push!(ht, (T[length(T)], H[length(H)]))
            # while length(c) > 0
            #     push!(H, pop!(c))
            #     push!(h, (H[length(H)-1], H[length(H)]))
            # end
            # push!(h, (H[length(H)], H[1]))
            @debug "1 - h = ", h
            @debug "1 - ht = ", ht
            @debug "generating combs...done"
            con = @build_constraint(
                sum(lazy_model[:x][i, j] for (i, j) in h) + sum(lazy_model[:x][i, j] for (i, j) in t) <= length(H) + length(T) - ceil(3*length(T)/2),
            )
            # con = @build_constraint(
            #     sum(lazy_model[:x][i, j] for (i, j) in h) + sum(lazy_model[:x][i, j] for (i, j) in ht) <= length(H) + floor(length(T)/2),
            # )
            MOI.submit(lazy_model, MOI.LazyConstraint(cb_data), con)
            return
        end

    end

    MOI.set(lazy_model, MOI.LazyConstraintCallback(), subtour_elimination_callback)

    @debug "optimizing..."
    optimize!(lazy_model)
    @debug "optimizing...done"

    status = MOI.get(lazy_model, MOI.TerminationStatus())
    if status == MOI.OPTIMAL
        # Uhuu! we solved the problem!
        @debug "Problem solved: MOI.OPTIMAL"
        show_info("Problem solved!")
        show_info(string("Objective value found:",objective_value(lazy_model)))
        plot(exp_id,lazy_model[:x])
    else
        @debug "Problem remains: $status"
    end

    show_result(exp_id, "lazy", "status ", status)
    show_result(exp_id, "lazy", "objective_value ", objective_value(lazy_model))

    save_step(exp_id,"lazy","finish","algorithm")

    @debug "status = ", status
    
end