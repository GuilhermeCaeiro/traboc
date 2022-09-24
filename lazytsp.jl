using JuMP, Gurobi
include("utils.jl")

function build_tsp_model(d, n, lower = -Inf, upper = +Inf)
    model = Model(Gurobi.Optimizer)
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
    return model
end

function subtour(edges::Vector{Tuple{Int,Int}}, n)
    shortest_subtour, unvisited = collect(1:n), Set(collect(1:n))
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
        end
    end
    return shortest_subtour
end

function selected_edges(x::Matrix{Float64}, n)
    return Tuple{Int,Int}[(i, j) for i in 1:n, j in 1:n if x[i, j] > 0.5]
end

subtour(x::Matrix{Float64}) = subtour(selected_edges(x, size(x, 1)), size(x, 1))
subtour(x::AbstractMatrix{VariableRef}) = subtour(value.(x))

function lazy(file, lower = -Inf, upper = Inf)
    testdata = String(read(file));
    xml_graph = parsexml(testdata)
    cost, n = graph_to_cost_matrix(xml_graph)
    # lower, upper = lagrangean_relaxation(900, false, file)
    # println("lower = ", lower)
    # println("upper = ", upper)
    lazy_model = build_tsp_model(cost, n, lower, upper)
    function subtour_elimination_callback(cb_data)
        status = callback_node_status(cb_data, lazy_model)
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            return  # Only run at integer solutions
        end
        cycle = subtour(callback_value.(cb_data, lazy_model[:x]))
        if !(1 < length(cycle) < n)
            return  # Only add a constraint if there is a cycle
        end
        @debug "Found cycle of length $(length(cycle))"
        S = [(i, j) for (i, j) in Iterators.product(cycle, cycle) if i < j]
        con = @build_constraint(
            sum(lazy_model[:x][i, j] for (i, j) in S) <= length(cycle) - 1,
        )
        MOI.submit(lazy_model, MOI.LazyConstraint(cb_data), con)
        return
    end
    MOI.set(lazy_model, MOI.LazyConstraintCallback(), subtour_elimination_callback)
    optimize!(lazy_model)
    status = MOI.get(lazy_model, MOI.TerminationStatus())
    if status == MOI.OPTIMAL
        # Ok, we solved the problem!
        # println(objective_value(lazy_model))
        show_info(objective_value(lazy_model))
    end
    @debug "status = ", status
    
end