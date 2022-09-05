using JuMP, Gurobi

function integer_problem(instance, n)
    # code from https://github.com/PiotrZakrzewski/julia-tsp/blob/main/TSP.
    # that source is related to https://www.youtube.com/watch?v=IYae_sfzKMs
    # and based on https://co-enzyme.fr/blog/traveling-salesman-problem-tsp-in-cplex-opl-with-miller-tucker-zemlin-mtz-formulation/

    model = Model(Gurobi.Optimizer)
    # route is an adjence matrix representing a route traveled
    route=@variable(model, route[1:n, 1:n], Bin)
    # mtzu is a helper variable to ensure no subtours are allowed (only one continous tour)
    # see MTZ constraint
    mtzu = @variable(model, mtzu[1:n], Int)

    # ensure all events are planned
    @constraint(model, [i = 1:n], sum(route[i, :]) == 1.0)
    # ensure there is just one route
    @constraint(model, [c = 1:n], sum(route[:, c]) == 1.0)
    # disallow traveling to itself
    @constraint(model, [j = 1:n], route[j, j] == 0)

    # MTZ constraints for removing subtours
    @constraint(model, [ui = 1:n, uj = 2:n], mtzu[ui] + route[ui, uj] <= mtzu[uj]+ (n - 1) * (1 - route[ui, uj]) )

    traveltime = instance .* route 
    @objective(model, Min, sum(traveltime))
    optimize!(model)
    @show termination_status(model)
    @show primal_status(model)
    @show dual_status(model)
    @show objective_value(model)
    @show value.(route)

    return objective_value(model)

end

function linear_relaxation(instance, n)
    model = Model(Gurobi.Optimizer)
    # route is an adjence matrix representing a route traveled
    route=@variable(model, 0 <= route[1:n, 1:n] <= 1)
    # mtzu is a helper variable to ensure no subtours are allowed (only one continous tour)
    # see MTZ constraint
    mtzu = @variable(model, mtzu[1:n])

    # ensure all events are planned
    @constraint(model, [i = 1:n], sum(route[i, :]) == 1.0)
    # ensure there is just one route
    @constraint(model, [c = 1:n], sum(route[:, c]) == 1.0)
    # disallow traveling to itself
    @constraint(model, [j = 1:n], route[j, j] == 0)

    # MTZ constraints for removing subtours
    @constraint(model, [ui = 1:n, uj = 2:n], mtzu[ui] + route[ui, uj] <= mtzu[uj]+ (n - 1) * (1 - route[ui, uj]) )

    traveltime = instance .* route 
    @objective(model, Min, sum(traveltime))
    optimize!(model)
    @show termination_status(model)
    @show primal_status(model)
    @show dual_status(model)
    @show objective_value(model)
    @show value.(route)

    return objective_value(model)
end

# discards the subtour removal constraints.
function ap_linear_relaxation(instance, n)
    model = Model(Gurobi.Optimizer)
    # route is an adjence matrix representing a route traveled
    route=@variable(model, 0 <= route[1:n, 1:n] <= 1)
    
    # ensure all events are planned
    @constraint(model, [i = 1:n], sum(route[i, :]) == 1.0)
    # ensure there is just one route
    @constraint(model, [c = 1:n], sum(route[:, c]) == 1.0)
    # disallow traveling to itself
    @constraint(model, [j = 1:n], route[j, j] == 0)

    traveltime = instance .* route 
    @objective(model, Min, sum(traveltime))
    optimize!(model)
    @show termination_status(model)
    @show primal_status(model)
    @show dual_status(model)
    @show objective_value(model)
    @show value.(route)

    return objective_value(model)
end

# discards the subtour removal constraints.
function ap_integer_relaxation(instance, n)
    model = Model(Gurobi.Optimizer)
    # route is an adjence matrix representing a route traveled
    route=@variable(model, route[1:n, 1:n], Bin)
    
    # ensure all events are planned
    @constraint(model, [i = 1:n], sum(route[i, :]) == 1.0)
    # ensure there is just one route
    @constraint(model, [c = 1:n], sum(route[:, c]) == 1.0)
    # disallow traveling to itself
    @constraint(model, [j = 1:n], route[j, j] == 0)

    traveltime = instance .* route 
    @objective(model, Min, sum(traveltime))
    optimize!(model)
    @show termination_status(model)
    @show primal_status(model)
    @show dual_status(model)
    @show objective_value(model)
    @show value.(route)

    return objective_value(model)
end

function check_cycle(edges, candidate_edge)
    println("rawedges ", edges)
    if length(edges) == 0
        println("variable edgs empty")
        return false
    end
    
    edges = reduce(hcat, edges)'
    println("check_cycle ", edges, " ", candidate_edge, " ", typeof(edges), " ", candidate_edge[1] in edges, " ", candidate_edge[2] in edges)

    if candidate_edge[1] in edges && candidate_edge[2] in edges
        println(edges, " ", candidate_edge[1], " ", candidate_edge[2], " ", true)
        return true
    else
        println(edges, " ", candidate_edge[1], " ", candidate_edge[2], " ", false)
        return false
    end
end

# based on http://web.tecnico.ulisboa.pt/~mcasquilho/CD_Casquilho/TSP_Raffensperger.ppsx
# implements the Prim's algorithm
function one_tree_construction(cost_edge, n, one_value=1, start_node=2)
    # sorting edges by the cost
    sorted_cost_edge = sortslices(cost_edge, dims=1)
    removed_edges = []
    spanning_tree = []
    covered_nodes = Set([start_node])

    println(sorted_cost_edge)

    # removing edges conected to the first node (1)
    for i in 1:size(sorted_cost_edge, 1)
        if one_value in sorted_cost_edge[i, 2:end]
            push!(removed_edges, sorted_cost_edge[i, 2:end])
            println(sorted_cost_edge[i])
            sorted_cost_edge[i, :] = [-1, -1, -1]
        end
    end

    # this while keeps trying to add new edges until all nodes are covered
    while length(covered_nodes) < (n - 1) # one node was removed earlier
        println("Starting while")
        for i in 1:size(sorted_cost_edge, 1)
            println(sorted_cost_edge[i, :])
            if -1 in sorted_cost_edge[i, :]
                continue
            end

            # passing a copy to avoid error complaining about "resizing shared data"
            if check_cycle(deepcopy(spanning_tree), sorted_cost_edge[i, 2:end]) 
                continue
            end

            if sorted_cost_edge[i, 2] in covered_nodes || sorted_cost_edge[i, 3] in covered_nodes
                push!(spanning_tree, sorted_cost_edge[i, 2:end])
                # don't know which node is new, so both are added to the set (sets don't repeat values)
                push!(covered_nodes, sorted_cost_edge[i, 2]) 
                push!(covered_nodes, sorted_cost_edge[i, 3])
                sorted_cost_edge[i, :] = [-1, -1, -1]
            end
        end
    end
    println(covered_nodes)
    println(spanning_tree)


    # adding the node removed node and its two cheapest edges
    println("removed_edges", removed_edges)
    push!(covered_nodes, one_value) 
    push!(spanning_tree, removed_edges[1])
    push!(spanning_tree, removed_edges[2])


    println(covered_nodes)
    println(spanning_tree)

    return spanning_tree

end

function one_tree_lagrangean_relaxation(instance, n, one_tree)

end

instance = [
    0   2   5   6  12  10
    2   0   3   4  12  12
    5   3   0   1   9  11
    6   4   1   0   8  10
    12 12   9   8   0   2
    10 12  11  10   2   0
]

n = 6

optimal_value = integer_problem(instance, n)
println(optimal_value)

linr_value = linear_relaxation(instance, n)
println(linr_value)

ap_linr_value = ap_linear_relaxation(instance, n)
println(ap_linr_value)

ap_int_value = ap_integer_relaxation(instance, n)
println(ap_int_value)


cost_edge = [
    2 1 2
    3 2 3
    5 2 4
    7 1 3
    1 3 4
    10 3 5
    8 4 5
    2 5 6
    10 1 6
]

one_tree_construction(cost_edge, n, 1, 2)