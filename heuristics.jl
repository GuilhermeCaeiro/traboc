using Graphs
using GraphPlot
include("utils.jl")
# Assumes that the graph is complete and that vertex i 
# isn't reacheable from i itself
# Parameter
# graph -> 
function nearest_neighbor(cost_matrix, n)
    greedy_path = [1]
    
    for i in 1:(n-1)
        candidates = sortperm(cost_matrix[greedy_path[end], :])
        for candidate in candidates                
            if !(candidate in greedy_path)
                push!(greedy_path, candidate)
                break
            end
        end
    end
    
    # building graph
    greedy_hamiltonian_cycle = SimpleGraph(n)

    for i in 1:n
        if i != n
            add_edge!(greedy_hamiltonian_cycle, greedy_path[i], greedy_path[i + 1])
        else
            add_edge!(greedy_hamiltonian_cycle, greedy_path[i], 1)
        end
    end

    return greedy_hamiltonian_cycle
end

function build_graph_from_tour(tour)
    # builds the proper graph to be returned
    graph = SimpleGraph(length(tour))

    for i in 1:n
        if i != n
            add_edge!(graph, tour[i], tour[i + 1])
        else
            add_edge!(graph, tour[i], tour[1])
        end
    end

    return graph
end

function farthest_insertion(cost_matrix, n)
    ordered_tour = []
    remaining_nodes = collect(1:n)

    # retrieving initial tour 
    initial_edge = argmax(cost_matrix)
    push!(ordered_tour, initial_edge[1])
    deleteat!(remaining_nodes, findall(x->x == initial_edge[1], remaining_nodes))# pop!(remaining_nodes, initial_edge[1])
    push!(ordered_tour, initial_edge[2])
    deleteat!(remaining_nodes, findall(x->x == initial_edge[2], remaining_nodes))

    while length(remaining_nodes) > 0
        farthest_node = 0
        farthest_distance = 0

        # retreives the farthest node
        for remaining_node in remaining_nodes
            closest_distance = Inf
            for tour_node in ordered_tour
                distance = cost_matrix[remaining_node, tour_node]
                if distance < closest_distance
                    closest_distance = distance
                end
            end

            if closest_distance >= farthest_distance # >= to replace the node 0, that doesn't exist
                farthest_node = remaining_node
                farthest_distance = closest_distance
            end
        end

        # finds where to insert the farthest node
        closest_edge_index = 0
        closest_edge_distance = Inf

        for idx in eachindex(ordered_tour)
            distance_i = 0
            distance_j = 0
            node_i = ordered_tour[idx]
            node_j = 0

            if idx == length(ordered_tour)
                node_j = ordered_tour[1]
            else
                node_j = ordered_tour[idx + 1]
            end

            distance_i = cost_matrix[farthest_node, node_i]
            distance_j = cost_matrix[farthest_node, node_j]

            distance = distance_i + distance_j - cost_matrix[node_i, node_j]

            if distance <= closest_edge_distance
                closest_edge_index = idx
                closest_edge_distance = distance
            end
        end

        # adds the selected node to the tour
        # for some reason, julia doesnt raise an error 
        # if ordered_tour[(closest_edge_index + 1):end] ends up
        # being out of bounds. For that reason, no special treatment
        # for that case was given.
        ordered_tour = vcat(ordered_tour[1:closest_edge_index], [farthest_node], ordered_tour[(closest_edge_index + 1):end])

        # removes the selected node from remaining_nodes
        deleteat!(remaining_nodes, findall(x->x == farthest_node, remaining_nodes))
    end

    graph = build_graph_from_tour(ordered_tour)

    return graph
end

function get_route_from_graph(graph, n)
    route = []
    available_nodes = collect(1:n)
    current_node = 1

    while (length(available_nodes) > 0)
        push!(route, current_node)
        deleteat!(available_nodes, findall(x->x == current_node, available_nodes))

        for neighbor in neighbors(graph, current_node)
            if !(neighbor in route)
                current_node = neighbor
                break
            end
        end
    end

    return route    
end

function two_opt_cost_variation(route, cost_matrix, i, j)
    cost_new_edge_1 = cost_matrix[route[i - 1], route[j - 1]]
    cost_new_edge_2 = cost_matrix[route[i], route[j]]
    cost_old_edge_1 = cost_matrix[route[i - 1], route[i]]
    cost_old_edge_2 = cost_matrix[route[j - 1], route[j]]
    return cost_new_edge_1 + cost_new_edge_2 - cost_old_edge_1 - cost_old_edge_2
end

# [1 3] -> 1
# [2 4] -> 1
# [1 2] -> 5
# [3 4] -> 7 

function two_opt_iteration(route, cost, cost_matrix)
    for i in 2:(n-1)
        for j in i+1:n
            if j - i == 1
                continue
            end

            cost_variation = two_opt_cost_variation(route, cost_matrix, i, j)

            if cost_variation < 0
                route[i:j-1] = reverse(route[i:j-1])
                cost = cost + cost_variation
                return route, cost
            end

        end
    end
    return route, cost
end

# based on the available solution at
# https://stackoverflow.com/questions/53275314/2-opt-algorithm-to-solve-the-travelling-salesman-problem-in-python
# from user Fradge
function two_opt(graph, cost_matrix, n; max_iterations = 1000, cost = Inf)
    # gets the cost of the graph, if not provided
    if cost == Inf
        cost = calculate_graph_cost(graph, cost_matrix, n)
    end

    if n <= 3
        return graph, cost
    end

    # retrieves the route/path (non cyclical)
    route = get_route_from_graph(graph, n)

    for iteration in 1:max_iterations
        route, current_cost = two_opt_iteration(route, cost, cost_matrix)
        
        if current_cost < cost
            cost = current_cost
        else
            break
        end
    end
    graph = build_graph_from_tour(route)

    return graph, cost
end
