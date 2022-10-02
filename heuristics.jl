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
        println(candidates)
        for candidate in candidates                
            if !(candidate in greedy_path)
                push!(greedy_path, candidate)
                break
            end
        end
    end

    println("Node order ", greedy_path)
    
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

#instance = [
#    0 5 1 2 
#    5 0 4 1
#    1 4 0 7
#    2 1 7 0
#]

#n = 4

#solution = nearest_neighbor(instance, n)
#println(calculate_graph_cost(solution, instance, n))

#gplot(solution, nodelabel=1:n)

function farthest_insertion(cost_matrix, n)
    ordered_tour = []
    remaining_nodes = collect(1:n)
    println("rn ", remaining_nodes)

    # retrieving initial tour 
    initial_edge = argmax(cost_matrix)
    push!(ordered_tour, initial_edge[1])
    deleteat!(remaining_nodes, findall(x->x == initial_edge[1], remaining_nodes))# pop!(remaining_nodes, initial_edge[1])
    push!(ordered_tour, initial_edge[2])
    deleteat!(remaining_nodes, findall(x->x == initial_edge[2], remaining_nodes))

    println("ot ", ordered_tour)
    println("rn ", remaining_nodes)

    while length(remaining_nodes) > 0
        farthest_node = 0
        farthest_distance = 0

        println("wrn ", remaining_nodes)

        # retreives the farthest node
        for remaining_node in remaining_nodes
            closest_distance = Inf
            for tour_node in ordered_tour
                println("rct ", remaining_node, " ", closest_distance, " ", tour_node)
                distance = cost_matrix[remaining_node, tour_node]
                println(distance)
                if distance < closest_distance
                    closest_distance = distance
                end
            end

            if closest_distance >= farthest_distance # >= to replace the node 0, that doesn't exist
                farthest_node = remaining_node
                farthest_distance = closest_distance
            end

            println("ft ", farthest_node, " ", farthest_distance)
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

            println("$farthest_node => $farthest_node, $node_i -> $distance_i; $farthest_node, $node_j -> $distance_j")

            distance = distance_i + distance_j - cost_matrix[node_i, node_j]

            println("subd ", distance)

            if distance <= closest_edge_distance
                closest_edge_index = idx
                closest_edge_distance = distance
            end
        end

        println("cei ", closest_edge_index)

        # adds the selected node to the tour
        # for some reason, julia doesnt raise an error 
        # if ordered_tour[(closest_edge_index + 1):end] ends up
        # being out of bounds. For that reason, no special treatment
        # for that case was given.
        ordered_tour = vcat(ordered_tour[1:closest_edge_index], [farthest_node], ordered_tour[(closest_edge_index + 1):end])

        println("nortour ", ordered_tour)

        # removes the selected node from remaining_nodes
        deleteat!(remaining_nodes, findall(x->x == farthest_node, remaining_nodes))
    end

    # builds the proper graph to be returned
    graph = SimpleGraph(n)

    for i in 1:n
        if i != n
            add_edge!(graph, ordered_tour[i], ordered_tour[i + 1])
        else
            add_edge!(graph, ordered_tour[i], ordered_tour[1])
        end
    end

    return graph
end

function 2opt()

end

instance = [
    0 5 1 2 
    5 0 4 1
    1 4 0 7
    2 1 7 0
]
n = 4
solution = nearest_neighbor(instance, n)
println(calculate_graph_cost(solution, instance, n))
#gplot(solution, nodelabel=1:n)

solution = farthest_insertion(instance, n)
println(calculate_graph_cost(solution, instance, n))
gplot(solution, nodelabel=1:n)