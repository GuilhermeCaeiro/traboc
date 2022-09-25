using Logging
using Compose
import Cairo, Fontconfig
include("christofides.jl")
include("onetree.jl")
include("utils.jl")
include("commons.jl")

#println(dirname(pwd()))

# Function to check if a graph consitutes a hamiltonian cycle.
# It assumes that the graph provided is symatric and connected.
# 
# Parameters
# graph -> a symetric SimpleGraph object
# Returns
# true if a hamiltonian cycle, or false otherwise
function check_hamiltonian_cycle(graph)
    nodes_checked = 0
    number_of_nodes = nv(graph)
    number_of_edges = ne(graph)
    previous_node = 0
    starting_node = 1
    next_node = starting_node

    # if there are only 2 nodes, there is no hamiltonian cycle
    if number_of_nodes <= 2
        #println("Two or less nodes.")
        return false
    end

    # attempts to follow a cycle where nodes_checked shout be equal number_of_nodes
    while true
        node = next_node
        node_degree = degree(graph, node)
        # all nodes must have degree == 2
        if node_degree != 2
            #println("Degree other than 2: node $node, degree $node_degree")
            return false
        end

        neighbor_nodes = neighbors(graph, node)
        next_node = deleteat!(neighbor_nodes, findall(x->x==previous_node, neighbor_nodes))[1] # probably inneficient, but there are only two values

        nodes_checked += 1
        previous_node = node

        if next_node == starting_node
            break
        end
    end

    #println(number_of_nodes, " ", nodes_checked)

    if number_of_nodes == nodes_checked
        return true
    else
        return false
    end
end


function factor_2_approximation(cost_matrix, n)
    graph = create_complete_graph(n)
    mst, c = minimum_spanning_tree(graph, cost_matrix, n) 
    # Converting the minimum spanning tree to a digraph.
    # The duplicated edges are generated during the conversion.
    dgraph = SimpleDiGraph(mst)
    ham_cycle_from_euler = euler_path(dgraph)
    return ham_cycle_from_euler
end



function lagrangean_relaxation(exp_id::String, testdatafile::String, max_iterations::Int64, epsilon::Float64)
    save_step(exp_id,"lagrangean_relaxation","START")

    upper_bounds = Array{Float64}(undef, 0)
    lower_bounds = Array{Float64}(undef, 0)
    epsilons = Array{Float64}(undef, 0)
    gaps = Array{Float64}(undef, 0)

    current_upper_bound = Inf
    current_lower_bound = -Inf
    best_upper_bound = Inf
    best_lower_bound = -Inf
    best_ub_sol = undef
    best_lb_sol = undef
    last_upper_bound_update = 0
    last_lower_boud_update = 0

    optimality_gap = Inf
    gap_threshold = 0.00001

    testdata = String(read(testdatafile));
    xml_graph = parsexml(testdata)
    original_cost_matrix, n = graph_to_cost_matrix(xml_graph)
    #original_cost_matrix = [
    #    0  30 26 50 40 
    #    30  0 24 40 50
    #    26 24  0 24 26
    #    50 40 24  0 30
    #    40 50 26 30  0
    #]
    #n = 5

    #current_cost_matrix = deepcopy(original_cost_matrix)
    u = zeros(1, n)

    # getting upper bound, christofides produces feasible solutions
    #christofides_sol = unite_and_hamilton(original_cost_matrix, n)
    #current_upper_bound = calculate_graph_cost(christofides_sol, original_cost_matrix, n) 

    total_iterations = 0

    save_step(exp_id,"lagrangean_relaxation","ITERATIONS")

    for iteration = 1:max_iterations
        total_iterations = iteration
        # updating cost matrix based on lagrangian multipliers
        current_cost_matrix = zeros(n, n)

        for i in 1:n
            for j in 1:n
                if i == j
                    continue
                end

                current_cost_matrix[i, j] = original_cost_matrix[i, j] - u[1, i] - u[1, j]
            end
        end

        @debug "current_cost_matrix " current_cost_matrix

        #ub_algorithm = "christofides"
        ub_algorithm = "factor2approximation"
        ub_solution = undef
        current_upper_bound = Inf

        # getting upper bound
        if ub_algorithm == "christofides"
            # christofides produces feasible solutions
            ub_solution = unite_and_hamilton(current_cost_matrix, n) 
        elseif ub_algorithm == "factor2approximation"
            # factor 2 approximation produces feasible solutions
            ub_solution = factor_2_approximation(current_cost_matrix, n) 
        else
            println("Invalid upper bound algorithm \"$ub_algorithm\".")
        end

        current_upper_bound = calculate_graph_cost(ub_solution, original_cost_matrix, n)

        
        

        #draw(PDF(string("christofides_", i, ".pdf"), 16cm, 16cm), gplot(christofides_sol))

        # getting lower bound
        one_tree_sol = one_tree_graph(current_cost_matrix, n)
        current_lower_bound = calculate_graph_cost(one_tree_sol, current_cost_matrix, n)
        #current_lower_bound = calculate_graph_cost(one_tree_sol, original_cost_matrix, n) 
        current_lower_bound = current_lower_bound + 2 * (sum(u))# ^ 2)

        if current_upper_bound < best_upper_bound
            best_upper_bound = current_upper_bound
            best_ub_sol = ub_solution
            last_upper_bound_update = iteration
        end

        if current_lower_bound > best_lower_bound
            best_lower_bound = current_lower_bound
            best_lb_sol = one_tree_sol
            last_lower_boud_update = iteration 
        end

        #println(u)

        #optimality_gap = (current_upper_bound - current_lower_bound) / current_upper_bound
        optimality_gap = (best_upper_bound - best_lower_bound) / best_upper_bound
        
        # lagrangian vars update
        edges_per_vertex = count_vertex_edges(one_tree_sol, n)
        G = 2 .- edges_per_vertex
        denominator = sum(G .^ 2)
        #mi = epsilon * (current_upper_bound - current_lower_bound) / denominator
        #mi = epsilon * (1.05 * current_upper_bound - current_lower_bound) / denominator
        #mi = epsilon * (best_upper_bound - current_lower_bound) / denominator
        mi = epsilon * (1.01 * best_upper_bound - current_lower_bound) / denominator
        u = u + mi * G

        #draw(PDF(string("christofides_", i, ".pdf"), 16cm, 16cm), gplot(christofides_sol, nodelabel=1:nv(one_tree_sol)))
        #draw(PDF(string("onetree_", i, ".pdf"), 16cm, 16cm), gplot(one_tree_sol, nodelabel=1:nv(one_tree_sol)))

        push!(gaps, optimality_gap)
        push!(upper_bounds, current_upper_bound)
        push!(lower_bounds, current_lower_bound)
        push!(epsilons, epsilon)


        if iteration - last_lower_boud_update > 20
            epsilon = epsilon / 2
        end

        if mod(iteration, 100) == 0
            print("Iteration ", iteration, " ")
            print(" ub ", current_upper_bound, " ")
            print(" lb ", current_lower_bound, " ")
            print(" gap ", optimality_gap, " ")
            println(" best boundaries (u/l) ", best_upper_bound, "/", best_lower_bound, " ")
            println("edges_per_vertex ", edges_per_vertex)
            println(" onetree cost on original costs ", calculate_graph_cost(one_tree_sol, original_cost_matrix, n), " ")
            println(" ub cost on original costs ", calculate_graph_cost(ub_solution, original_cost_matrix, n), " ")
            println(" onetree cost on current costs ", calculate_graph_cost(one_tree_sol, current_cost_matrix, n), " ")
            println(" ub cost on current costs ", calculate_graph_cost(ub_solution, current_cost_matrix, n), " ")
        end

        if optimality_gap == 0
            println("Optimal solution found. GAP ", optimality_gap, " BLB ", best_lower_bound, " BUB ", best_upper_bound, " CLB ", current_lower_bound, " CUB ", current_upper_bound, " UB ", calculate_graph_cost(christofides_sol, original_cost_matrix, n))
            println("min ub ", minimum(upper_bounds))
            println("max lb ", maximum(lower_bounds))
            println("min gap ", minimum(gaps))
            #return optimality_gap, best_lower_bound, best_upper_bound, best_ub_sol, best_lb_sol
            break
        end

        if optimality_gap <= gap_threshold
            println("Acceptable solution found. GAP ", optimality_gap, " BLB ", best_lower_bound, " BUB ", best_upper_bound, " CLB ", current_lower_bound, " CUB ", current_upper_bound, " UB ", calculate_graph_cost(christofides_sol, original_cost_matrix, n))
            println("min ub ", minimum(upper_bounds))
            println("max lb ", maximum(lower_bounds))
            println("min gap ", minimum(gaps))
            #return optimality_gap, best_lower_bound, best_upper_bound, best_ub_sol, best_lb_sol
            break
        end   

        if denominator == 0
            println("iteration ", iteration, " denominator 0 ", optimality_gap, " BLB ", best_lower_bound, " BUB ", best_upper_bound, " CLB ", current_lower_bound, " CUB ", current_upper_bound, " UB ")
            feasible = check_hamiltonian_cycle(one_tree_sol)

            # if the solution is feasible and the denominator is zero, that solution is optimal
            if feasible
                println("Optimal solution")
            else
                println("Stopping. Yet another lower bound found, but the denominator is zero and the solution is unfeasible.")
            end

            break
        end
        #break
        save_step(exp_id,"lagrangean_relaxation","IT_$iteration")

    end

    #println("ubs ", upper_bounds)
    #println("lbs ", lower_bounds)
    #println("gaps ", gaps)
    println("min ub ", minimum(upper_bounds))
    println("max lb ", maximum(lower_bounds))
    println("gap ", optimality_gap)
    println("iterations run ", total_iterations)

    save_step(exp_id,"lagrangean_relaxation","FINISH")
end

#lagrangean_relaxation("1", "", 1000, 1.0)
