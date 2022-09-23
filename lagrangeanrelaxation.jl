using Compose
import Cairo, Fontconfig
include("cristofides.jl")
include("onetree.jl")
println(dirname(pwd()))

function lagrangean_relaxation()

    max_iterations = 10000
    upper_bounds = Array{Float64}(undef, 0)
    lower_bounds = Array{Float64}(undef, 0)
    gaps = Array{Float64}(undef, 0)

    current_upper_bound = Inf
    current_lower_bound = -Inf
    best_upper_bound = Inf
    best_lower_bound = -Inf
    best_ub_sol = undef
    best_lb_sol = undef

    optimality_gap = Inf
    gap_threshold = 0.00001
    epsilon = 1


    testdata = String(read("/home/guilherme/Documentos/workspace/traboc/test_cristofides.xml"));
    xml_graph = parsexml(testdata)
    #original_cost_matrix, n = graph_to_cost_matrix(xml_graph)
    original_cost_matrix = [
        0  30 26 50 40 
        30  0 24 40 50
        26 24  0 24 26
        50 40 24  0 30
        40 50 26 30  0
    ]
    n = 5

    #current_cost_matrix = deepcopy(original_cost_matrix)
    u = zeros(1, n)

    # getting upper bound, christofides produces feasible solutions
    #christofides_sol = unite_and_hamilton(original_cost_matrix, n)
    #current_upper_bound = calculate_graph_cost(christofides_sol, original_cost_matrix, n) 

    total_iterations = 0

    for i = 1:max_iterations
        # updating cost matrix based on lagrangian multipliers
        total_iterations = i
        current_cost_matrix = zeros(n, n)

        for i in 1:n
            for j in 1:n
                if i == j
                    continue
                end

                current_cost_matrix[i, j] = original_cost_matrix[i, j] - u[1, i] - u[1, j]
            end
        end

        #println(current_cost_matrix)

        #println("current_cost_matrix ", current_cost_matrix)

        # getting upper bound, christofides produces feasible solutions
        christofides_sol = unite_and_hamilton(current_cost_matrix, n)
        #christofides_sol = unite_and_hamilton(original_cost_matrix, n)
        #current_upper_bound = calculate_graph_cost(christofides_sol, current_cost_matrix, n) 
        current_upper_bound = calculate_graph_cost(christofides_sol, original_cost_matrix, n) 
        

        #draw(PDF(string("christofides_", i, ".pdf"), 16cm, 16cm), gplot(christofides_sol))

        # getting lower bound
        one_tree_sol = one_tree_graph(current_cost_matrix, n)
        current_lower_bound = calculate_graph_cost(one_tree_sol, current_cost_matrix, n)
        #current_lower_bound = calculate_graph_cost(one_tree_sol, original_cost_matrix, n) 
        current_lower_bound = current_lower_bound + 2 * (sum(u))# ^ 2)

        if current_upper_bound < best_upper_bound
            best_upper_bound = current_upper_bound
            best_ub_sol = christofides_sol
        end

        if current_lower_bound > best_lower_bound
            best_lower_bound = current_lower_bound
            best_lb_sol = one_tree_sol
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
        mi = epsilon * (1.05 * best_upper_bound - current_lower_bound) / denominator
        u = u + mi * G

        #draw(PDF(string("christofides_", i, ".pdf"), 16cm, 16cm), gplot(christofides_sol, nodelabel=1:nv(one_tree_sol)))
        #draw(PDF(string("onetree_", i, ".pdf"), 16cm, 16cm), gplot(one_tree_sol, nodelabel=1:nv(one_tree_sol)))

        push!(gaps, optimality_gap)
        push!(upper_bounds, current_upper_bound)
        push!(lower_bounds, current_lower_bound)

        if mod(i, 100) == 0
            print("Iteration ", i, " ")
            print(" ub ", current_upper_bound, " ")
            print(" lb ", current_lower_bound, " ")
            print(" gap ", optimality_gap, " ")
            println(" best boundaries (u/l) ", best_upper_bound, "/", best_lower_bound, " ")
            println("edges_per_vertex ", edges_per_vertex)
            println(" onetree cost on original costs ", calculate_graph_cost(one_tree_sol, original_cost_matrix, n), " ")
            println(" christofides cost on original costs ", calculate_graph_cost(christofides_sol, original_cost_matrix, n), " ")
            println(" onetree cost on current costs ", calculate_graph_cost(one_tree_sol, current_cost_matrix, n), " ")
            println(" christofides cost on current costs ", calculate_graph_cost(christofides_sol, current_cost_matrix, n), " ")
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
            println("iteration ", i, " denominator 0 ", optimality_gap, " BLB ", best_lower_bound, " BUB ", best_upper_bound, " CLB ", current_lower_bound, " CUB ", current_upper_bound, " UB ")
            # TODO implementar verificacao de viabilidade da solucao
            feasible = false

            # if the solution is feasible and the denominator is zero, that solution is optimal
            if feasible
                println("Optimal solution")
            else
                println("Stopping. Yet another lower bound found, but the denominator is zero and the solution is unfeasible.")
            end

            break
        end
        #break
    end

    #println("ubs ", upper_bounds)
    #println("lbs ", lower_bounds)
    #println("gaps ", gaps)
    println("min ub ", minimum(upper_bounds))
    println("max lb ", maximum(lower_bounds))
    println("gap ", optimality_gap)
    println("iterations run ", total_iterations)
end

lagrangean_relaxation()