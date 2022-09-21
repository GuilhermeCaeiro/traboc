using Compose
import Cairo, Fontconfig
include("cristofides.jl")
include("onetree.jl")
println(dirname(pwd()))

function lagrangean_relaxation()

    max_iterations = 100
    upper_bounds = Array{Float64}(undef, 0)
    lower_bounds = Array{Float64}(undef, 0)
    gaps = Array{Float64}(undef, 0)

    current_upper_bound = Inf
    current_lower_bound = -Inf
    gap_threshold = 0.00001
    epsilon = 1


    testdata = String(read("/home/guilherme/Documentos/workspace/traboc/test_cristofides.xml"));
    xml_graph = parsexml(testdata)
    original_cost_matrix, n = graph_to_cost_matrix(xml_graph)

    #current_cost_matrix = deepcopy(original_cost_matrix)
    u = zeros(1, n)

    # getting upper bound, christofides produces feasible solutions
    #christofides_sol = unite_and_hamilton(original_cost_matrix, n)
    #current_upper_bound = calculate_graph_cost(christofides_sol, original_cost_matrix, n) 

    for i = 1:max_iterations
        # updating cost matrix based on lagrangian multipliers
        current_cost_matrix = zeros(n, n)
        for j in 1:n
            current_cost_matrix[:, j] = original_cost_matrix[:, j] .- u[1, j]
        end

        println("current_cost_matrix ", current_cost_matrix)

        # getting upper bound, christofides produces feasible solutions
        christofides_sol = unite_and_hamilton(current_cost_matrix, n)
        #current_upper_bound = calculate_graph_cost(christofides_sol, current_cost_matrix, n) 
        current_upper_bound = calculate_graph_cost(christofides_sol, original_cost_matrix, n) 

        draw(PDF(string("christofides_", i, ".pdf"), 16cm, 16cm), gplot(christofides_sol))

        # getting lower bound
        one_tree_sol = one_tree_graph(current_cost_matrix, n)
        current_lower_bound = calculate_graph_cost(one_tree_sol, current_cost_matrix, n)

        optimality_gap = (current_upper_bound - current_lower_bound) / current_upper_bound
        
        # lagrangia vars update
        edges_per_vertex = count_vertex_edges(one_tree_sol, n)
        denominador = sum((2 .- edges_per_vertex) .^ 2)
        if denominador != 0
            mi = epsilon * (1.06*current_upper_bound - current_lower_bound) / (denominador)
            u = u + mi * (2 .- edges_per_vertex)
        else
            println("DENOMINADOR IGUAL A 0, se a solução atual for viavel então é otima")
            draw(PDF(string("./output/one_tree_sol_SolDoDenominador0", ".pdf"), 16cm, 16cm), gplot(one_tree_sol, nodelabel=1:nv(one_tree_sol)))
            draw(PDF(string("./output/christofides_sol_SolDoDenominador0", ".pdf"), 16cm, 16cm), gplot(christofides_sol, nodelabel=1:nv(christofides_sol)))
            break
        end

        println("ub ", current_upper_bound)
        println("lb ", current_lower_bound)
        println("gap ", optimality_gap)
        println("edges_per_vertex ", edges_per_vertex)
        println("mi ", mi)
        println("u ", u)

        push!(gaps, optimality_gap)
        push!(upper_bounds, current_upper_bound)
        push!(lower_bounds, current_lower_bound)

        #break
    end

    println("ubs ", upper_bounds)
    println("lbs ", lower_bounds)
    println("gaps ", gaps)
    println("min ub ", minimum(upper_bounds))
    println("max lb ", maximum(lower_bounds))
    println("min gap ", minimum(gaps))
end

lagrangean_relaxation()
