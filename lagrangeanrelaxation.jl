using Graphs
using EzXML
using GraphsMatching
using GraphPlot
using Logging
using Compose
import Cairo
import Fontconfig

include("utils.jl")
include("logs.jl")
include("commons.jl")

include("christofides.jl")
include("factor2approximation.jl")
include("onetree.jl")

#
# Function lagrangean_relaxation
# obtains the lagrangean_relaxation with christofides method 
#  
# Parameters
# exp_id -> experiment id
# testdatafile -> tesd data file
# max_iterations -> max iterations for lagrangean relaxation
# epsilon -> epsilon gap
# mi_option -> mi function option
# Returns 
# void
function lagrangean_relaxation(exp_id::String, testdatafile::String, ub_algorithm::String, max_iterations::Int64, gap_threshold::Float64, epsilon::Float64, mi_option::String, epsilon_strategy::String, check_point::Int64)

    save_step(exp_id,"lagrangean_relaxation","start","method")

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

    last_lower_boud_update = 0
    epsilon_min = 0.0001

    testdata = String(read(testdatafile));
    xml_graph = parsexml(testdata)
    
    save_step(exp_id,"lagrangean_relaxation:graph_to_cost_matrix","start","step")
    original_cost_matrix, n = graph_to_cost_matrix(xml_graph)
    save_step(exp_id,"lagrangean_relaxation:graph_to_cost_matrix","finish","step")

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

    save_step(exp_id,"lagrangean_relaxation:iterations","start","step")

    for iteration = 1:max_iterations
        save_step(exp_id,"lagrangean_relaxation:iterations","start","iteration_$iteration")

        # updating cost matrix based on lagrangian multipliers
        total_iterations = iteration
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
        #ub_algorithm = "factor2approximation"
        ub_solution = undef
        second_ub_solution = undef # only used if ub_algorithm is "christofidesandfactor2"
        current_upper_bound = Inf

        # getting upper bound
        if ub_algorithm == "christofides"
            # christofides produces feasible solutions
            ub_solution = unite_and_hamilton(exp_id, iteration, current_cost_matrix, n)
        elseif ub_algorithm == "factor2approximation"
            # factor 2 approximation produces feasible solutions
            ub_solution = factor_2_approximation(exp_id, iteration, current_cost_matrix, n) 
        elseif ub_algorithm == "christofidesandfactor2"
            # 
            ub_solution = unite_and_hamilton(exp_id, iteration, current_cost_matrix, n)
            second_ub_solution = factor_2_approximation(exp_id, iteration, current_cost_matrix, n) 
        else
            println("Invalid upper bound algorithm \"$ub_algorithm\".")
        end

        if ub_algorithm == "christofidesandfactor2"
            first_upper_bound = calculate_graph_cost(ub_solution, original_cost_matrix, n)
            second_upper_bound = calculate_graph_cost(second_ub_solution, original_cost_matrix, n)

            if second_upper_bound < first_upper_bound
                ub_solution = second_ub_solution
            end

            current_upper_bound = minimum([first_upper_bound, second_upper_bound])

        else
            current_upper_bound = calculate_graph_cost(ub_solution, original_cost_matrix, n)
        end

        # getting upper bound, christofides produces feasible solutions
        #save_step(exp_id,"lagrangean_relaxation:unite_and_hamilton","start","iteration_$iteration")
        #christofides_sol = unite_and_hamilton(exp_id, iteration, current_cost_matrix, n)
        #save_step(exp_id,"lagrangean_relaxation:unite_and_hamilton","finish","iteration_$iteration")

        #save_step(exp_id,"lagrangean_relaxation:calculate_graph_cost:upper_bound","start","iteration_$iteration")
        #current_upper_bound = calculate_graph_cost(christofides_sol, original_cost_matrix, n) 
        #save_step(exp_id,"lagrangean_relaxation:calculate_graph_cost:upper_bound","finish","iteration_$iteration")

        # getting lower bound
        save_step(exp_id,"lagrangean_relaxation:one_tree_graph","start","iteration_$iteration")
        one_tree_sol = one_tree_graph(exp_id, iteration, current_cost_matrix, n)
        save_step(exp_id,"lagrangean_relaxation:one_tree_graph","finish","iteration_$iteration")

        save_step(exp_id,"lagrangean_relaxation:calculate_graph_cost:lower_bound","start","iteration_$iteration")
        current_lower_bound = calculate_graph_cost(one_tree_sol, current_cost_matrix, n)
        save_step(exp_id,"lagrangean_relaxation:calculate_graph_cost:lower_bound","finish","iteration_$iteration")
 
        current_lower_bound = current_lower_bound + 2 * (sum(u))

        if current_upper_bound < best_upper_bound
            best_upper_bound = current_upper_bound
            best_ub_sol = ub_solution
        end

        if current_lower_bound > best_lower_bound
            best_lower_bound = current_lower_bound
            best_lb_sol = one_tree_sol
            last_lower_boud_update = iteration
        end

        save_result(exp_id, "lagrangean_relaxation:$iteration", "iteration", iteration)
        save_result(exp_id, "lagrangean_relaxation:$iteration", "upper_bound", best_upper_bound)
        save_result(exp_id, "lagrangean_relaxation:$iteration", "lower_bound", best_lower_bound)
        save_result(exp_id, "lagrangean_relaxation:$iteration", "upper_bound_sol", best_ub_sol)
        save_result(exp_id, "lagrangean_relaxation:$iteration", "lower_bound_sol", best_lb_sol)

        @debug u

        #optimality_gap = (current_upper_bound - current_lower_bound) / current_upper_bound
        optimality_gap = (best_upper_bound - best_lower_bound) / best_upper_bound
        
        save_result(exp_id, "lagrangean_relaxation:$iteration", "optimality_gap", optimality_gap)

        # lagrangian vars update
        edges_per_vertex = count_vertex_edges(one_tree_sol, n)
        save_result(exp_id, "lagrangean_relaxation:$iteration", "edges_per_vertex", edges_per_vertex)

        G = 2 .- edges_per_vertex
        denominator = sum(G .^ 2)

        #mi function selector
        if mi_option == "current" 
            mi = epsilon * (current_upper_bound - current_lower_bound) / denominator
        elseif mi_option == "best"
            mi = epsilon * (best_upper_bound - current_lower_bound) / denominator
        elseif mi_option == "5pct"
            mi = epsilon * (1.05 * best_upper_bound - current_lower_bound) / denominator
        elseif mi_option == "1pct"
            mi = epsilon * (1.01 * best_upper_bound - current_lower_bound) / denominator
        else
            @error "Mi function option não informada."
        end

        save_result(exp_id, "lagrangean_relaxation:$iteration", "mi", mi)

        u = u + mi * G

        #draw(PDF(string("christofides_", iteration, ".pdf"), 16cm, 16cm), gplot(christofides_sol, nodelabel=1:nv(one_tree_sol)))
        #draw(PDF(string("onetree_", iteration, ".pdf"), 16cm, 16cm), gplot(one_tree_sol, nodelabel=1:nv(one_tree_sol)))

        push!(gaps, optimality_gap)
        push!(upper_bounds, current_upper_bound)
        push!(lower_bounds, current_lower_bound)

        

        if mod(iteration, check_point) == 0            
            show_result(exp_id, "lagrangean_relaxation:$iteration:check_point", "Iteration ", iteration)
            show_result(exp_id, "lagrangean_relaxation:$iteration:check_point", "current_upper_bound", current_upper_bound)
            show_result(exp_id, "lagrangean_relaxation:$iteration:check_point", "current_lower_bound", current_lower_bound)
            show_result(exp_id, "lagrangean_relaxation:$iteration:check_point", "epsilon_strategy", epsilon_strategy)
            show_result(exp_id, "lagrangean_relaxation:$iteration:check_point", "epsilon", epsilon)
            show_result(exp_id, "lagrangean_relaxation:$iteration:check_point", "optimality_gap", optimality_gap)
            #show_result(exp_id, "lagrangean_relaxation:$iteration:check_point", "best boundaries (upper/lower) ", string(best_upper_bound) + " / " + string(best_lower_bound))
            show_result(exp_id, "lagrangean_relaxation:$iteration:check_point", "edges_per_vertex ", edges_per_vertex)
            show_result(exp_id, "lagrangean_relaxation:$iteration:check_point", "onetree cost on original costs ", calculate_graph_cost(one_tree_sol, original_cost_matrix, n))
            show_result(exp_id, "lagrangean_relaxation:$iteration:check_point", "christofides cost on original costs ", calculate_graph_cost(ub_solution, original_cost_matrix, n))
            show_result(exp_id, "lagrangean_relaxation:$iteration:check_point", "onetree cost on current costs ", calculate_graph_cost(one_tree_sol, current_cost_matrix, n))
            show_result(exp_id, "lagrangean_relaxation:$iteration:check_point", "christofides cost on current costs ", calculate_graph_cost(ub_solution, current_cost_matrix, n))
        end

        if optimality_gap == 0
            show_result(exp_id, "lagrangean_relaxation:$iteration:optimal_sol", "optimal_sol", "true" )
            show_result(exp_id, "lagrangean_relaxation:$iteration:optimal_sol", "optimality_gap", 0 )
            show_result(exp_id, "lagrangean_relaxation:$iteration:optimal_sol", "best_lower_bound", best_lower_bound)
            show_result(exp_id, "lagrangean_relaxation:$iteration:optimal_sol", "best_upper_bound", best_upper_bound)
            show_result(exp_id, "lagrangean_relaxation:$iteration:optimal_sol", "current_lower_bound", current_lower_bound)
            show_result(exp_id, "lagrangean_relaxation:$iteration:optimal_sol", "current_upper_bound", current_upper_bound)
            show_result(exp_id, "lagrangean_relaxation:$iteration:optimal_sol", "epsilon_strategy", epsilon_strategy)
            show_result(exp_id, "lagrangean_relaxation:$iteration:optimal_sol", "epsilon", epsilon)
            show_result(exp_id, "lagrangean_relaxation:$iteration:optimal_sol", "graph_cost", calculate_graph_cost(ub_solution, original_cost_matrix, n))
            show_result(exp_id, "lagrangean_relaxation:$iteration:optimal_sol", "min ub ", minimum(upper_bounds))
            show_result(exp_id, "lagrangean_relaxation:$iteration:optimal_sol", "max lb ", maximum(lower_bounds))
            show_result(exp_id, "lagrangean_relaxation:$iteration:optimal_sol", "min gap ", minimum(gaps))
            #return optimality_gap, best_lower_bound, best_upper_bound, best_ub_sol, best_lb_sol
            break
        end

        if optimality_gap <= gap_threshold
            show_result(exp_id, "lagrangean_relaxation:$iteration:acceptable_sol", "optimal_sol", "false" )
            show_result(exp_id, "lagrangean_relaxation:$iteration:acceptable_sol", "optimality_gap", optimality_gap )
            show_result(exp_id, "lagrangean_relaxation:$iteration:acceptable_sol", "best_lower_bound", best_lower_bound)
            show_result(exp_id, "lagrangean_relaxation:$iteration:acceptable_sol", "best_upper_bound", best_upper_bound)
            show_result(exp_id, "lagrangean_relaxation:$iteration:acceptable_sol", "current_lower_bound", current_lower_bound)
            show_result(exp_id, "lagrangean_relaxation:$iteration:acceptable_sol", "current_upper_bound", current_upper_bound)
            show_result(exp_id, "lagrangean_relaxation:$iteration:acceptable_sol", "epsilon_strategy", epsilon_strategy)
            show_result(exp_id, "lagrangean_relaxation:$iteration:acceptable_sol", "epsilon", epsilon)
            show_result(exp_id, "lagrangean_relaxation:$iteration:acceptable_sol", "graph_cost", calculate_graph_cost(ub_solution, original_cost_matrix, n))
            show_result(exp_id, "lagrangean_relaxation:$iteration:acceptable_sol", "min ub ", minimum(upper_bounds))
            show_result(exp_id, "lagrangean_relaxation:$iteration:acceptable_sol", "max lb ", maximum(lower_bounds))
            show_result(exp_id, "lagrangean_relaxation:$iteration:acceptable_sol", "min gap ", minimum(gaps))
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

        # updating epsilon
        if epsilon_strategy == "static"
            # does nothing, because static means 1 * epsilon
        elseif epsilon_strategy == "lbdecay"
            if (iteration - last_lower_boud_update) > (max_iterations / 20)
                epsilon = epsilon / 2
    
                if epsilon < epsilon_min
                    println("Stopping. epsilon < epsilon_min")
                    break
                end
            end
        elseif epsilon_strategy == "itdecay"
            epsilon = epsilon * 0.999 #0.9999999
        elseif epsilon_strategy == "itincrease"
            epsilon = epsilon * 1.001 #1.0000001
        else
            println("Invalid epsilon strategy \"$epsilon_strategy\".")
            break
        end

        if mi < epsilon_min
            println("Stopping. mi < epsilon_min.")
            break
        end
        
        #break
        save_step(exp_id,"lagrangean_relaxation:iterations","finish","iteration_$iteration")

    end

    save_step(exp_id,"lagrangean_relaxation:iterations","finish","step")

    #println("ubs ", upper_bounds)
    #println("lbs ", lower_bounds)
    #println("gaps ", gaps)
    show_result(exp_id, "lagrangean_relaxation", "min_gap ", minimum(gaps))
    show_result(exp_id, "lagrangean_relaxation", "min_upper_bound ", minimum(upper_bounds))
    show_result(exp_id, "lagrangean_relaxation", "max_lower_bound ", maximum(lower_bounds))
    show_result(exp_id, "lagrangean_relaxation", "optimality_gap ", optimality_gap)
    show_result(exp_id, "lagrangean_relaxation", "iterations ran ", total_iterations)

    save_step(exp_id,"lagrangean_relaxation","finish","method")
    return maximum(lower_bounds), minimum(upper_bounds)
end

#lagrangean_relaxation("1", "", 1000, 1.0)
