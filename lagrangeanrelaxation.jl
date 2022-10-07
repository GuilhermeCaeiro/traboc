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
include("heuristics.jl")

#
# Function lagrangean_relaxation
# obtains the lagrangean_relaxation with christofides method 
#  
# Parameters
# exp_params -> dictionary containing the experiment parameters
# Returns 
# void
function lagrangean_relaxation(exp_params)

    exp_id = exp_params["exp_id"]
    testdatafile = exp_params["testdatafile"]
    ub_algorithm = exp_params["ub_algorithm"]
    local_search = exp_params["local_search"]
    primal_input = exp_params["primal_input"]
    max_iterations = parse(Int64, exp_params["max_iterations"])
    gap_threshold = parse(Float64, exp_params["gap_threshold"])
    epsilon = parse(Float64, exp_params["epsilon"])
    mi_function = exp_params["mi_function"]
    epsilon_strategy = exp_params["epsilon_strategy"]
    epsilon_decay_interval = parse(Int64, exp_params["epsilon_decay_interval"])
    epsilon_decay_multiplier = parse(Float64, exp_params["epsilon_decay_multiplier"])

    check_point = parse(Int64,exp_params["check_point"])

    save_step(exp_id,"lagrangean_relaxation","start","algorithm")

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

    last_lb_or_epsilon_update = 0
    epsilon_min = 0.0001

    iteration_data = undef

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
    experiment_start_time = get_time_in_ms()
    is_optimal = false
    stop_condition = ""

    save_step(exp_id,"lagrangean_relaxation:iterations","start","step")

    for iteration = 1:max_iterations
        save_step(exp_id,"lagrangean_relaxation:iterations","start","iteration_$iteration")
        iteration_start_time = get_time_in_ms()

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

        # getting lower bound
        one_tree_sol = one_tree_graph(exp_id, iteration, current_cost_matrix, n)
        current_lower_bound = calculate_graph_cost(one_tree_sol, current_cost_matrix, n)
        current_lower_bound = current_lower_bound + 2 * sum(u)

        ub_solution = undef
        second_ub_solution = undef # only used if ub_algorithm is "christofidesandfactor2"
        current_upper_bound = Inf
        primal_input_cost_matrix = undef

        if primal_input == "lagrangean"
            primal_input_cost_matrix = current_cost_matrix
        elseif primal_input === "complementary"
            incidence = zeros(n, n)

            for edge in edges(one_tree_sol)
                s, d = src(edge), dst(edge)
            
                incidence[s, d] = 1
                incidence[d, s] = 1
            end

            primal_input_cost_matrix =  original_cost_matrix .* (1 .- incidence)
        else
            show_error("Invalid local search algorithm \"local_search\".")
        end

        # getting upper bound
        if ub_algorithm == "christofides"
            # christofides produces feasible solutions
            ub_solution = unite_and_hamilton(exp_id, iteration, primal_input_cost_matrix, n)
        elseif ub_algorithm == "factor2approximation"
            # factor 2 approximation produces feasible solutions
            ub_solution = factor_2_approximation(exp_id, iteration, primal_input_cost_matrix, n) 
        elseif ub_algorithm == "christofidesandfactor2"
            # 
            ub_solution = unite_and_hamilton(exp_id, iteration, primal_input_cost_matrix, n)
            second_ub_solution = factor_2_approximation(exp_id, iteration, primal_input_cost_matrix, n) 
        elseif ub_algorithm == "nearestneighbor"
            ub_solution = nearest_neighbor(primal_input_cost_matrix, n)
        elseif ub_algorithm == "farthestinsertion"
            ub_solution = farthest_insertion(primal_input_cost_matrix, n)
        else
            show_error("Invalid upper bound algorithm \"$ub_algorithm\".")
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

        if local_search == "2opt"
            #ub_solution, current_upper_bound = two_opt(ub_solution, primal_input_cost_matrix, n)
            ub_solution, current_upper_bound = two_opt(ub_solution, original_cost_matrix, n)
        elseif local_search == "none"
            # does nothing
        else
            show_error("Invalid local search algorithm \"local_search\".")
        end

        

        if current_upper_bound < best_upper_bound
            best_upper_bound = current_upper_bound
            best_ub_sol = ub_solution
        end

        if current_lower_bound > best_lower_bound
            best_lower_bound = current_lower_bound
            best_lb_sol = one_tree_sol
            last_lb_or_epsilon_update = iteration
        end

        #optimality_gap = (current_upper_bound - current_lower_bound) / current_upper_bound
        optimality_gap = (best_upper_bound - best_lower_bound) / best_upper_bound
        
        # lagrangian vars update
        edges_per_vertex = count_vertex_edges(one_tree_sol, n)

        G = 2 .- edges_per_vertex
        denominator = sum(G .^ 2)

        #mi function selector
        if mi_function == "current" 
            mi = epsilon * (current_upper_bound - current_lower_bound) / denominator
        elseif mi_function == "best"
            mi = epsilon * (best_upper_bound - current_lower_bound) / denominator
        elseif mi_function == "5pct"
            mi = epsilon * (1.05 * best_upper_bound - current_lower_bound) / denominator
        elseif mi_function == "1pct"
            mi = epsilon * (1.01 * best_upper_bound - current_lower_bound) / denominator
        else
            @error "Mi function option não informada."
        end

        current_u = u
        u = u + mi * G

        #draw(PDF(string("christofides_", iteration, ".pdf"), 16cm, 16cm), gplot(christofides_sol, nodelabel=1:nv(one_tree_sol)))
        #draw(PDF(string("onetree_", iteration, ".pdf"), 16cm, 16cm), gplot(one_tree_sol, nodelabel=1:nv(one_tree_sol)))

        push!(gaps, optimality_gap)
        push!(upper_bounds, current_upper_bound)
        push!(lower_bounds, current_lower_bound)

        iteration_finish_time = get_time_in_ms()

        iteration_data = Dict(
            "iteration" => iteration,
            "current_upper_bound" => current_upper_bound,
            "current_lower_bound" => current_lower_bound,
            "best_lower_bound" => best_lower_bound,
            "best_upper_bound" => best_upper_bound,
            "epsilon_strategy" => epsilon_strategy,
            "epsilon" => epsilon,
            "epsilon_decay_interval" => epsilon_decay_interval, 
            "epsilon_decay_multiplier" => epsilon_decay_multiplier,
            "optimality_gap" => optimality_gap,
            "experiment_id" => exp_id,
            "instance" => testdatafile,
            "upper_bound_algorithm" => ub_algorithm,
            "local_search" => local_search,
            "primal_input" => primal_input,
            "max_iterations" => max_iterations,
            "gap_threshold" => gap_threshold,
            "mi_function" => mi_function, 
            "mi" => mi,
            "current_u" => current_u,
            "new_u" => u,
            "G" => G,
            "is_optimal" => is_optimal,
            "stop_condition" => stop_condition,
            "iteration_start_time" => iteration_start_time,
            "iteration_finish_time" => iteration_finish_time,
            "total_iteration_time" => iteration_finish_time - iteration_start_time,
            "experiment_start_time" => experiment_start_time,
            "experiment_finish_time" => undef,
            "total_experiment_time" => undef
        )

        # stopping criteria
        if optimality_gap == 0
            is_optimal = true
            stop_condition = "Optimal solution found. UB = LB."
        elseif optimality_gap <= gap_threshold
            is_optimal = true
            stop_condition = "Acceptable solution found. Gap threshold reached."
        elseif denominator == 0
            feasible = check_hamiltonian_cycle(one_tree_sol)

            # if the solution is feasible and the denominator is zero, that solution is optimal
            if feasible
                is_optimal = true
                stop_condition = "Optimal solution found. The denominator do the subgradient is zero and the solution is feasible."
            else
                is_optimal = false
                stop_condition = "Lower bound found. Gap threshold reached. The denominator do the subgradient is zero, but the solution is unfeasible."
            end
        elseif mi < epsilon_min
            stop_condition = "Lower bound found. mi too small (mi < epsilon_min)."
        elseif epsilon < epsilon_min
            stop_condition = "Lower bound found. epsilon too small (epsilon < epsilon_min)."
        elseif iteration == max_iterations
            stop_condition = "Maximum number of iterations reached"
        end

        iteration_data["is_optimal"] = is_optimal
        iteration_data["stop_condition"] = stop_condition

        stop_iterating = false
        if iteration_data["stop_condition"] != ""
            stop_iterating = true
            iteration_data["experiment_finish_time"] = iteration_finish_time
            iteration_data["total_experiment_time"] = iteration_finish_time - experiment_start_time
        end

        print_iteration_data(
            iteration_data, 
            cli_only_checkpoint = true, 
            checkpoint = check_point, 
            to_csv = true, 
            output_csv = string("./work/", exp_id, "/experiment_iterations.csv")
        )

        if stop_iterating
            break
        end

        # updating epsilon
        if epsilon_strategy == "static"
            # does nothing, because static means 1 * epsilon
        elseif epsilon_strategy == "lbdecay1"
            if (iteration - last_lb_or_epsilon_update) > epsilon_decay_interval #(max_iterations / 20)
                epsilon = epsilon * epsilon_decay_multiplier
                last_lb_or_epsilon_update = iteration
            end
        elseif epsilon_strategy == "lbdecay2"
            if (iteration - last_lb_or_epsilon_update) > epsilon_decay_interval #(max_iterations / 20)
                epsilon = epsilon * epsilon_decay_multiplier
            end
        elseif epsilon_strategy == "itdecay"
            epsilon = epsilon * epsilon_decay_multiplier #0.9999999
        else
            show_error("Invalid epsilon strategy \"$epsilon_strategy\".")
            break
        end
        
        #break
        save_step(exp_id,"lagrangean_relaxation:iterations","finish","iteration_$iteration")

    end

    save_step(exp_id,"lagrangean_relaxation:iterations","finish","step")

    @debug "ubs ", upper_bounds
    @debug "lbs ", lower_bounds
    @debug "gaps ", gaps

    results = Dict(
        "exp_id" => exp_id,
        "method" => "lagrangean_relaxation",
        "status" => "-",
        "objective_value" => "-",
        "min_gap" => minimum(gaps),
        "current_lower_bound" => current_lower_bound,
        "best_lower_bound" => best_lower_bound,
        "min_upper_bound" => minimum(upper_bounds),
        "max_lower_bound" => maximum(lower_bounds),
        "optimality_gap" => optimality_gap,
        "iterations_ran" => total_iterations,
        "is_optimal" => is_optimal,
        "stop_condition" => stop_condition
    )

    for (key, value) in results
        show_result(exp_id, "lagrangean_relaxation", string(key," "), string(value))
    end

    save_step(exp_id,"lagrangean_relaxation","finish","algorithm")
    return results
end

#lagrangean_relaxation("1", "", 1000, 1.0)
