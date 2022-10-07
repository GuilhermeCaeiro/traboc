
#
# JGV Solver For Travelling Sallesman problem
# 
# Trabalho para a disciplina de Otimização Combinatória - 2022p2
#
# Grupo:    João Pedro
#           Guilherme Caieiros
#           Victor Xavier
#

using Logging
using DataFrames
using CSV

include("setupenv.jl")
include("commons.jl")
include("logs.jl")
include("utils.jl")
include("lagrangeanrelaxation.jl")
include("lazytsp.jl")

function save_params(exp_params)
    parameters = ["exp_id" "strategy" "testdatafile" "max_iterations" "ub_algorithm" "gap_threshold" "epsilon" "mi_function" "epsilon_strategy" "epsilon_decay_interval" "epsilon_decay_multiplier" "local_search" "primal_input" "check_point" "log_level" "rac_strategy" "pre_solver" "lower_bound" "upper_bound"]
    params_data = [ string(get!(exp_params,param_key,"-")) for param_key in parameters ]
    if !isfile("./work/experiments.csv")
        hdf = DataFrame(parameters, :auto)
        CSV.write("./work/experiments.csv", hdf; append=true)    
    end
    df = DataFrame(params_data, :auto)
    CSV.write("./work/experiments.csv", df; append=true)
end

function save_experiment_summary(results)
    summary_items = ["exp_id" "method" "status" "objective_value" "min_gap" "current_lower_bound" "best_lower_bound" "min_upper_bound" "max_lower_bound" "optimality_gap" "iterations_ran" "is_optimal" "stop_condition"]
    summary_data = [ string(get!(results,key,"-")) for key in summary_items ]
     
    if !isfile("./work/experiments_summary.csv")
        hdf = DataFrame(summary_items, :auto)
        CSV.write("./work/experiments_summary.csv", hdf; append=true)    
    end
    df = DataFrame(summary_data, :auto)
    CSV.write("./work/experiments_summary.csv", df; append=true)
end


function main(args)
 
    exp_id = new_experiment() 
    exp_params = Dict{String}{String}()
    
    set_logging(Logging.Info,exp_id)

    show_info("********************************************************************")
    show_info("JGV Solver for the Travelling Salesman Problem ")
    show_info("********************************************************************")
    show_info("Trabalho para a disciplina de Otimização Combinatória - 2022p2")
    show_info("João Pedro Souza - joaosouza@cos.ufrj.br")
    show_info("Guilherme Caeiro de Mattos - mattosgc@cos.ufrj.br")
    show_info("Victor de Almeida Xavier - victorx@cos.ufrj.br")
    show_info("********************************************************************")

    strategy = args[1]
    if strategy == "-help"
        show_info( " " )
        show_info( "Lagrangean-relaxation based strategies:" )
        show_info( "jgvtspsolver.jl <christofides|factor2approximation|christofidesandfactor2> <testdatafile> <max_iterations> <gap_threshold> <epsilon> <current|best|5pct|1pct> <static|lbdecay|itdecay|itincrease> <epsilon_decay_interval> <epsilon_decay_multiplier> <2opt|none> <lagrangean|complementary> [<check_point>] [<debug|info|error>]" )
        show_info( "<christofides|factor2approximation|christofidesandfactor2|nearestneighbor|farthestinsertion> strategy engine to be used in experiment." )
        show_info( "<testdatafile> path to test data file to be used" )
        show_info( "<max_iterations> maximn number of iterations" )
        show_info( "<gap_threshold> gap gap_threshold to be used in lagrangean_relaxation/christofides" )
        show_info( "<epsilon> epsilon value to be used in lagrangean_relaxation/christofides" )
        show_info( "<current|best|5pct|1pct> mi_function to be used to calculate mi parameter" )
        show_info( "<static|lbdecay1|lbdecay2|itdecay|itincrease> strategy to update epsilon" )
        show_info( "<epsilon_decay_interval> number of iterations without update to the lb that must pass before epsilon is updated" )
        show_info( "<epsilon_decay_multiplier> epsilon decay multiplier to be used by lbdecay and itdecay (can be > 1, but should be ideally between 0 and 1." )
        show_info( "<2opt|none> local search algorithm" )
        show_info( "<lagrangean|complementary> costs to be provided to the primal algorithms" )
        show_info( "[<check_point>] iteration check-point to show progress (default: max_iterations)" )
        show_info( "[<debug|info|error>] log level (default: info)" )
        show_info( " " )
        show_info( "Solver based strategies:" )
        show_info( "jgvtspsolver.jl <gurobi|glpk> <testdatafile> <max_iterations> [<combs|secs|both>] [<none|christofides|factor2approximation|christofidesandfactor2>[,<gap_threshold>,<epsilon>,<current|best|5pct|1pct>,<static|lbdecay|itdecay|itincrease>,<epsilon_decay_interval>,<epsilon_decay_multiplier>]] [<check_point>] [<debug|info|error>]" )
        show_info( "<gurobi|glpk> solver engine to be used in experiment." )
        show_info( "<testdatafile> path to test data file to be used" )
        show_info( "<max_iterations> maximn number of iterations" )
        show_info( "[<combs|secs|both>] relax-and-cut option (default: both)" )
        show_info( "[<none|christofides|factor2approximation|christofidesandfactor2>] lagrangean_relaxation preprocess (default: none)" )
        show_info( "[<check_point>] iteration check-point to show progress (default: max_iterations)" )
        show_info( "[<debug|info|error>] log level (default: info)" )
        show_info( " " )
        exit()
    end

    if length(args) < 3
        show_error("Por favor informe os parâmetros para execução ou -help para a lista completa de opções.")
        exit()
    end 

    exp_params["exp_id"] = exp_id
    exp_params["strategy"] = strategy
    exp_params["testdatafile"] = args[2]
    exp_params["max_iterations"] = args[3] 
    exp_params["check_point"] = args[3]
    exp_params["log_level"] = "info"             

    show_info("Experimento: ", exp_id)
    show_info("Parâmetros: " )
    show_info("strategy: ", strategy)
    show_info("testdatafile: ", exp_params["testdatafile"])
    show_info("max_iterations: ", exp_params["max_iterations"])

    save_step(exp_id,"main","start","experiment")

    if strategy in ["christofides", "factor2approximation", "christofidesandfactor2", "nearestneighbor", "farthestinsertion"]
        if length(args) < 7
            @error "Por favor informe os parâmetros para execução do experimento $strategy ou -help para a lista completa de opções."
            exit()
        end     
        exp_params["ub_algorithm"] = strategy 
        exp_params["gap_threshold"] = args[4] 
        if !check_float_param( exp_params["gap_threshold"], 0.000001, 0.1 )
            show_error("O parâmetro gap_threshold precisa estar entre os valores 0.000001 e 0.1" )
            exit()
        end

        exp_params["epsilon"] = args[5] 
        if !check_int_param( exp_params["epsilon"], 1, 2 )
            show_error("O parâmetro epsilon precisa estar entre os valores inteiros 1 e 2" )
            exit()
        end

        exp_params["mi_function"] = args[6]
        if !check_opt_param( exp_params["mi_function"], "current|best|5pct|1pct" )
            show_error("O parâmetro mi_function precisa ser uma das opções: current|best|5pct|1pct" )
            exit()
        end

        exp_params["epsilon_strategy"] = args[7] 
        if !check_opt_param( exp_params["epsilon_strategy"], "static|lbdecay1|lbdecay2|itdecay|itincrease" )
            show_error("O parâmetro epsilon_strategy precisa ser uma das opções: static|lbdecay1|lbdecay2|itdecay|itincrease" )
            exit()
        end

        exp_params["epsilon_decay_interval"] = args[8]
        if !check_int_param( exp_params["epsilon_decay_interval"], 1, 1000 )
            show_error("O parâmetro epsilon_decay_interval precisa estar entre os valores 1 e 1000")
            exit()
        end

        exp_params["epsilon_decay_multiplier"] = args[9] 
        if !check_float_param( exp_params["epsilon_decay_multiplier"], 0.000001, 2.0)
            show_error("O parâmetro epsilon_decay_multiplier precisa estar entre os valores 0.000001 e 1, mas pode receber valores até 2" )
            exit()
        end

        exp_params["local_search"] = args[10] 
        if !check_opt_param( exp_params["local_search"], "2opt|none" )
            show_error("O parâmetro local_search precisa ser uma das opções: 2opt|none")
            exit()
        end

        exp_params["primal_input"] = args[11] 
        if !check_opt_param( exp_params["primal_input"], "lagrangean|complementary" )
            show_error("O parâmetro primal_input precisa ser uma das opções: lagrangean|complementary")
            exit()
        end

        show_info("gap_threshold: ", exp_params["gap_threshold"])
        show_info("epsilon: ", exp_params["epsilon"])
        show_info("mi_function: ", exp_params["mi_function"])
        show_info("epsilon_strategy: ", exp_params["epsilon_strategy"])
        show_info("epsilon_decay_interval: ", exp_params["epsilon_decay_interval"])
        show_info("epsilon_decay_multiplier: ", exp_params["epsilon_decay_multiplier"])
        show_info("local_search: ", exp_params["local_search"])
        show_info("primal_input: ", exp_params["primal_input"])
    
        if length(args) > 11
            exp_params["check_point"] = args[12] 
            show_info("check_point: ", exp_params["check_point"])
        end
        if length(args) > 12
            exp_params["log_level"] = args[13]             
            show_info("log_level: ", exp_params["log_level"])
            set_loglevel(exp_params["log_level"],exp_id)
        end

        show_info("********************************************************************")
        show_info("Running experiment $exp_id...")

        @debug "lagrangean_relaxation..." 
        #execute lagrangean_relaxation with exp_params 
        results = lagrangean_relaxation(exp_params)
        @debug "lagrangean_relaxation...done"

        @debug "saving summary..." 
        save_experiment_summary(results)
        @debug "saving summary...done" 

        show_info("Running experiment $exp_id...done")
        show_info("********************************************************************")

    elseif strategy in [ "gurobi", "glpk" ]

        if length(args) < 3
            @error "Por favor informe os parâmetros para execução do experimento $strategy ou -help para a lista completa de opções."
            exit()
        end 

        exp_params["rac_strategy"] = "both"
        if length(args) > 3
            exp_params["rac_strategy"] = args[4] 

            if !check_opt_param( exp_params["rac_strategy"], "combs|secs|both" )
                show_error("O parâmetro rac_strategy precisa estar entre as opções: combs|secs|both" )
                exit()
            end

            show_info("rac_strategy: ", exp_params["rac_strategy"])
        end

        exp_params["pre_solver"] = "none"
        if length(args) > 4
            pre_solver_params = split(args[5],",")
            if pre_solver_params[1] != "none"
                if length(pre_solver_params) < 7
                    show_error("Por favor informe os parâmetros adicionais para o pre-solver seperados por vírgula: <gap_threshold>,<epsilon>,<current|best|5pct|1pct>,<static|lbdecay|itdecay|itincrease>,<epsilon_decay_interval>,<epsilon_decay_multiplier>")
                    exit()
                end

                exp_params["pre_solver"] = pre_solver_params[1]
                pre_solver = exp_params["pre_solver"]

                exp_params["gap_threshold"] = pre_solver_params[2]
                if !check_float_param( exp_params["gap_threshold"], 0.000001, 0.1 )
                    show_error("O parâmetro $pre_solver:gap_threshold precisa estar entre os valores 0.000001 e 0.1" )
                    exit()
                end
        
                exp_params["epsilon"] = pre_solver_params[3]
                if !check_int_param( exp_params["epsilon"], 1, 2 )
                    show_error("O parâmetro $pre_solver:epsilon precisa estar entre os valores inteiros 1 e 2" )
                    exit()
                end
        
                exp_params["mi_function"] = pre_solver_params[4]
                if !check_opt_param( exp_params["mi_function"], "current|best|5pct|1pct" )
                    show_error("O parâmetro $pre_solver:mi_function precisa ser uma das opções: current|best|5pct|1pct" )
                    exit()
                end
        
                exp_params["epsilon_strategy"] = pre_solver_params[5]
                if !check_opt_param( exp_params["epsilon_strategy"], "static|lbdecay1|lbdecay2|itdecay|itincrease" )
                    show_error("O parâmetro $pre_solver:epsilon_strategy precisa ser uma das opções: static|lbdecay1|lbdecay2|itdecay|itincrease" )
                    exit()
                end
        
                exp_params["epsilon_decay_interval"] = pre_solver_params[6]
                if !check_float_param( exp_params["epsilon_decay_interval"], 0.00001, 0.9 )
                    show_error("O parâmetro $pre_solver:epsilon_strategy precisa estar entre os valores 0.000001 e 0.9" )
                    exit()
                end
        
                exp_params["epsilon_decay_multiplier"] = pre_solver_params[7]
                if !check_int_param( exp_params["epsilon_decay_multiplier"], 2, 50 )
                    show_error("O parâmetro $pre_solver:epsilon_decay_multiplier precisa estar entre os valores 2 e 50" )
                    exit()
                end
        
                show_info("pre_solver: ", exp_params["pre_solver"])
                show_info("$pre_solver:gap_threshold: ", exp_params["gap_threshold"])
                show_info("$pre_solver:epsilon: ", exp_params["epsilon"])
                show_info("$pre_solver:mi_function: ", exp_params["mi_function"])
                show_info("$pre_solver:epsilon_strategy: ", exp_params["epsilon_strategy"])
                show_info("$pre_solver:epsilon_decay_interval: ", exp_params["epsilon_decay_interval"])
                show_info("$pre_solver:epsilon_decay_multiplier: ", exp_params["epsilon_decay_multiplier"])

            else
                exp_params["pre_solver"] = args[5]
                show_info("pre_solver: ", exp_params["pre_solver"])
            end
            
        end

        exp_params["check_point"] = exp_params["max_iterations"] 
        if length(args) > 5
            exp_params["check_point"] = args[6] 
            show_info("check_point: ", exp_params["check_point"])
        end

        exp_params["log_level"] = "info" 
        if length(args) > 6
            exp_params["log_level"] = args[7] 
            show_info("log_level: ", exp_params["check_point"])
            set_loglevel(exp_params["log_level"],exp_id)
        end


        show_info("********************************************************************")
        show_info("Running experiment $exp_id...")

        #pre-solvers using lagrangean_relaxation
        exp_params["lower_bound"] = string(-Inf)
        exp_params["upper_bound"] = string(Inf)

        if exp_params["pre_solver"] in ["christofides", "factor2approximation", "christofidesandfactor2"]
            exp_params["ub_algorithm"] = exp_params["pre_solver"]
            @debug "lagrangean_relaxation as pre solver..." 
            results = lagrangean_relaxation(exp_params)            
            @debug "lagrangean_relaxation as pre solver...done"
            exp_params["lower_bound"] = results["lower_bound"]
            exp_params["upper_bound"] = results["upper_bound"]
        end

        #execute solver with lazy constraints passing exp_params
        @debug "lazy..." 
        results = lazy(exp_params)
        @debug "lazy...done"

        @debug "saving summary..." 
        save_experiment_summary(results)
        @debug "saving summary...done" 

        show_info("Running experiment $exp_id...done")
        show_info("********************************************************************")
    
    else
        show_error("Por favor informe uma estratégia de execução válida: christofides, factor2approximation, christofidesandfactor2, gurobi, glpk.")
    end

    save_params(exp_params)

    save_step(exp_id,"main","finish","experiment")

end

main(ARGS)