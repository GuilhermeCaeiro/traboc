
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
    parameters = ["exp_id" "strategy" "testdatafile" "max_iterations" "ub_algorithm" "gap_threshold" "epsilon" "mi_function" "epsilon_strategy" "check_point" "log_level" "use_combs" "pre_solver" "lower_bound" "upper_bound"]
    params_data = [ string(get!(exp_params,param_key,"-")) for param_key in parameters ]
    if !isfile("./work/experiments.csv")
        hdf = DataFrame(parameters, :auto)
        CSV.write("./work/experiments.csv", hdf; append=true)    
    end
    df = DataFrame(params_data, :auto)
    CSV.write("./work/experiments.csv", df; append=true)
end

function main(args)
 
    exp_id = new_experiment() 
    exp_params = Dict{String}{String}()
    
    set_logging(Logging.Info,exp_id)

    show_info("********************************************************************")
    show_info("JGV Solver For Travelling Sallesman problem ")
    show_info("********************************************************************")
    show_info("Trabalho para a disciplina de Otimização Combinatória - 2022p2")
    show_info("João Pedro")
    show_info("Guilherme Caeiro")
    show_info("Victor Xavier")
    show_info("********************************************************************")

    strategy = args[1]
    if strategy == "-help"
        show_info( " " )
        show_info( "Lagrangean-relaxation based strategies:" )
        show_info( "jgvtspsolver.jl <christofides|factor2approximation|christofidesandfactor2> <testdatafile> <max_iterations> <gap_threshold> <epsilon> <current|best|5pct|1pct> <static|lbdecay|itdecay|itincrease> [<check_point>] [<debug|info|error>]" )
        show_info( "<christofides|factor2approximation|christofidesandfactor2> strategy engine to be used in experiment." )
        show_info( "<testdatafile> path to test data file to be used" )
        show_info( "<max_iterations> maximn number of iterations" )
        show_info( "<gap_threshold> gap gap_threshold to be used in lagrangean_relaxation/christofides" )
        show_info( "<epsilon> epsilon value to be used in lagrangean_relaxation/christofides" )
        show_info( "<current|best|5pct|1pct> mi_function to be used to calculate mi parameter" )
        show_info( "<static|lbdecay|itdecay|itincrease> strategy to update epsilon" )
        show_info( "[<check_point>] iteration check-point to show progress (default: max_iterations)" )
        show_info( "[<debug|info|error>] log level (default: info)" )
        show_info( " " )
        show_info( "Solver based strategies:" )
        show_info( "jgvtspsolver.jl <gurobi|glpk> <testdatafile> <max_iterations> [<use_combs>] [<none|christofides|factor2approximation|christofidesandfactor2>] [<check_point>] [<debug|info|error>]" )
        show_info( "<gurobi|glpk> solver engine to be used in experiment." )
        show_info( "<testdatafile> path to test data file to be used" )
        show_info( "<max_iterations> maximn number of iterations" )
        show_info( "[<use_combs>] true indicates to use comb (relax and cut) (default: false)" )
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

    show_info("Exeperimento: ", exp_id)
    show_info("Parâmetros: " )
    show_info("strategy: ", strategy)
    show_info("testdatafile: ", exp_params["testdatafile"])
    show_info("max_iterations: ", exp_params["max_iterations"])

    save_step(exp_id,"main","start","experiment")

    if strategy in ["christofides", "factor2approximation", "christofidesandfactor2"]
        if length(args) < 7
            @error "Por favor informe os parâmetros para execução do experimento $strategy ou -help para a lista completa de opções."
            exit()
        end     
        exp_params["ub_algorithm"] = strategy 
        exp_params["gap_threshold"] = args[4] 
        exp_params["epsilon"] = args[5] 
        exp_params["mi_option"] = args[6] 
        exp_params["epsilon_strategy"] = args[7] 

        show_info("gap_threshold: ", exp_params["gap_threshold"])
        show_info("epsilon: ", exp_params["epsilon"])
        show_info("mi_function: ", exp_params["mi_option"])
        show_info("epsilon strategy: ", exp_params["epsilon_strategy"])
    
        if length(args) > 7
            exp_params["check_point"] = args[8] 
            show_info("check_point: ", exp_params["check_point"])
        end
        if length(args) > 8
            exp_params["log_level"] = args[9]             
            show_info("log_level: ", exp_params["check_point"])
            set_loglevel(exp_params["log_level"],exp_id)
        end

        show_info("********************************************************************")
        show_info("Running experiment $exp_id...")

        @debug "lagrangean_relaxation..." 
        #execute lagrangean_relaxation with exp_params 
        lagrangean_relaxation(exp_params)
        @debug "lagrangean_relaxation...done"

        show_info("Running experiment $exp_id...done")
        show_info("********************************************************************")

    elseif strategy in [ "gurobi", "glpk" ]

        if length(args) < 3
            @error "Por favor informe os parâmetros para execução do experimento $strategy ou -help para a lista completa de opções."
            exit()
        end 

        exp_params["use_combs"] = "false"
        if length(args) > 3
            exp_params["use_combs"] = args[4] 
            show_info("use_combs: ", exp_params["use_combs"])
        end

        exp_params["pre_solver"] = "none"
        if length(args) > 4
            exp_params["pre_solver"] = args[5] 
            show_info("pre_solver: ", exp_params["pre_solver"])
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

        if exp_params["pre_solver"] != "none"
            exp_params["ub_algorithm"] = exp_params["pre_solver"]
            @debug "lagrangean_relaxation as pre solver..." 
            lower, upper = lagrangean_relaxation(exp_params)
            @debug "lagrangean_relaxation as pre solver...done"
            exp_params["lower_bound"] = string(lower)
            exp_params["upper_bound"] = string(upper)
        end

        #execute solver with lazy constraints passing exp_params
        @debug "lazy..." 
        lazy(exp_params)
        @debug "lazy...done"

        show_info("Running experiment $exp_id...done")
        show_info("********************************************************************")
    
    else
        show_error("Por favor informe uma estratégia de execução válida: christofides, factor2approximation, christofidesandfactor2, gurobi, glpk.")
    end

    save_params(exp_params)

    save_step(exp_id,"main","finish","experiment")

end

main(ARGS)