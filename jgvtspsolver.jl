
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

function save_params(exp_id, strategy, testdatafile, max_iterations, gap_threshold, epsilon, mi_option, check_point)
    df = DataFrame([exp_id strategy testdatafile max_iterations gap_threshold epsilon mi_option check_point], :auto)
    CSV.write("./work/experiments.csv", df; append=true)
end

function main(args)
 
    exp_id = new_experiment() 
    
    set_logging(Logging.Info,exp_id)

    if length(args) < 3
        @error "Por favor informe os parâmetros para execução: testdatafile max_iterations epsilon"
    end 

    strategy = args[1]
    if strategy == "-help"
        show_info( "jgvtspsolver.jl <christofides|gurobi> <testdatafile> <max_iterations> <gap_threshold> <epsilon> <current|best|5pct|1pct> [<check_point>]" )
        show_info( "<christofides|gurobi> solver engine to be used in experiment" )
        show_info( "<testdatafile> path to test data file to be used" )
        show_info( "<max_iterations> maximn number of iterations" )
        show_info( "<gap_threshold> gap gap_threshold to be used in lagrangean_relaxation/christofides" )
        show_info( "<epsilon> epsilon value to be used in lagrangean_relaxation/christofides" )
        show_info( "<current|best|5pct|1pct> mi_function to be used to calculate mi parameter" )
        show_info( "[<check_point>] iteration check-point to show progress (default: max_iterations)" )
        exit()
    end

    testdatafile = args[2]
    max_iterations = parse(Int64,args[3])
    gap_threshold = parse(Float64,args[4])
    epsilon = parse(Float64,args[5])
    mi_option = args[6]
    check_point = max_iterations
    if length(args) > 6
        check_point = parse(Int64,args[7])
    end

    show_info("********************************************************************")
    show_info("JGV Solver For Travelling Sallesman problem ")
    show_info("********************************************************************")
    show_info("Trabalho para a disciplina de Otimização Combinatória - 2022p2")
    show_info("João Pedro")
    show_info("Guilherme Caeiro")
    show_info("Victor Xavier")
    show_info("********************************************************************")
    show_info("Exeperimento: ", exp_id)
    show_info("Parâmetros: " )
    show_info("strategy: ", strategy)
    show_info("testdatafile: ", testdatafile)
    show_info("max_iterations: ", max_iterations)
    show_info("gap_threshold: ", gap_threshold)
    show_info("epsilon: ", epsilon)
    show_info("mi_function: ", mi_option)
    show_info("check_point: ", check_point)

    save_step(exp_id,"main","start","experiment")
    save_params(exp_id,strategy,testdatafile,max_iterations,gap_threshold,epsilon,mi_option,check_point)

    if strategy in ["christofides", "factor2approximation", "christofidesandfactor2"]
        lagrangean_relaxation(exp_id, testdatafile, strategy, max_iterations, gap_threshold, epsilon, mi_option, check_point)
    elseif strategy == "gurobi"
        lazy(testdatafile)
    elseif strategy == "lazyandchris"
        lower, upper = lagrangean_relaxation(exp_id, testdatafile, "christofides", max_iterations, gap_threshold, epsilon, mi_option, check_point)
        lazy(testdatafile, lower, upper)
    else
        @error "Por favor informe uma estratégia de execução válida: cristofides ou solvers"
    end

    save_step(exp_id,"main","finish","experiment")

end

main(ARGS)