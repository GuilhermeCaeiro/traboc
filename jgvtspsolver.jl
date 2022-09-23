
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

function save_params(exp_id, strategy, testdatafile,max_iterations,epsilon)
    df = DataFrame([exp_id strategy testdatafile max_iterations epsilon], :auto)
    CSV.write("./work/experiments.csv", df; append=true)
end

function main(args)
 
    exp_id = new_experiment() 
    
    set_logging(Logging.Debug,exp_id)

    if length(args) < 3
        @error "Por favor informe os parâmetros para execução: testdatafile max_iterations epsilon"
    end 

    strategy = args[1]
    testdatafile = args[2]
    max_iterations = parse(Int64,args[3])
    epsilon = parse(Float64,args[4])

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
    show_info("epsilon: ", epsilon)

    save_step(exp_id,"MAIN","START")
    save_params(exp_id,strategy,testdatafile,max_iterations,epsilon)

    if strategy == "christofides"
        lagrangean_relaxation(exp_id,testdatafile,max_iterations,epsilon)
    elseif strategy == "gurobi"
        @info "Implementar chamada a solvers aqui..."
    else
        @error "Por favor informe uma estratégia de execução válida: cristofides ou solvers"
    end

end

main(ARGS)