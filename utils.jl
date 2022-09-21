
using Logging

function set_logging(loglevel,experiment_id)

    # Open a textfile for writing
    io = open(string("./logs/logs_",experiment_id,".txt"), "w+")

    # Create a simple logger
    logger = SimpleLogger(io, loglevel)

    # Set the global logger to logger
    global_logger(logger)
 
end

function show_info(info)
    show_info(info,"")
end

function show_info(info, detail)
    @info string(info,detail)
    println(info, detail)
end