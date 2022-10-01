
using Logging

function set_logging(loglevel,experiment_id)

    # Open a textfile for writing
    io = open(string("./logs/logs_",experiment_id,".txt"), "w+")

    # Create a simple logger
    logger = SimpleLogger(io, loglevel)

    # Set the global logger to logger
    global_logger(logger)
 
end

function set_loglevel(log_str,exp_id)
    if log_str == "debug"
        set_logging(Logging.Debug,exp_id)        
    elseif log_str == "info"
        set_logging(Logging.Info,exp_id)        
    elseif log_str == "error"
        set_logging(Logging.Info,exp_id)        
    end
end

function show_info(info)
    show_info(info,"")
end

function show_info(info, detail)
    @info string(info,detail)
    println(info, detail)
end

function show_error(error_msg)
    @error error_msg
    show_info(error_msg)
end
