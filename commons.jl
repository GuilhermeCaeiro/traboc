
using Logging
using Dates
using DataFrames
using CSV


function new_folder(folder)
    mkdir(folder)
end

function get_ISODateTime()
    local_datatime = Dates.now()
    return Dates.format(local_datatime, "yyyymmddHHMMSS")
end

function new_experiment()
    exp_id = get_ISODateTime()
    #println(exp_id)
    new_folder(string("./work/",exp_id))
    df_steps = DataFrame(["timestamp" "step_id" "action" "target"], :auto)
    CSV.write(string("./work/", exp_id, "/experiment_steps.csv"), df_steps; append=true)
    df_results = DataFrame(["timestamp" "step_id" "key" "value"], :auto)
    CSV.write(string("./work/", exp_id, "/experiment_steps.csv"), df_results; append=true)
    return exp_id 
end

function save_step(exp_id,step_id,action,target)
    moment=now()
    timestamp = Dates.value(moment)
    df = DataFrame([timestamp step_id action target], :auto)
    CSV.write(string("./work/", exp_id, "/experiment_steps.csv"), df; append=true)
end

function save_result(exp_id,step_id,key,value)
    moment=now()
    timestamp = Dates.value(moment)
    df = DataFrame([timestamp step_id key value], :auto)
    CSV.write(string("./work/", exp_id, "/experiment_results.csv"), df; append=true)
end

function show_result(exp_id,step_id,key,value)
    show_info(string(key,value),"")
    save_result(exp_id,step_id,key,value)
end
