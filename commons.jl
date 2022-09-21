
using Logging
using Dates


function new_folder(folder)
    mkdir(folder)
end

function get_ISODateTime()
    local_datatime = Dates.now()
    return Dates.format(local_datatime, "yyyymmddHHMMSS")
end

function new_experiment()
    exp_id = get_ISODateTime()
    new_folder(string("./work/",exp_id))
    df = DataFrame(["step_id" "action" "timestamp"], :auto)
    CSV.write(string("./work/",exp_id,"/experiment_steps.csv"), df; append=true)
    return exp_id 
end

function save_step(exp_id,step_id,action)
    moment=now()
    timestamp = Dates.value(moment)
    df = DataFrame([step_id action timestamp], :auto)
    CSV.write(string("./work/",exp_id,"/experiment_steps.csv"), df; append=true)
end