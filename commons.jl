using Logging
using Dates
using DataFrames
using CSV

function get_next_id()
    local_datatime = Dates.now()
    return Dates.format(local_datatime, "HHMMSSsss")
end

function new_folder(folder)
    mkdir(folder)
end

function get_ISODateTime()
    local_datatime = Dates.now()
    return Dates.format(local_datatime, "yyyymmddHHMMSS")
end

function get_time_in_ms()
    return convert(Dates.Millisecond, Dates.now())
end

function new_experiment()
    exp_id = get_ISODateTime()
    #println(exp_id)
    new_folder(string("./work/",exp_id))
    df_steps = DataFrame(["timestamp" "step_id" "action" "target"], :auto)
    CSV.write(string("./work/", exp_id, "/experiment_steps.csv"), df_steps; append=true)
    df_results = DataFrame(["timestamp" "step_id" "key" "value"], :auto)
    CSV.write(string("./work/", exp_id, "/experiment_results.csv"), df_results; append=true)
    return exp_id 
end

function save_step(exp_id,step_id,action,target)
    moment=now()
    timestamp = Dates.value(moment)
    df = DataFrame([timestamp step_id action target], :auto)
    CSV.write(string("./work/", exp_id, "/experiment_steps.csv"), df; append=true)
    return
end

function save_result(exp_id,step_id,key,value)
    moment=now()
    timestamp = Dates.value(moment)
    df = DataFrame([timestamp step_id key value], :auto)
    CSV.write(string("./work/", exp_id, "/experiment_results.csv"), df; append=true)
    return
end

function show_result(exp_id,step_id,key,value)
    show_info(string(key,value),"")
    save_result(exp_id,step_id,key,value)
end

function format_cli_print(data_dict)
    iteration = data_dict["iteration"]
    template = "#################### Iteration $iteration ####################\n[DATA]############################################################\n"
    data = ""

    for key in sort(collect(keys(data_dict)))
        if key != "iteration"
            value = data_dict[key]
            line = "$key: $value\n"
            data = data * line
        end
    end

    #println(template)
    #println(data)

    template = replace(template, "[DATA]" => data)

    return template
end

function format_csv_print(data_dict)
    iteration = data_dict["iteration"]
    line = "$iteration;"

    for key in sort(collect(keys(data_dict)))
        if key != "iteration"
            value = data_dict[key]
            line = line * "$value;"
        end
    end

    return line * "\n"    
end

function format_csv_header(data_dict)
    line = "iteration;"

    for key in sort(collect(keys(data_dict)))
        if key != "iteration"
            line = line * "$key;"
        end
    end

    return line * "\n"    
end

function save_to_csv(filename, data_dict)
    line = format_csv_print(data_dict)

    if !isfile(filename)
        header = format_csv_header(data_dict)
        line = header * line
    end

    open(filename, "a") do file
        write(file, line)
    end
end

function print_iteration_data(data_dict; cli_only_checkpoint = false, checkpoint = 0, to_csv = false, output_csv = "")
    iteration = data_dict["iteration"]

    if cli_only_checkpoint
        if mod(iteration, checkpoint) == 0 || data_dict["stop_condition"] != ""
            print(format_cli_print(data_dict))
        end
    else
        print(format_cli_print(data_dict))
    end

    # will save to csv at every call to print_iteration_data
    if to_csv
        save_to_csv(output_csv, data_dict)
    end
end

function check_int_param( param_value::String, min_value::Int64, max_value::Int64)
    try
        return ( min_value <= parse(Int64, param_value) <= max_value )
    catch
        return false
    end
    return true
end

function check_float_param( param_value::String, min_value::Float64, max_value::Float64)
    try
        return ( min_value <= parse(Float64, param_value) <= max_value )
    catch
        return false
    end
    return true
end

function check_opt_param(param_value::String, param_options::String)
    try
        return (param_value in split(param_options, "|"))
    catch
        return false
    end
    return true
end