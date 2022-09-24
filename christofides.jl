using Graphs
using EzXML
using GraphsMatching
using GraphPlot

include("utils.jl")
include("logs.jl")
include("commons.jl")

#
# Function minimum_spanning_tree_from_complete_graph 
# created to decouple graph creation from spanning tree creation. 
# 
# Parameters
# exp_id -> experiment identification
# step_id -> step identification
# cost_matrix -> cost matrix of the graph
# num_nodes -> number of nodes in the graph
# Returns
# mst -> the minimum spanning tree
# cost -> cost matrix used
function minimum_spanning_tree_from_complete_graph(exp_id, step_id, cost_matrix, num_nodes)

    @debug "minimum_spanning_tree_from_complete_graph..."

    save_step(exp_id,"minimum_spanning_tree_from_complete_graph:create_complete_graph","start","step") 
    graph = create_complete_graph(num_nodes)
    save_step(exp_id,"minimum_spanning_tree_from_complete_graph:create_complete_graph","finish","step") 

    save_step(exp_id,"minimum_spanning_tree_from_complete_graph:minimum_spanning_tree","start","step") 
    mst, cost = minimum_spanning_tree(exp_id, step_id, graph, cost_matrix, num_nodes)
    save_step(exp_id,"minimum_spanning_tree_from_complete_graph:minimum_spanning_tree","finish","step") 

    @debug "minimum_spanning_tree_from_complete_graph...done"

    return mst, cost
end


#
# Function unite_and_hamilton 
# base function to obtain hamilton cicle for TSP using cost_matrix 
# 
# Parameters
# exp_id -> experiment identification
# step_id -> step identification
# cost_matrix -> cost matrix of the graph
# num_vertices -> number of nodes in the graph
# Returns
# hamilton -> a hamilton circle

function unite_and_hamilton(exp_id, step_id, cost_matrix, num_vertices)

    @debug "unite_and_hamilton..."

    save_step(exp_id,"unite_and_hamilton:minimum_spanning_tree_from_complete_graph","start","step") 
    graph, cost = minimum_spanning_tree_from_complete_graph(exp_id, step_id, cost_matrix, num_vertices)
    save_step(exp_id,"unite_and_hamilton:minimum_spanning_tree_from_complete_graph","finish","step") 
    
    @debug cost

    save_step(exp_id,"unite_and_hamilton:count_degree_and_form_subgraph","start","step")    
    subg, subgc, index = count_degree_and_form_subgraph(exp_id, step_id, graph, cost)
    save_step(exp_id,"unite_and_hamilton:count_degree_and_form_subgraph","finish","step") 

    mw = min_weight_perfect_matching(exp_id, subg, subgc)
    # g U mw
    for e in edges(mw)
        add_edge!(graph, index[src(e)], index[dst(e)])
    end
    
    @debug degree(graph)

    save_step(exp_id,"unite_and_hamilton:euler_path","start","step") 
    hamilton = euler_path(exp_id, step_id, graph)
    save_step(exp_id,"unite_and_hamilton:euler_path","finish","step") 

    @debug "unite_and_hamilton...done"

    return hamilton
end