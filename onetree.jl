using Graphs
using EzXML
using GraphsMatching
using GraphPlot

include("commons.jl")
include("logs.jl")
include("utils.jl")

#
# Function partial_graph 
# returns a partial graph from number of nodes and igore nodes i and j 
#
# Parameters
# num_vertices -> number of nodes in the graph
# ignore -> nodes to be ignored
# Returns
# new_graph -> a partial graph

function partial_graph(num_vertices, ignore)

    @debug "partial_graph..."

    new_graph = SimpleGraph(num_vertices)

	for i in 1:num_vertices
		for j in i+1:num_vertices
            if i == j
                continue
            end

            if i == ignore || j == ignore
                continue
            end

            add_edge!(new_graph, i, j)
		end
	end

    @debug new_graph

    @debug "partial_graph...done"

    return new_graph
end

#
# Function one_tree_graph 
# return a one-tree based on give number of nodes and cost matrix 
# 
# Parameters 
# exp_id -> experiment identification
# step_id -> step identification
# cost_matrix -> cost_matrix
# num_nodes -> number of nodes in the graph
# Returns
# mst -> the minimum spanning tree turned into a 1-tree.
function one_tree_graph(exp_id, step_id, cost_matrix, num_nodes)

    @debug "one_tree_graph..."

    graph = partial_graph(num_nodes, 1)

    save_step(exp_id,"one_tree_graph:minimum_spanning_tree","start","step") 
    mst, c = minimum_spanning_tree(exp_id, step_id, graph, cost_matrix, num_nodes)
    save_step(exp_id,"one_tree_graph:minimum_spanning_tree","finish","step") 

    # Adding the two cheapest edges from node 1
    sorted_nodes = sortperm(cost_matrix[1, 2:end]) # must ignore the first column
    add_edge!(mst, 1, sorted_nodes[1] + 1) # + 1 to compensate for the indice of the first column that is ignored in the previous line 
    add_edge!(mst, 1, sorted_nodes[2] + 1)

    #println("bbbbbbb ", sorted_nodes)
    #println("cccccccc" , sorted_nodes[1], sorted_nodes[2])
    #println("aaaaaaaa ", cost_matrix[1, sorted_nodes[1]], cost_matrix[1, sorted_nodes[2]])
    #println(cost_matrix)
    #println(nv(mst))

    @debug begin
        plt = gplot(mst, nodelabel=1:nv(mst))
        draw(PNG(string("./work/",exp_id,"/",step_id,"_one_tree.png"), 16cm, 16cm), plt)
        display(plt)
    end

    @debug "one_tree_graph...done"

    return mst
end