using Graphs
using EzXML
using GraphsMatching
using GraphPlot

include("cristofides.jl")

function partial_graph(n, ignore)
    new_graph = SimpleGraph(n)

	for i in 1:n
		for j in i+1:n
            if i == j
                continue
            end

            if i == ignore || j == ignore
                continue
            end

            add_edge!(new_graph, i, j)
		end
	end

    println("partial_graph ", new_graph)

    return new_graph
end

# Parameters 
# cost_matrix -> cost_matrix
# num_nodes -> number of nodes in the graph
# Returns
# mst -> the minimum spanning tree turned into a 1-tree.
function one_tree_graph(cost_matrix, num_nodes)
    graph = partial_graph(num_nodes, 1)
    mst, c = minimum_spanning_tree(graph, cost_matrix, num_nodes)
    # Adding the two cheapest edges from node 1
    sorted_nodes = sortperm(cost_matrix[1, 2:end]) # must ignore the first column
    add_edge!(mst, 1, sorted_nodes[1] + 1) # + 1 to compensate for the indice of the first column that is ignored in the previous line 
    add_edge!(mst, 1, sorted_nodes[2] + 1)
    #println("bbbbbbb ", sorted_nodes)
    #println("cccccccc" , sorted_nodes[1], sorted_nodes[2])
    #println("aaaaaaaa ", cost_matrix[1, sorted_nodes[1]], cost_matrix[1, sorted_nodes[2]])
    #println(cost_matrix)
    #println(nv(mst))
    gplot(mst, nodelabel=1:nv(mst))

    return mst
end