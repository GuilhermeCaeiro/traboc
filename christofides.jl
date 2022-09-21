using Graphs
using EzXML
using GraphsMatching
using GraphPlot
using Logging
using Compose
import Cairo
import Fontconfig

include("utils.jl")
include("logs.jl")
include("commons.jl")

#
# Function minimum_spanning_tree_from_complete_graph 
# created to decouple graph creation from spanning tree creation. 
# 
# Parameters
# exp_id -> experiment_id
# cost_matrix -> cost matrix of the graph
# num_nodes -> number of nodes in the graph
# Returns
# mst -> the minimum spanning tree
# cost -> cost matrix used
function minimum_spanning_tree_from_complete_graph(exp_id,cost_matrix, num_nodes)
    graph = create_complete_graph(num_nodes)
    mst, cost = minimum_spanning_tree(exp_id, graph, cost_matrix, num_nodes)
    return mst, cost
end


#
# Function unite_and_hamilton 
# base function to obtain hamilton cicle for TSP using cost_matrix 
# 
# Parameters
# exp_id -> experiment_id
# cost_matrix -> cost matrix of the graph
# num_vertices -> number of nodes in the graph
# Returns
# hamilton -> a hamilton circle
function unite_and_hamilton(exp_id,cost_matrix, num_vertices)

    @debug "unite_and_hamilton..."

    g, c = minimum_spanning_tree_from_complete_graph(exp_id, cost_matrix, num_vertices)
    
    subg, subgc, index = count_degree_and_form_subgraph(exp_id, g, c)

    mw = min_weight_perfect_matching(exp_id, subg, subgc)
    # g U mw
    for e in edges(mw)
        add_edge!(g, index[src(e)], index[dst(e)])
    end

    plt = gplot(g, nodelabel=1:nv(g), nodelabeldist=1.5, nodelabelangleoffset=π/4)
    draw(PNG(string("./work/",exp_id,"/gUmw_graph.png"), 16cm, 16cm), plt)
    @debug begin
        display(plt)
    end

    deg = degree(g)
    @debug "degree: " deg

    hamilton = SimpleGraph(nv(g))
    df = dfs_tree(g, articulation(g)[1])

    plt = gplot(df, nodelabel=1:nv(df))
    draw(PNG(string("./work/",exp_id,"/dfs_tree.png"), 16cm, 16cm), plt)
    @debug begin
        display(plt)
    end

    @debug "neighbors: " neighbors(df, articulation(g)[1])

    @debug "processing articulation..."

    for a in articulation(g)
        df = dfs_tree(g, a)
        nei = neighbors(df, a)

        @debug "nei = " nei

        allN = neighbors(g, a)

        @debug "allN = " allN

        for n in 1:2

            @debug "nei[" n "] = " nei[n]

            add_edge!(hamilton, a, nei[n])
            deleteat!(allN, findall( x -> x == nei[n], allN ))
        end

        @debug "allN = " allN

        for ind in 1:length(allN)-1
            add_edge!(hamilton,  allN[ind], allN[ind+1])
        end
    end

    @debug "processing articulation...done"
    @debug "processing edges..."

    for e in edges(g)
        i = src(e)
        j = dst(e)
        if deg[i] == 2 && deg[j] == 2
            add_edge!(hamilton, i, j)
        end
        @debug "edge: " i " -> " j
    end

    @debug "processing edges...done"

    plt = gplot(hamilton, nodelabel=1:nv(g))
    draw(PNG(string("./work/",exp_id,"/hamilton_circle.png"), 16cm, 16cm), plt)
    @debug begin
        display(plt)
    end

    @debug "unite_and_hamilton...done"

    return hamilton
end
