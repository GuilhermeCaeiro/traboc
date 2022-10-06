using Graphs
using EzXML
using GraphsMatching
using GraphPlot
using Logging
using Compose
import Cairo
import Fontconfig

#
# Function calculate_graph_cost
# calculates graph total cost based on cost matrix
#  
# Parameters
# graph -> graph to be measured
# cost_matrix -> cost matrix used
# Returns 
# graph_cost -> graph cost
function calculate_graph_cost(graph, cost_matrix, n)
    graph_cost = 0

    @debug "calculate_graph_cost..."
    for edge in edges(graph)
        u, v = src(edge), dst(edge)
        edge_cost = cost_matrix[u, v]
        
        @debug "edge $u -> $v cost = " edge_cost

        graph_cost = graph_cost + edge_cost
    end

    @debug "calculate_graph_cost...done"

    return graph_cost
end

#
# Function count_vertex_edges
# calculates number of edges for each vertex in given graph
#  
# Parameters
# graph -> graph to be measured
# n -> number of vertices
# Returns 
# result -> vector containig vertices and each edge count
function count_vertex_edges(graph, n)
    result = zeros(1, n)

    @debug "count_vertex_edges..."

    for edge in edges(graph)
        u, v = src(edge), dst(edge)
        result[1, u] = result[1, u] + 1
        result[1, v] = result[1, v] + 1
    end

    @debug "count_vertex_edges...done"

    return result
end

#
# Function graph_to_cost_matrix
# obtains cost matrix based on graph structued
#  
# Parameters
# xml_graph -> graph defined in XML format
# Returns
# c -> cost matrix
# n -> number of elements in the graph
function graph_to_cost_matrix(xml_graph)

    @debug "graph_to_cost_matrix..."

    graph = root(xml_graph)
    n = countelements(graph)

    @debug "countelements: " n

    c = zeros(n,n) # cost matrix
    i = 1
    for vertex in eachelement(graph)
        for edge in eachelement(vertex)

            @debug "nodecontent(edge): " nodecontent(edge)
            @debug "edge cost: " edge["cost"]

            j = parse(Int64, nodecontent(edge))+1

            @debug "edge index: " j

            cost = parse(Float64,edge["cost"])
            
            @debug "edge cost as float64: " cost

            c[i,j] = cost
        end
        i=i+1
    end

    @debug "graph_to_cost_matrix...done"

    return c, n
end

#
# Function create_complete_graph
# creates a complete graph with n vertices
#  
# Parameters
# n -> number of elements in the graph
# Returns 
# new_graph -> complete graph with n nodes
function create_complete_graph(n)
    @debug "create_complete_graph..."

    new_graph = SimpleGraph(n)

	for i in 1:n
		for j in i+1:n
            if i == j
                continue
            end
            add_edge!(new_graph, i, j)
		end
	end

    @debug begin
        plt = gplot(new_graph, nodelabel=1:nv(new_graph))
        draw(PNG(string("./temp/complete_graph",get_next_id(),".png"), 16cm, 16cm), plt)
        display(plt)
    end

    @debug new_graph
    @debug "create_complete_graph...done"

    return new_graph
end

#
# Function minimum_spanning_tree
# returns a minimum spnning tree version of a given graph and cost matrix
# using Prim's or Kruskal's algorithms
#  
# Parameters
# exp_id -> experiment identification
# step_id -> step identification
# graph -> graph 
# c -> cost matrix
# n -> number of vertices
# Returns
# mst -> the minimum spanning tree of the provided graph
# c -> cost matrix (same as the one provided)
function minimum_spanning_tree(exp_id, step_id, graph, c, n) 
    graph_name = "mst"

    @debug "minimum_spanning_tree..."

    if nv(graph) > ne(graph) # necessary if "graph" complete?
        e = prim_mst(graph, c)
        graph_name = string(step_id,"_prim_",graph_name)
        save_result(exp_id,step_id,"mst","prim")
    else
        e = kruskal_mst(graph, c)
        graph_name = string(step_id,"_kruskal_",graph_name)
        save_result(exp_id,step_id,"mst","kruskal")
    end

    @debug "edges_vector: "
    @debug e

    mst = SimpleGraph(n)
    for edge in e
        add_edge!(mst, src(edge), dst(edge))
    end
    
    @debug begin
        plt = gplot(mst, nodelabel=1:nv(mst))
        draw(PNG(string("./work/",exp_id,"/",graph_name,".png"), 16cm, 16cm), plt)    
        display(plt)
    end

    @debug "minimum_spanning_tree...done"

    return mst, c

end

#
# Function count_degree_and_form_subgraph
# 
# Parameters
# exp_id -> experiment identification
# step_id -> step identification
# graph -> SimpleGraph object
# cost_matrix -> cost matrix
# Returns
# subg
# subgc
# index

function count_degree_and_form_subgraph(exp_id, step_id, graph, cost_matrix) 

    @debug "count_degree_and_form_subgraph..."

    # g, c = minimum_spanning_tree()
    deg = degree(graph)

    @debug "degree: " deg

    index = zeros(Int64, nv(graph))
    subg = SimpleGraph()
    k = 1
    for d in 1:nv(graph)
        if (deg[d] & 1) == 1
            index[k] = d
            add_vertex!(subg)
            k+=1
        end
    end

    # indices dos vertices com grau negativo
    @debug "negative indexes: " index

    subgc = zeros(nv(subg),nv(subg))
    i = 1
    while i <= nv(graph) && index[i] != 0
        j=i+1
        while j <= nv(graph) && index[j] != 0
            add_edge!(subg, i, j)
            subgc[i,j] = cost_matrix[index[i],index[j]]
            j+=1
        end
        i+=1
    end

    @debug begin
        plt = gplot(subg, nodelabel=1:nv(subg)) 
        draw(PNG(string("./work/",exp_id,"/",step_id,"_subgraph.png"), 16cm, 16cm), plt)    
        display(plt)
    end

    @debug "count_degree_and_form_subgraph...done"

    return subg, subgc, index

end

#
# Function min_weight_perfect_matching 
# obtain a perfect matching with mininum weigth for given graph based on cost_matrix 
# 
# Parameters
# exp_id -> experiment_id
# graph -> graph from which to extract perfect matching
# cost_matrix -> cost matrix of the graph
# num_vertices -> number of nodes in the graph
# Returns
# match_graph -> graph with match vertices
function min_weight_perfect_matching(exp_id, step_id, graph, cost_matrix) 

    @debug "min_weight_perfect_matching..."

    if length(edges(graph)) <= 1
        return graph
    end

    #prepare weigths parameter
    w = Dict{Graphs.SimpleGraphs.SimpleEdge{Int64}, Float64}()
    for e in edges(graph)
        w[Edge(Int(src(e)),Int(dst(e)))] = cost_matrix[src(e),dst(e)]
    end

    match = minimum_weight_perfect_matching(graph, w)

    @debug begin
        for i in 1:nv(graph)
            @debug match.mate[i]
        end
    end

    match_graph = SimpleGraph(nv(graph))
    for i in 1:nv(graph)
        add_edge!(match_graph, i, match.mate[i])
    end

    @debug begin
        plt = gplot(match_graph, nodelabel=1:nv(match_graph))
        draw(PNG(string("./work/",exp_id,"/",step_id,"_match_graph.png"), 16cm, 16cm), plt)    
        display(plt)
    end

    @debug "min_weight_perfect_matching...done"

    return match_graph
end

#
# Function euler_path 
# return the euler path of a graph 
# 
# Parameters
# exp_id -> experiment identification
# step_id -> step identification
# graph -> graph from which to extract perfect matching
# Returns
# path -> euler path
function euler_path(exp_id, step_id, graph)
   
    @debug "euler_path..."

    visited = Array{Int64}(undef, 0)
    while length(visited) < nv(graph)
        tempg = graph
        u = 0
        for i in 1:nv(graph)
            if (length(filter( x -> x == i, visited )) == 0)
                u = i
                break
            end
        end

        @debug u

        v = 0
        k = Array{Int64}(undef, 0)
        j1 = neighbors(tempg, u)
        push!(k, u)
        while v != u

            @debug j1

            for j in j1

                @debug j

                v = j
                rem_edge!(tempg, k[length(k)], v)
                push!(k, v)
                break
            end
            j1 = neighbors(tempg, v)
            if length(j1) == 0
                j1 = neighbors(tempg, k[length(k)-1])
                if length(j1) == 0
                    break
                end
            end

            @debug k
            # break
        end
        for k1 in 1:length(k)
            if (length(filter( x -> x == k[k1], visited )) == 0)
                push!(visited, Int(k[k1]))
            end
        end

        @debug visited

    end

    visited = [visited;visited[1]]
    path = SimpleGraph(nv(graph))
    for ind in 1:length(visited)-1
        add_edge!(path,  visited[ind], visited[ind+1])
    end
   
    @debug path

    @debug "euler_path...done"

    return path
end

# Function to check if a graph consitutes a hamiltonian cycle.
# It assumes that the graph provided is symatric and connected.
# 
# Parameters
# graph -> a symetric SimpleGraph object
# Returns
# true if a hamiltonian cycle, or false otherwise
function check_hamiltonian_cycle(graph)
    nodes_checked = 0
    number_of_nodes = nv(graph)
    number_of_edges = ne(graph)
    previous_node = 0
    starting_node = 1
    next_node = starting_node

    # if there are only 2 nodes, there is no hamiltonian cycle
    if number_of_nodes <= 2
        #println("Two or less nodes.")
        return false
    end

    # attempts to follow a cycle where nodes_checked shout be equal number_of_nodes
    while true
        node = next_node
        node_degree = degree(graph, node)
        # all nodes must have degree == 2
        if node_degree != 2
            #println("Degree other than 2: node $node, degree $node_degree")
            return false
        end

        neighbor_nodes = neighbors(graph, node)
        next_node = deleteat!(neighbor_nodes, findall(x->x==previous_node, neighbor_nodes))[1] # probably inneficient, but there are only two values

        nodes_checked += 1
        previous_node = node

        if next_node == starting_node
            break
        end
    end

    #println(number_of_nodes, " ", nodes_checked)

    if number_of_nodes == nodes_checked
        return true
    else
        return false
    end
end