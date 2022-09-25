using Graphs
using EzXML
using GraphsMatching
using GraphPlot

#include("utils.jl")
#include("logs.jl")
#include("commons.jl")

function calculate_graph_cost(graph, cost_matrix, n)
    graph_cost = 0
    for edge in edges(graph)
        u, v = src(edge), dst(edge)
        edge_cost = cost_matrix[u, v]
        # println("edge $u -> $v = ", edge_cost)

        graph_cost = graph_cost + edge_cost
    end

    return graph_cost
end

function count_vertex_edges(graph, n)
    result = zeros(1, n)

    for edge in edges(graph)
        u, v = src(edge), dst(edge)
        result[1, u] = result[1, u] + 1
        result[1, v] = result[1, v] + 1
    end

    return result
end

# Parameters
# xml_graph -> graph defined in XML format
# Returns
# c -> cost matrix
# n -> number of elements in the graph
function graph_to_cost_matrix(xml_graph)
    graph = root(xml_graph)
    n = countelements(graph)
    # println(n)
    c = zeros(Float64,n,n) # cost matrix
    i = 1
    for vertex in eachelement(graph)
        for edge in eachelement(vertex)
            # println(nodecontent(edge))
            # println(edge["cost"])
            j = parse(Int64, nodecontent(edge))+1
            # println(j)
            # println(edge["cost"])
            custo = parse(Float64,edge["cost"])
            # println(custo)
            c[i,j] = custo
        end
        i=i+1
    end

    return c, n
end

# Parameters
# n -> number of elements in the graph
# Returns 
# new_graph -> complete graph with n nodes
function create_complete_graph(n, is_digraph = false)
    new_graph = undef
    if is_digraph
        new_graph = DiGraph(n)
    else
        new_graph = SimpleGraph(n)
    end

	for i in 1:n
		for j in i+1:n
            if i == j
                continue
            end
            add_edge!(new_graph, i, j)
		end
	end

    # println("create_complete_graph ", new_graph)

    return new_graph
end

# Parameters
# graph -> graph 
# c -> cost matrix
# n -> number of vertices
# Returns
# mst -> the minimum spanning tree of the provided graph
# c -> cost matrix (same as the one provided)
function minimum_spanning_tree(graph, c, n) 
    # if nv(graph) > ne(graph) # necessary if "graph" complete?
    #     e = prim_mst(graph, c)
    # else
    e = kruskal_mst(graph, c)
    # end
    # println("A ", e)

    mst = SimpleGraph(n)
    for edge in e
        add_edge!(mst, src(edge), dst(edge))
    end

    # draw(PDF(string("./output/mst", ".pdf"), 16cm, 16cm), gplot(mst, nodelabel=1:nv(mst)))
    # graphplot(mst, nodelabel=1:nv(mst), curves=false)

    return mst, c

end

# Function created to decouple graph creation from 
# spanning tree creation.
# Parameters
# cost_matrix -> cost matrix of the graph
# n -> number of nodes in the graph
# Returns
# mst -> the minimum spanning tree
function minimum_spanning_tree_from_complete_graph(cost_matrix, n)
    graph = create_complete_graph(n)
    mst, c = minimum_spanning_tree(graph, cost_matrix, n)
    return mst, c
end

# g -> SimpleGraph object
# c -> cost matrix
function count_degree_and_form_subgraph(g, c) 

    # g, c = minimum_spanning_tree()
    deg = degree(g)
    # println("deg = ", deg)
    index = zeros(Int64, nv(g))
    subg = SimpleGraph()
    k = 1
    for d in 1:nv(g)
        if (deg[d] & 1) == 1
            index[k] = d
            add_vertex!(subg)
            k+=1
        end
    end
    # println("index = ", index) # indice dos vertices com grau negativo
    subgc = zeros(nv(subg),nv(subg))
    i = 1
    while index[i] != 0 && i <= nv(g)
        j=i+1
        if j >= nv(g)
            break
        end
        while j <= nv(g) && index[j] != 0
            add_edge!(subg, i, j)
            subgc[i,j] = c[index[i],index[j]]
            j+=1
        end
        i+=1
    end

    # graphplot(subg, nodelabel=1:nv(subg), curves=false)
    return subg, subgc, index
end

# subg -> subgraph (SimpleGraph object)
# subgc -> subgraph's cost matrix.
function mwpm(subg, subgc) 

    # draw(PDF(string("./output/subg", ".pdf"), 16cm, 16cm), gplot(subg, nodelabel=1:nv(subg)))
    if length(edges(subg)) <= 1
        return subg
    end
    w = Dict{Graphs.SimpleGraphs.SimpleEdge{Int64}, Float64}()
    for e in edges(subg)
        # println("mwpm ", src(e), " ", dst(e), " ", subgc[src(e),dst(e)], " ", typeof(subgc[src(e),dst(e)]))
        w[Edge(Int(src(e)),Int(dst(e)))] = subgc[src(e),dst(e)]
    end
    # println("w = ", w)
    # draw(PDF(string("./output/subg", ".pdf"), 16cm, 16cm), gplot(subg, nodelabel=1:nv(subg)))
    match = minimum_weight_perfect_matching(subg, w, tmaxscale=45.)
    temp = SimpleGraph(nv(subg))
    for i in 1:nv(subg)
        add_edge!(temp, i, match.mate[i])
    end
    return temp
end


function euler_path(g)
     
    vizitados = Array{Int64}(undef, 0)
    while length(vizitados) < nv(g)
        tempg = g
        u = 0
        for i in 1:nv(g)
            if (length(filter( x -> x == i, vizitados )) == 0)
                u = i
                break
            end
        end
        # println("u = ", u)
        v = 0
        k = Array{Int64}(undef, 0)
        j1 = neighbors(tempg, u)
        push!(k, u)
        while v != u
            # println("j1 = ", j1)
            for j in j1
                # println("j = ", j)
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
            # println("k = ", k)
            # break
        end
        for k1 in 1:length(k)
            if (length(filter( x -> x == k[k1], vizitados )) == 0)
                push!(vizitados, Int(k[k1]))
            end
        end
        # println("vizitados = ", vizitados)
    end
    vizitados = [vizitados;vizitados[1]]
    hamilton = SimpleGraph(nv(g))
    for ind in 1:length(vizitados)-1
        add_edge!(hamilton,  vizitados[ind], vizitados[ind+1])
    end
    return hamilton
end

function unite_and_hamilton(cost_matrix, num_vertices)
    g, c = minimum_spanning_tree_from_complete_graph(cost_matrix, num_vertices)
    # println("c = ", c)
    subg, subgc, index = count_degree_and_form_subgraph(g, c)
    mw = mwpm(subg, subgc)
    # g U mw
    for e in edges(mw)
        add_edge!(g, index[src(e)], index[dst(e)])
    end
    # deg = degree(g)
    hamilton = euler_path(g)
    return hamilton
end