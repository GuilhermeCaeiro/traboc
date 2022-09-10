using Graphs
using EzXML
using GraphsMatching
using GraphPlot

testdata = String(read("test_cristofides.xml"));
doc = parsexml(testdata)

function minimum_spanning_tree() 

    graph = root(doc)
    n = countelements(graph)
    # println(n)
    c = zeros(n,n)
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

    mst = SimpleGraph(n)

	for i in 1:n
		for j in i+1:n
            if i == j
                continue
            end
            add_edge!(mst, i, j)
		end
	end

    if nv(mst) > ne(mst)
        e = prim_mst(mst, c)
    else
        e = kruskal_mst(mst, c)
    end
    println(e)

    temp = SimpleGraph(n)
    for edge in e
        add_edge!(temp, src(edge), dst(edge))
    end
    
    # graphplot(temp, nodelabel=1:nv(temp), curves=false)

    return temp, c

end

function count_degree_and_form_subgraph(g, c) 

    # g, c = minimum_spanning_tree()
    deg = degree(g)
    println(deg)
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
    println(index) # indice dos vertices com grau negativo
    subgc = zeros(nv(subg),nv(subg))
    i = 1
    while index[i] != 0 && i <= nv(g)
        j=i+1
        while index[j] != 0 && j <= nv(g)
            add_edge!(subg, i, j)
            subgc[i,j] = c[index[i],index[j]]
            j+=1
        end
        i+=1
    end

    # graphplot(subg, nodelabel=1:nv(subg), curves=false)
    return subg, subgc, index
end

function mwpm(subg, subgc) 

    w = Dict{Graphs.SimpleGraphs.SimpleEdge{Int64}, Int64}()
    for e in edges(subg)
        w[Edge(Int(src(e)),Int(dst(e)))] = subgc[src(e),dst(e)]
    end
    match = minimum_weight_perfect_matching(subg, w)
    # for i in 1:nv(subg)
    #     println(match.mate[i])
    # end
    temp = SimpleGraph(nv(subg))
    for i in 1:nv(subg)
        add_edge!(temp, i, match.mate[i])
    end

    # graphplot(temp, nodelabel=1:nv(subg), curves=false)
    return temp
end


function unite_and_hamilton()
    g, c = minimum_spanning_tree()
    subg, subgc, index = count_degree_and_form_subgraph(g, c)
    mw = mwpm(subg, subgc)
    # g U mw
    for e in edges(mw)
        add_edge!(g, index[src(e)], index[dst(e)])
    end

    # graphplot(g, nodelabel=1:nv(g), curves=false, nodelabeldist=1.5, nodelabelangleoffset=π/4)
    deg = degree(g)
#     println(deg)
    hamilton = SimpleGraph(nv(g))
    df = dfs_tree(g, articulation(g)[1])
#     println(neighbors(df, articulation(g)[1]))
    for a in articulation(g)
        df = dfs_tree(g, a)
        nei = neighbors(df, a)
#         println("nei = ", nei)
        allN = neighbors(g, a)
#         println("allN = ", allN)
        for n in 1:2
#             println(nei[n])
            add_edge!(hamilton, a, nei[n])
            deleteat!(allN, findall( x -> x == nei[n], allN ))
        end
#         println("allN = ", allN)
        for ind in 1:length(allN)-1
            add_edge!(hamilton,  allN[ind], allN[ind+1])
        end
    end
    for e in edges(g)
        i = src(e)
        j = dst(e)
        if deg[i] == 2 && deg[j] == 2
            add_edge!(hamilton, i, j)
        end
        # println(i, " -> ", j)
    end
    graphplot(hamilton, nodelabel=1:nv(g), curves=false)
end
