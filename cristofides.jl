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

    e = prim_mst(mst, c)
    println(e)

    temp = SimpleGraph(n)
    for edge in e
        add_edge!(temp, src(edge), dst(edge))
    end

    gplot(temp, nodelabel=1:nv(temp))

    return temp, c

end

function count_degree_and_form_subgraph(g, c)

    # g, c = minimum_spanning_tree()
    deg = zeros(Int64, nv(g))
    for e in edges(g)
        deg[dst(e)]+=1
        deg[src(e)]+=1
    end
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
    println(deg)
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

    # glot(subg, nodelabel=1:nv(subg), curves=false)
    return subg, subgc, index
end

function mwpm(subg, subgc)

    w = Dict{Graphs.SimpleGraphs.SimpleEdge{Int64}, Int64}()
    for e in edges(subg)
        w[Edge(Int(src(e)),Int(dst(e)))] = subgc[src(e),dst(e)]
    end
    match = minimum_weight_perfect_matching(subg, w)
    for i in 1:nv(subg)
        println(match.mate[i])
    end
    temp = SimpleGraph(nv(subg))
    for i in 1:nv(subg)
        add_edge!(temp, i, match.mate[i])
    end

    # gplot(temp, nodelabel=1:nv(subg), curves=false)
    return temp
end


function unite()
    g, c = minimum_spanning_tree()
    subg, subgc, index = count_degree_and_form_subgraph(g, c)
    mw = mwpm(subg, subgc)
    # g U mw
    for e in edges(mw)
        add_edge!(g, index[src(e)], index[dst(e)])
    end
    gplot(g, nodelabel=1:nv(g))

end