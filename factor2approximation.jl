using Graphs
using EzXML
using GraphsMatching
using GraphPlot

include("commons.jl")
include("logs.jl")
include("utils.jl")

# Function to produce a feasible solution with a factor 2 approximation
#
# Parameters
# cost_matrix -> cost matrix
# n -> number of vertices in the graph
# Returns
# ham_cycle_from_euler -> hamiltonian cycle (graph) resulting from the factor 2 approximation
function factor_2_approximation(exp_id, step_id, cost_matrix, n)
    graph = create_complete_graph(n)
    mst, c = minimum_spanning_tree(exp_id, step_id, graph, cost_matrix, n) 
    # Converting the minimum spanning tree to a digraph.
    # The duplicated edges are generated during the conversion.
    dgraph = SimpleDiGraph(mst)
    ham_cycle_from_euler = euler_path(exp_id, step_id, dgraph)
    return ham_cycle_from_euler
end