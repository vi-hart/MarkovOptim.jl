function LightGraphs.add_edge!(sys::MarkovSystem{<:RandomParameter}, e::SimpleGraphEdge)
    graph = get_graph(sys)
    add_edge!(graph, e) || return false
    e_set = Set(Tuple(e))
    edges_set = Set.(Tuple.((edges(graph))))
    index = findfirst(x->isequal(x, e_set), edges_set)
    @inbounds insert!(edge_vals(sys), index)
    return true
end
function LightGraphs.add_edge!(sys::MarkovSystem{RandomParameter}, x)
    graph = get_graph(sys)
    return add_edge!(sys, edgetype(graph)(x))
end
LightGraphs.add_edge!(sys::MarkovSystem{RandomParameter}, x, y) = add_edge!(M, edgetype(get_graph(M))(x, y))

function LightGraphs.rem_edge!(sys::MarkovSystem{<:RandomParameter}, e::SimpleGraphEdge)
    graph = get_graph(sys)
    e_set = Set(Tuple(e))
    edges_set = Set.(Tuple.((edges(graph))))
    index = findfirst(x->isequal(x, e_set), edges_set)
    rem_edge!(graph, e) || return false
    @inbounds popat!(edge_vals(sys), index)
    return true
end
function LightGraphs.rem_edge!(sys::MarkovSystem{<:RandomParameter}, u::Integer, v::Integer)
    return rem_edge!(M, edgetype(get_graph(M))(u, v))
end

function LightGraphs.add_vertex!(sys::MarkovSystem{<:RandomParameter})
	LightGraphs.add_vertex!(get_graph(sys)) || return false
    push!(node_vals(sys))
    return true
end

function LightGraphs.rem_vertex!(sys::MarkovSystem{<:RandomParameter}, v::Integer)
    g = get_graph(sys)
    v in vertices(g) || return false
    n = nv(g)
    self_loop_n = false  # true if n is self-looped (see #820)

    # remove the in_edges from v
    srcs = copy(inneighbors(g, v))
    @inbounds for s in srcs
        rem_edge!(sys, edgetype(g)(s, v))
    end
    # remove the in_edges from the last vertex
    neigs = copy(inneighbors(g, n))
    @inbounds for s in neigs
        rem_edge!(sys, edgetype(g)(s, n))
    end
    if v != n
        # add the edges from n back to v
        @inbounds for s in neigs
            if s != n  # don't add an edge to the last vertex - see #820.
                add_edge!(sys, edgetype(g)(s, v))
            else
                self_loop_n = true
            end
        end
    end
    if self_loop_n
        add_edge!(sys, edgetype(g)(v, v))
    end
    pop!(g.fadjlist)
    popat!(node_vals(sys), v)
    return true
end

LightGraphs.rem_vertex!(sys::MarkovSystem{<:RandomParameter}) = rem_vertex!(sys, nv(get_graph(sys)))

function connect!(sys::MarkovSystem{<:RandomParameter})
    # Connect disconnected clusters
	clusters = connected_components(get_graph(sys))
	if size(clusters,1) > 1
		v1 = rand(clusters[1])
		v2 = rand(clusters[2])
		add_edge!(sys, v1, v2)
		connect!(sys)
	end
end
