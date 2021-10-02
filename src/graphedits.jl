function LightGraphs.add_edge!(M::MarkovModel{T, P}, e::SimpleGraphEdge{T}) where {T, P}
    graph = get_graph(M)
    add_edge!(graph, e) || return false
    e_set = Set(Tuple(e))
    edges_set = Set.(Tuple.((edges(graph))))
    index = findfirst(x->isequal(x, e_set), edges_set)
    @inbounds insert!(edge_parameters(M), index)
    return true
end
function LightGraphs.add_edge!(M::AbstractMarkovModel, x)
    graph = get_graph(M)
    return add_edge!(M, edgetype(graph)(x))
end
LightGraphs.add_edge!(M::AbstractMarkovModel, x, y) = add_edge!(M, edgetype(get_graph(M))(x, y))

function LightGraphs.rem_edge!(M::MarkovModel{T, P}, e::SimpleGraphEdge{T}) where {T, P}
    graph = get_graph(M)
    e_set = Set(Tuple(e))
    edges_set = Set.(Tuple.((edges(graph))))
    index = findfirst(x->isequal(x, e_set), edges_set)
    rem_edge!(graph, e) || return false
    @inbounds popat!(edge_parameters(M), index)
    return true
end
function LightGraphs.rem_edge!(M::MarkovModel{T, P}, u::Integer, v::Integer) where {T, P}
    return rem_edge!(M, edgetype(get_graph(M))(T(u), T(v)))
end

function LightGraphs.add_vertex!(M::MarkovModel)
    LightGraphs.add_vertex!(get_graph(M)) || return false
    push!(node_parameters(M))
    return true
end

function LightGraphs.rem_vertex!(M::MarkovModel, v::Integer)
    g = get_graph(M)
    v in vertices(g) || return false
    n = nv(g)
    self_loop_n = false  # true if n is self-looped (see #820)

    # remove the in_edges from v
    srcs = copy(inneighbors(g, v))
    @inbounds for s in srcs
        rem_edge!(M, edgetype(g)(s, v))
    end
    # remove the in_edges from the last vertex
    neigs = copy(inneighbors(g, n))
    @inbounds for s in neigs
        rem_edge!(M, edgetype(g)(s, n))
    end
    if v != n
        # add the edges from n back to v
        @inbounds for s in neigs
            if s != n  # don't add an edge to the last vertex - see #820.
                add_edge!(M, edgetype(g)(s, v))
            else
                self_loop_n = true
            end
        end
    end
    if self_loop_n
        add_edge!(M, edgetype(g)(v, v))
    end
    pop!(g.fadjlist)
    popat!(node_parameters(M), v)
    return true
end

function connect!(M::MarkovModel)
    # Connect disconnected clusters
	clusters = connected_components(get_graph(M))
	if size(clusters,1) > 1
		v1 = rand(clusters[1])
		v2 = rand(clusters[2])
		add_edge!(M, v1, v2)
		connect!(M)
	end
end
