function _incidence(graph::SimpleGraph)
    N = nv(graph)
    E = ne(graph)
    ic = zeros(Int, N, E)
    for i in 1:N
        for (j,edge) in enumerate(edges(graph))
            if edge.src == i
                ic[i,j] = -1
            elseif edge.dst == i
                ic[i,j] = 1
            end
        end
    end
    return ic
end

function _directed_incidence(graph::SimpleGraph)
    N = nv(graph)
    E = ne(graph)
    ic = zeros(Int, N, 2E)
    for i in 1:N
        for (j,edge) in enumerate(edges(graph))
            if edge.src == i
                ic[i,2j-1] = -1
                ic[i,2j] = 1
            elseif edge.dst == i
                ic[i,2j-1] = 1
                ic[i,2j] = -1
            end
        end
    end
    return ic
end

_vecOfvec(x::M) where M<:AbstractMatrix = [x[1,:] for i in 1:size(x,1)]
_vecOfvec(x::P) where P<:RandomParameter = getfield(x, :ps)
