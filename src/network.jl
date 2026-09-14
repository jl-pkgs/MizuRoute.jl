"""Directed acyclic river network.

`downstream[i]` is the internal downstream reach index; `0` means outlet.
`upstream[i]` stores immediate upstream internal indices. `order` is a
headwater-to-outlet topological order.
"""
struct RiverNetwork
    reach_id::Vector{Int}
    downstream_id::Vector{Int}
    downstream::Vector{Int}
    upstream::Vector{Vector{Int}}
    order::Vector{Int}
    id_to_index::Dict{Int,Int}
end

function RiverNetwork(reach_id::AbstractVector{<:Integer}, downstream_id::AbstractVector{<:Integer})
    n = length(reach_id)
    n > 0 || throw(ArgumentError("river network cannot be empty"))
    length(downstream_id) == n || throw(ArgumentError("reach_id and downstream_id lengths differ"))
    rid = Int.(reach_id)
    dsid = Int.(downstream_id)
    length(unique(rid)) == n || throw(ArgumentError("reach IDs must be unique"))
    index = Dict(id => i for (i, id) in pairs(rid))

    downstream = zeros(Int, n)
    upstream = [Int[] for _ in 1:n]
    for i in 1:n
        d = dsid[i]
        if d <= 0
            downstream[i] = 0
        else
            haskey(index, d) || throw(ArgumentError("downstream reach ID $d for reach $(rid[i]) is not present"))
            j = index[d]
            j == i && throw(ArgumentError("reach $(rid[i]) drains to itself"))
            downstream[i] = j
            push!(upstream[j], i)
        end
    end

    indegree = length.(upstream)
    queue = Int[i for i in 1:n if indegree[i] == 0]
    order = Int[]
    first = 1
    while first <= length(queue)
        i = queue[first]
        first += 1
        push!(order, i)
        j = downstream[i]
        if j != 0
            indegree[j] -= 1
            indegree[j] == 0 && push!(queue, j)
        end
    end
    length(order) == n || throw(ArgumentError("river network contains a directed cycle"))

    return RiverNetwork(rid, dsid, downstream, upstream, order, index)
end

Base.length(net::RiverNetwork) = length(net.reach_id)
upstream_indices(net::RiverNetwork, i::Integer) = net.upstream[i]
downstream_index(net::RiverNetwork, i::Integer) = net.downstream[i]
headwaters(net::RiverNetwork) = findall(isempty, net.upstream)
outlets(net::RiverNetwork) = findall(==(0), net.downstream)
