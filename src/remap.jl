"""Convert runoff depth rate [m/s] and area [m²] to volume rate [m³/s]."""
runoff_depth_to_volume(runoff_mps::Real, area_m2::Real) = Float64(runoff_mps) * Float64(area_m2)

"""Map gridded/HRU runoff depth [m/s] to reach lateral inflow [m³/s].

`cell_to_reach[k]` is the internal reach index receiving cell `k`; use `0` to
ignore a cell. Optional `weights[k]` can represent fractional overlap or any
precomputed conservative remapping fraction.
"""
function map_runoff!(qlat::AbstractVector, runoff_depth::AbstractVector,
    cell_area::AbstractVector, cell_to_reach::AbstractVector{<:Integer}; weights=nothing)
    ncell = length(runoff_depth)
    length(cell_area) == ncell || throw(ArgumentError("cell_area length mismatch"))
    length(cell_to_reach) == ncell || throw(ArgumentError("cell_to_reach length mismatch"))
    weights === nothing || length(weights) == ncell || throw(ArgumentError("weights length mismatch"))
    fill!(qlat, zero(eltype(qlat)))
    for k in 1:ncell
        r = Int(cell_to_reach[k])
        r == 0 && continue
        1 <= r <= length(qlat) || throw(BoundsError(qlat, r))
        w = weights === nothing ? 1.0 : Float64(weights[k])
        qlat[r] += Float64(runoff_depth[k]) * Float64(cell_area[k]) * w
    end
    return qlat
end

function map_runoff(nreach::Integer, runoff_depth::AbstractVector,
    cell_area::AbstractVector, cell_to_reach::AbstractVector{<:Integer}; weights=nothing)
    qlat = zeros(Float64, nreach)
    map_runoff!(qlat, runoff_depth, cell_area, cell_to_reach; weights=weights)
end
