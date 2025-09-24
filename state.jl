# state.jl

module State2D

using ..Grid2D: Grid   # import our Grid type

export State, initialize_dam_break!

"""
    State(grid::Grid)

Hold shallow-water state variables (water height, x/y momentum, ground elevation).
"""
mutable struct State
    grid::Grid
    water_height::Matrix{Float64}
    x_momentum::Matrix{Float64}
    y_momentum::Matrix{Float64}
    ground_elevation::Matrix{Float64}
end

function State(grid::Grid)
    nx, ny = grid.nx, grid.ny
    h  = zeros(nx, ny)
    hu = zeros(nx, ny)
    hv = zeros(nx, ny)
    zb = zeros(nx, ny)
    return State(grid, h, hu, hv, zb)
end

"""
    initialize_dam_break!(state::State, h_left, h_right, x_split)

Set up a dam-break initial condition: water height is `h_left` for `x < x_split`,
and `h_right` otherwise. Momenta are set to zero.
"""
function initialize_dam_break!(state::State, h_left::Real, h_right::Real, x_split::Real)
    g = state.grid
    for j in 1:g.ny, i in 1:g.nx
        state.water_height[i,j] = g.x[i] < x_split ? h_left : h_right
        state.x_momentum[i,j]   = 0.0
        state.y_momentum[i,j]   = 0.0
        state.ground_elevation[i,j] = 0.0
    end
    return state
end

end # module

