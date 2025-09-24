# grid.jl

module Grid2D

export Grid, generate_2d_grids, get_grid_coords

"""
    Grid2D(xmin, xmax, ymin, ymax, dx; dy=dx)

A 2D grid spanning `[xmin, xmax] × [ymin, ymax]` with spacing `dx` and optionally `dy`.
"""
struct Grid
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    dx::Float64
    dy::Float64
    nx::Int
    ny::Int
    x::Vector{Float64}
    y::Vector{Float64}
    x2d::Matrix{Float64}
    y2d::Matrix{Float64}
end

# Constructor that builds 1D and 2D coordinates
function Grid(xmin::Real, xmax::Real, ymin::Real, ymax::Real, dx::Real; dy::Real=dx)
    dx = float(dx)
    dy = float(dy)

    nx = Int(round((xmax - xmin)/dx)) + 1
    ny = Int(round((ymax - ymin)/dy)) + 1

    x = collect(xmin:dx:xmax)
    y = collect(ymin:dy:ymax)

    x2d = repeat(x, 1, ny)   
    y2d = repeat(y', nx, 1) 

    return Grid(xmin, xmax, ymin, ymax, dx, dy, nx, ny, x, y, x2d, y2d)
end

"""
    generate_2d_grids(grid::Grid2D)

Return `(x2d, y2d)` meshgrids from a Grid2D object.
"""
generate_2d_grids(grid::Grid) = (grid.x2d, grid.y2d)

"""
    get_grid_coords(grid::Grid2D, i, j)

Return the `(x, y)` coordinates at indices `(i, j)`.
"""
get_grid_coords(grid::Grid, i::Int, j::Int) = (grid.x[i], grid.y[j])

end # module

