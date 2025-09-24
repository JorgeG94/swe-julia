include("grid.jl")
include("state.jl")
include("time_driver.jl")

using .Grid2D
using .State2D
using .TimeDriver

function print_simulation_setup(xmin, xmax, ymin, ymax, dx,
                                h_left, h_right, x_split, t_end, cfl)
    println("=== Simulation Setup ===")
    println(" Domain:")
    println("   x ∈ [$xmin, $xmax], y ∈ [$ymin, $ymax]")
    println("   dx = $dx, dy = $dx  (uniform grid)")

    println(" Initial condition (dam break):")
    println("   h_left  = $h_left")
    println("   h_right = $h_right")
    println("   x_split = $x_split")

    println(" Time integration:")
    println("   t_end = $t_end")
    println("   CFL   = $cfl")
    println("=========================")
end


function main()
    # Domain and discretization
    xmin, ymin = 0.0, 0.0
    xmax, ymax = 512.0, 512.0
    dx = 1.0

    # Simulation parameters
    h_left, h_right, x_split = 20.0, 7.0, 150.0
    t_end, cfl = 35.05, 0.45

    print_simulation_setup(xmin, xmax, ymin, ymax, dx,
                       h_left, h_right, x_split, t_end, cfl)

    # Grid and state
    grid = Grid(xmin, xmax, ymin, ymax, dx)
    state = State(grid)
    initialize_dam_break!(state, h_left, h_right, x_split)

    # Run solver
    evolve!(state, t_end, cfl)

    println("Done.")
end

main()

