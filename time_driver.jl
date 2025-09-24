# time_driver.jl

module TimeDriver

const gravity = 9.81

using ..State2D: State
using Printf


export evolve!, Flux, zero_flux!

using DelimitedFiles

"""
    write_water_height_to_csv(state, step)

Write the water height field to a CSV file named "height_step_<step>.csv".
"""
function write_water_height_to_csv(state::State, step::Int; outdir="output")
    mkpath(outdir)  
    filename = joinpath(outdir, "height_step_$(step).csv")
    writedlm(filename, state.water_height, ',')
end

elapsed(t0) = (time_ns() - t0) * 1e-9

# --------------------------
# Flux container
# --------------------------

mutable struct Flux
    flux_h::Matrix{Float64}
    flux_hu::Matrix{Float64}
    flux_hv::Matrix{Float64}
end

function Flux(nx::Int, ny::Int)
    return Flux(zeros(nx, ny), zeros(nx, ny), zeros(nx, ny))
end

function zero_flux!(flux::Flux)
    fill!(flux.flux_h,  0.0)
    fill!(flux.flux_hu, 0.0)
    fill!(flux.flux_hv, 0.0)
    return flux
end

# --------------------------
# Main evolve loop
# --------------------------

"""
    evolve!(state::State, t_end, cfl)

Advance the shallow water state until `t_end` using CFL condition `cfl`.
"""
function evolve!(state::State, t_end::Real, cfl::Real)
    t     = 0.0
    step  = 0
    dt    = 0.0
    print_interval = 10

    initial_mass = sum(state.water_height) * state.grid.dx * state.grid.dy
    
    println("Grid points: ", state.grid.nx * state.grid.ny)

    #println(rpad("Step",8), rpad("Time",12), rpad("dt",12),
    #        rpad("Mass",12), rpad("ΔMass",12), rpad("s/step",12))
    @printf("%-8s %-12s %-12s %-12s %-12s %-12s\n",
        "Step", "Time", "dt", "Mass", "ΔMass", "s/step")

    # Allocate flux containers
    flux_x = Flux(state.grid.nx+1, state.grid.ny)
    flux_y = Flux(state.grid.nx,   state.grid.ny+1)

    outdir = "output"
    
    # Remove old output directory if it exists
    if isdir(outdir)
        rm(outdir; recursive=true, force=true)
    end
    
    # Recreate clean directory
    mkpath(outdir)
    # Start overall timer
    t0_total = time_ns()

    t0_step = 0
    while t < t_end
        # Start step timer

        if step % print_interval == 0
        t0_step = time_ns()
        end

        dt = compute_dt(state, cfl)
        if t + dt > t_end
            dt = t_end - t
        end

        apply_reflective_boundaries!(state)
        compute_rusanov_fluxes_xy!(state, flux_x, flux_y)

        before_mass = sum(state.water_height) * state.grid.dx * state.grid.dy

        update_state!(state, flux_x, flux_y, dt)
        enforce_min_height!(state, 1.0e-5)

        after_mass = sum(state.water_height) * state.grid.dx * state.grid.dy

        t    += dt
        step += 1

        if step % print_interval == 0
            # Elapsed time per step interval
            write_water_height_to_csv(state, step; outdir="output")
            sec_per_step = elapsed(t0_step) / print_interval

            total_mass = sum(state.water_height) * state.grid.dx * state.grid.dy
            Δmass = after_mass - before_mass


            @printf("%-8d %-12.4f %-12.4f %-12.4e %-12.4e %-12.6f\n",
        step, t, dt, total_mass, Δmass, sec_per_step)
            #println(rpad(string(step),8),
            #        rpad(string(round(t,digits=4)),12),
            #        rpad(string(round(dt,digits=4)),12),
            #        rpad(string(round(total_mass,digits=4)),12),
            #        rpad(string(round(Δmass,digits=4)),12),
            #        rpad(string(round(sec_per_step,digits=6)),12))

        end
    end

    # Total timing
    total_time = elapsed(t0_total)

    final_mass = sum(state.water_height) * state.grid.dx * state.grid.dy
    println("Final mass: ", final_mass,
            ", Initial mass: ", initial_mass,
            ", Δmass: ", final_mass - initial_mass,
            ", Lost ", 100*(1 - final_mass/initial_mass), " %")

    println("Total evolve! time = $(round(total_time,digits=4)) s",
            ", Average time/step = $(round(total_time/step,digits=6)) s")

    return state
end


"""
    compute_dt(state::State, cfl)

Compute the timestep based on the CFL condition:
    dt = cfl * min(dx,dy) / max_wave_speed
"""
function compute_dt(state::State, cfl::Real)
    nx, ny = state.grid.nx, state.grid.ny
    dx, dy = state.grid.dx, state.grid.dy

    max_speed = 0.0

    @inbounds for j in 1:ny, i in 1:nx
        h = state.water_height[i,j]
        if h > eps()   # avoid division by 0
            u = state.x_momentum[i,j] / h
            v = state.y_momentum[i,j] / h
            a = max(abs(u) + sqrt(gravity*h),
                    abs(v) + sqrt(gravity*h))
            if a > max_speed
                max_speed = a
            end
        end
    end

    if max_speed > 0.0
        return cfl * min(dx, dy) / max_speed
    else
        return 1.0e-3
    end
end

"""
    apply_reflective_boundaries!(state)

Apply reflective boundary conditions to all edges of the domain.
"""
function apply_reflective_boundaries!(state::State)
    nx, ny = state.grid.nx, state.grid.ny

    # Corners
    state.water_height[1,1]   = state.water_height[2,2]
    state.x_momentum[1,1]     = -state.x_momentum[2,2]
    state.y_momentum[1,1]     = -state.y_momentum[2,2]

    state.water_height[1,ny]  = state.water_height[2,ny-1]
    state.x_momentum[1,ny]    = -state.x_momentum[2,ny-1]
    state.y_momentum[1,ny]    = -state.y_momentum[2,ny-1]

    state.water_height[nx,1]  = state.water_height[nx-1,2]
    state.x_momentum[nx,1]    = -state.x_momentum[nx-1,2]
    state.y_momentum[nx,1]    = -state.y_momentum[nx-1,2]

    state.water_height[nx,ny] = state.water_height[nx-1,ny-1]
    state.x_momentum[nx,ny]   = -state.x_momentum[nx-1,ny-1]
    state.y_momentum[nx,ny]   = -state.y_momentum[nx-1,ny-1]

    # Left/right edges (excluding corners)
    for j in 2:ny-1
        state.water_height[1,j]  = state.water_height[2,j]
        state.x_momentum[1,j]    = -state.x_momentum[2,j]
        state.y_momentum[1,j]    =  state.y_momentum[2,j]

        state.water_height[nx,j] = state.water_height[nx-1,j]
        state.x_momentum[nx,j]   = -state.x_momentum[nx-1,j]
        state.y_momentum[nx,j]   =  state.y_momentum[nx-1,j]
    end

    # Bottom/top edges (excluding corners)
    for i in 2:nx-1
        state.water_height[i,1]  = state.water_height[i,2]
        state.x_momentum[i,1]    =  state.x_momentum[i,2]
        state.y_momentum[i,1]    = -state.y_momentum[i,2]

        state.water_height[i,ny] = state.water_height[i,ny-1]
        state.x_momentum[i,ny]   =  state.x_momentum[i,ny-1]
        state.y_momentum[i,ny]   = -state.y_momentum[i,ny-1]
    end

    return state
end

"""
    update_state!(state, flux_x, flux_y, dt)

Update water height and momentum using the flux divergence.
"""
function update_state!(state::State, flux_x, flux_y, dt::Real)
    nx, ny = state.grid.nx, state.grid.ny
    dx, dy = state.grid.dx, state.grid.dy

    @inbounds for j in 1:ny, i in 1:nx
        dh  = -dt/dx * (flux_x.flux_h[i+1,j]  - flux_x.flux_h[i,j]) +
              -dt/dy * (flux_y.flux_h[i,j+1]  - flux_y.flux_h[i,j])

        dhu = -dt/dx * (flux_x.flux_hu[i+1,j] - flux_x.flux_hu[i,j]) +
              -dt/dy * (flux_y.flux_hu[i,j+1] - flux_y.flux_hu[i,j])

        dhv = -dt/dx * (flux_x.flux_hv[i+1,j] - flux_x.flux_hv[i,j]) +
              -dt/dy * (flux_y.flux_hv[i,j+1] - flux_y.flux_hv[i,j])

        state.water_height[i,j] += dh
        state.x_momentum[i,j]   += dhu
        state.y_momentum[i,j]   += dhv
    end
    return state
end


"""
    enforce_min_height!(state, h_min)

Enforce non-negative water height by zeroing small values and momenta.
"""
function enforce_min_height!(state::State, h_min::Real)
    nx, ny = state.grid.nx, state.grid.ny
    @inbounds for j in 1:ny, i in 1:nx
        if state.water_height[i,j] < h_min
            state.water_height[i,j] = 0.0
            state.x_momentum[i,j]   = 0.0
            state.y_momentum[i,j]   = 0.0
        end
    end
    return state
end

function compute_rusanov_fluxes_xy!(state::State, flux_x::Flux, flux_y::Flux)
    nx, ny = state.grid.nx, state.grid.ny
    dx, dy = state.grid.dx, state.grid.dy

    sqroot_gravity = sqrt(gravity)
    half_gravity   = 0.5 * gravity

    @inbounds for j in 1:ny, i in 1:nx
        # -------------------
        # X-direction fluxes
        # -------------------
        hL  = state.water_height[i, j]
        huL = state.x_momentum[i, j]
        hvL = state.y_momentum[i, j]

        hR  = state.water_height[i+1, j]
        huR = state.x_momentum[i+1, j]
        hvR = state.y_momentum[i+1, j]

        if hL < eps(); hL = eps(); end
        if hR < eps(); hR = eps(); end

        uL = hL > eps() ? huL / hL : 0.0
        vL = hL > eps() ? hvL / hL : 0.0
        uR = hR > eps() ? huR / hR : 0.0
        vR = hR > eps() ? hvR / hR : 0.0

        fluxL1 = huL
        fluxL2 = huL * uL + half_gravity * hL^2
        fluxL3 = huL * vL

        fluxR1 = huR
        fluxR2 = huR * uR + half_gravity * hR^2
        fluxR3 = huR * vR

        cL = abs(uL) + sqroot_gravity * sqrt(hL)
        cR = abs(uR) + sqroot_gravity * sqrt(hR)
        amax = max(cL, cR)

        flux_x.flux_h[i+1, j]  = 0.5*(fluxL1 + fluxR1) - 0.5*amax*(hR  - hL)
        flux_x.flux_hu[i+1, j] = 0.5*(fluxL2 + fluxR2) - 0.5*amax*(huR - huL)
        flux_x.flux_hv[i+1, j] = 0.5*(fluxL3 + fluxR3) - 0.5*amax*(hvR - hvL)

        # -------------------
        # Y-direction fluxes
        # -------------------
        hL  = state.water_height[i, j]
        huL = state.x_momentum[i, j]
        hvL = state.y_momentum[i, j]

        hR  = state.water_height[i, j+1]
        huR = state.x_momentum[i, j+1]
        hvR = state.y_momentum[i, j+1]

        if hL < eps(); hL = eps(); end
        if hR < eps(); hR = eps(); end

        uL = hL > eps() ? huL / hL : 0.0
        vL = hL > eps() ? hvL / hL : 0.0
        uR = hR > eps() ? huR / hR : 0.0
        vR = hR > eps() ? hvR / hR : 0.0

        fluxL1 = huL
        fluxL2 = huL * uL + half_gravity * hL^2
        fluxL3 = huL * vL

        fluxR1 = huR
        fluxR2 = huR * uR + half_gravity * hR^2
        fluxR3 = huR * vR

        cL = abs(vL) + sqroot_gravity * sqrt(hL)
        cR = abs(vR) + sqroot_gravity * sqrt(hR)
        amax = max(cL, cR)

        flux_y.flux_h[i, j+1]  = 0.5*(fluxL1 + fluxR1) - 0.5*amax*(hR  - hL)
        flux_y.flux_hu[i, j+1] = 0.5*(fluxL2 + fluxR2) - 0.5*amax*(huR - huL)
        flux_y.flux_hv[i, j+1] = 0.5*(fluxL3 + fluxR3) - 0.5*amax*(hvR - hvL)
    end

    # Zero out boundaries
    @inbounds for i in 1:size(flux_y.flux_h,1)
        flux_y.flux_h[i,1]     = 0.0
        flux_y.flux_h[i,ny+1]  = 0.0
        flux_y.flux_hu[i,1]    = 0.0
        flux_y.flux_hu[i,ny+1] = 0.0
        flux_y.flux_hv[i,1]    = 0.0
        flux_y.flux_hv[i,ny+1] = 0.0
    end

    @inbounds for j in 1:size(flux_x.flux_h,2)
        flux_x.flux_h[1,j]     = 0.0
        flux_x.flux_h[nx+1,j]  = 0.0
        flux_x.flux_hu[1,j]    = 0.0
        flux_x.flux_hu[nx+1,j] = 0.0
        flux_x.flux_hv[1,j]    = 0.0
        flux_x.flux_hv[nx+1,j] = 0.0
    end

    return nothing
end



end # module

