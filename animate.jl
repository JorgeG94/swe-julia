using DelimitedFiles
using Plots

# Extract step number from "height_step_123.csv"
function extract_step_number(filename::AbstractString)
    m = match(r"height_step_(\d+)\.csv", filename)
    return isnothing(m) ? -1 : parse(Int, m.captures[1])
end

"""
    animate_height(files; outname="evolution.gif", fps=10)

Create an animation from CSV files with water heights.
Files are sorted by their numeric step number.
"""
function animate_height(files; outname="evolution.gif", fps=10)
    # Sort numerically by step index
    files = sort(files, by=extract_step_number)
    @info "Animating $(length(files)) frames → $outname"

    anim = @animate for filename in files
        step = extract_step_number(filename)
        H = readdlm(filename, ',')
        heatmap(H,
            c=:blues,
            xlabel="x", ylabel="y",
            title="Water Height (step = $step)",
            colorbar_title="h",
            aspect_ratio=:equal
        )
    end

    gif(anim, outname, fps=fps)
end

# --- Run it ---
using Glob
files = Glob.glob("output/height_step_*.csv", ".")
animate_height(files; outname="dam_break.gif", fps=10)

