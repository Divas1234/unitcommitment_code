using DelimitedFiles
using Plots

# Define the path to the reformed_results file
data_path = "/Users/yuanyiping/Documents/GitHub/unit_commitment_code/master-3/out/trick_draw_3dsurface/.data/reformed_result.txt"

# Read the data
data_matrix = readdlm(data_path, Float64)

# Get dimensions
rows, cols = size(data_matrix)
x = 1:cols
y = 1:rows

# Create figure
fig = Figure(resolution = (800, 600))
ax = Axis3(fig[1, 1],
    xlabel = "Time Period",
    ylabel = "Generator Units",
    zlabel = "Power Output",
    title = "Unit Commitment Surface Plot",
)

# Create surface plot
surface(data_matrix,
    # colormap = :viridis,
    # shading = true,
    # transparency = false,
)

# Add colorbar
# Colorbar(fig[1, 2], limits = (minimum(data_matrix), maximum(data_matrix)),
#     colormap = :viridis,
#     label = "Power Level")

display(fig)

# Save the figure
save("/Users/yuanyiping/Documents/GitHub/unit_commitment_code/master-3/out/trick_draw_3dsurface/.fig/unit_commitment_surface_plot.png", fig)
