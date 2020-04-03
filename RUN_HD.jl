# Basic test script - generates the file sample_hd_traj.csv containing one hobby drone trajectory

# Include the files
include("HobbyDroneInterface.jl")
# Set the random seed
Random.seed!(1)
# Generate the trajectory file
generate_trajectory_file(1.0, "output/hd_traj.csv")