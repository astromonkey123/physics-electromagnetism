"""
    2D Electric Field

* Definitions
    Imports, defining constants, and setting up the figure and axes for plotting.
    Defining structs for point charges and the electric field

* Point Charges
    Instances of the PointCharge struct with customizable position functions

* Electric Field
    Generation of the electric field and its divergence and curl at a given time

* Plotting
    Animation and plotting of the generated electric field and point charges

astromonkey123 12/2024
"""

# ------------------------------ Definitions ------------------------------ #

using GLMakie
using CalculusWithJulia

# Physics parameters
const c = 299792458 # Speed of light

# Simulation parameters
dt = 5e-10 # Time step
steps = 250 # Number of steps
xs = [i for i in -10:0.5:10]; ys = [i for i in -10:0.5:10] # Spatial bounds

# Graphing parameters
set_theme!(theme_black())
fig = Figure(size=(650, 650))
ax = Axis(fig[1, 1], limits=(xs[begin], xs[end], ys[begin], ys[end]), 
    xgridvisible=false, ygridvisible=false, xtickcolor=:black, 
    xticklabelsvisible=false, ytickcolor=:black, yticklabelsvisible=false) # Simplified axis with the spatial bounds
points = [Point2f(x, y) for x in xs for y in ys] # Points for the vectors

# Point charge object with position and charge
mutable struct PointCharge
    xt::Function # Position as a function of time
    x::Observable{Float64}
    y::Observable{Float64}
    q::Float64 # Charge
end

# Electric field object containing the electric field with its divergence and curl at a specific time
Base.@kwdef mutable struct ElectricField
    time::Float64 # Timestamp
    field::Observable{Matrix} = Observable([zeros(2) for x in xs, y in ys]) # Electric field
    div::Observable{Matrix} = Observable([0.0 for x in xs, y in ys]) # Divergence
    curl::Observable{Matrix} = Observable([0.0 for x in xs, y in ys]) # Curl
end

# ------------------------------ Point Charges ------------------------------ #

# Position functions
function x1(t)
    t < 0 ? t = 0 : t = t
    [2*cos(0.25*c*t), 2*sin(0.25*c*t)]
end 

function x2(t)
    t < 0 ? t = 0 : t = t
    [2*cos(0.25*c*t + π), 2*sin(0.25*c*t + π)]
end 

# Point charge instances
pc_1 = PointCharge(x1, Observable(x1(0)[1]), Observable(x1(0)[2]), 1.0)
pc_2 = PointCharge(x2, Observable(x2(0)[1]), Observable(x2(0)[2]), 1.0)

charge_list = [pc_1, pc_2]

# ------------------------------ Electric Field ------------------------------ #

# Generate the electric field throughout space at a given time
function generate_electric_field(t::Float64, E_field::ElectricField)
    vector_field::Vector{Vector{Float64}} = []
    div_field::Vector{Float64} = []
    curl_field::Vector{Float64} = []

    # Iterate through all positions
    for x in xs, y in ys
        vector = zeros(3)
        div = 0.0
        curl = zeros(3)

        # Iterate through all point charges (pc)
        for pc in charge_list
            # Find the contribution to the electric field from pc at a location
            # Must be (x, y, z) for the curl operator to work, but it's always evaluated at z = 0
            function electric_field(x, y, z)
                r = hypot(x - pc.xt(t)[1], y - pc.xt(t)[2]) # Distance from location to the pc
                delay = r / c # Delay from the speed of light
                r_delay = hypot(x - pc.xt(t - delay)[1], y - pc.xt(t - delay)[2]) # Distance from location to the pc's previous position
    
                # Electric field strength from that pc by Coulomb's Law
                u = pc.q * (1/r_delay^2) * (x - pc.xt(t - delay)[1]) # X component
                v = pc.q * (1/r_delay^2) * (y - pc.xt(t - delay)[2]) # Y component
    
                return [u, v, 0]
            end
            E(v) = electric_field(v...)
    
            # Add the contributions to the field and its divergence and curl from that pc
            vector += electric_field(x, y, 0)
            div += (∇⋅E)(x, y, 0)
            curl += (∇×E)(x, y, 0)
        end

        # Append the field and its divergence and curl to the previously calculated locations
        append!(vector_field, [vector[1:2]])
        append!(div_field, div)
        append!(curl_field, [curl[3]])
    end

    E_lengths = [sum(abs.(vector)) for vector in vector_field]

    for i in eachindex(vector_field)
        field_vector = vector_field[i]
        size = hypot(field_vector[1], field_vector[2])
        vector_field[i] = field_vector * 1/size
    end

    # Reshape the electric field and curl into 2D matrices and return as an ElectricField instance
    E_field.time = t
    E_field.field[] = reshape(vector_field, length(xs), length(ys))
    E_field.div[] = transpose(reshape(div_field, length(xs), length(ys)))
    E_field.curl[] = transpose(reshape(curl_field, length(xs), length(ys)))

    return E_lengths
end

# ------------------------------ Plotting - ----------------------------- #

E_colors = Observable([RGBf(0, 0, 0) for _ in eachindex(xs) for _ in eachindex(ys)])

# Intial electric field state (all zeros)
E = ElectricField(time=0)

# Plot divergence and curl
divergence = heatmap!(ax, xs, ys, E.div, colormap=:tofino, colorrange=(-1, 1), alpha=1, interpolate=true)
# curl = heatmap!(ax, xs, ys, E.curl, colormap=:tofino, colorrange=(-1, 1), alpha=1, interpolate=true)

# Plot the electric field
vectors = arrows!(ax, points, E.field, lengthscale=0.4, arrowsize=5, color=E_colors)

# Plot the point charges with colors based on their charge
for pc in charge_list
    scatter!(ax, pc.x, pc.y, color=RGBf(0, 0.5-pc.q, 0.5+pc.q))
end

# Plot colorbars
Colorbar(fig[2, 1], colormap=:tofino, colorrange=(-1, 1), vertical=false, label="Divergence (∇⋅E)")
# Colorbar(fig[2, 1], colormap=:tofino, colorrange=(-1, 1), vertical=false, label="Curl (∇×E)")

record(fig, "2d_anim.mp4", 0:dt:(steps*dt)) do t
    E_lengths = generate_electric_field(t, E) # Generate the electric field

    # Calculate vector colors for plotting
    E_colors[] = [RGBf(2l, 2l, 2l) for l in E_lengths]

    # Update the positions of the point charges
    for pc in charge_list
        pc.x[] = pc.xt(t)[1]
        pc.y[] = pc.xt(t)[2]
    end
end