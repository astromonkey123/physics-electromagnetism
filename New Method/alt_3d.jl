"""
    3D Electric Field

* Definitions
    Imports, defining constants, and setting up the figure and axes for plotting.
    Defining structs for point charges and the electric field

* Point Charges
    Instances of the PointCharge struct with customizable position functions

* Electric Field
    Generation of the electric field and its divergence and curl at a given time

* Plotting
    Animation and plotting of the generated electric field and point charges

Ths script uses the Jefimenko–Feynman formula to calculate the electric and magnetic
fields for a configuration of moving point charges.

astromonkey123 4/2025
"""

# ------------------------------ Definitions ------------------------------ #

using GLMakie
using CalculusWithJulia
using BenchmarkTools
using Base.Threads
using Metal

# Physics parameters
const ϵ₀ = 8.854187817e-12 # Permittivity of space
const μ₀ = 1.25663706127e-6 # Permeability of space
const c = 299792458 # Speed of light

# Simulation parameters
const dt = 7.5e-10 # Time step 
const steps = 300 # Number of steps
max = 10
min = -10
step = 2
const xs = [i for i in min:step:max]; const ys = [i for i in min:step:max]; const zs = [i for i in min:step:max] # Spatial bounds


# Point charge object with position and charge
mutable struct pointCharge
    x::Function # Position as a function of time
    q::Float64 # Charge
end

# ------------------------------ Useful Math Functions ------------------------------ #

const Δ = 1e-5

function ddt(x, t)
    return ( x(t + Δ) - x(t) ) / Δ
end

function ddt2(x, t)
    return ( ddt(x, t + Δ) - ddt(x, t) ) / Δ
end

# ------------------------------ Electromagnetic Field ------------------------------ #

# Generate the electric and magnetic fields at a given time
function generate_fields(t₀, charge_list)
    E::Vector{Vector{Float64}} = []
    B::Vector{Vector{Float64}} = []

    # Iterate through all positions
    for x in xs, y in ys, z in zs
        E_vector = zeros(3)
        B_vector = zeros(3)

        # Iterate through all point charges (pc)
        for pc in charge_list
            k = -pc.q/(4*π*ϵ₀)
            r(t) = hypot( ([x, y, z] - pc.x(t) )...) # Distance from location to the pc
            delay(t) = t - ( r(t) / c )
            r_delay(t) = hypot(( [x, y, z] - pc.x(delay(t)) )...) # Distance from location to the pc's previous position
            e_r(t) = ( [x, y, z] - pc.x(delay(t)) ) / r_delay(t)

            t1 = e_r(t₀) / ( r_delay(t₀)^2 )

            t2 = (r_delay(t₀) / c) * ddt(t -> ( e_r(t) / r_delay(t)^2 ), t₀)

            t3 = (1 / ( c^2 )) * ddt2(e_r, t₀)

            E_vector += k * (t1 + t2 + t3)
        end

        # Iterate through all point charges (pc)
        for pc in charge_list
            r(t) = hypot( ([x, y, z] - pc.x(t) )...) # Distance from location to the pc
            delay(t) = t - ( r(t) / c )
            r_delay(t) = hypot(( [x, y, z] - pc.x(delay(t)) )...) # Distance from location to the pc's previous position
            e_r(t) = ( [x, y, z] - pc.x(delay(t)) ) / r_delay(t)

            B_vector += cross(-e_r(t₀), E_vector/c)
        end

        append!(E, [E_vector])
        append!(B, [B_vector])
    end

    E_reshaped = reshape(E, length(xs), length(ys), length(zs))
    B_reshaped = reshape(B, length(xs), length(ys), length(zs))

    return E_reshaped, B_reshaped
end

# ------------------------------ Plotting ------------------------------ #

function animate()
    # Graphing parameters
    set_theme!(theme_black())
    fig = Figure(size=(650, 650))
    ax = Axis3(fig[1, 1], limits=(xs[begin], xs[end], ys[begin], ys[end], zs[begin], zs[end])) # Axis with the spatial bounds
    points = [Point3f(x, y, z) for x in xs, y in ys, z in zs] # Points for the vectors

    # Point charge instances and charge list
    charge_list = [pointCharge(t -> [5*cos(0.15*c*t), 5*sin(0.15*c*t), 0], 1.0)]

    for t in 0:dt:(steps * dt)
        E, B = generate_fields(t, charge_list) # Generate the electromagnetic field

        # Colors for E
        E_lengths = [hypot(vector...) for vector in E]
        avg_E_length = median(E_lengths) / 1
        E_lengths = [length / avg_E_length for length in E_lengths]
        E_colors = [RGBAf(l/100, 0, 0, l - 50) for l in E_lengths]
        E_display = [Vec3f(v/hypot(v...)) for v in E]

        # Colors for B
        B_lengths = [hypot(vector...) for vector in B]
        avg_B_length = median(B_lengths) / 20
        B_lengths = [length / avg_B_length for length in B_lengths]
        B_colors = [RGBAf(0, l/500, l/100, l - 50) for l in B_lengths]
        B_display = [Vec3f(v/hypot(v...)) for v in B]

        empty!(ax) # Clear the axes each animation cycle
        
        # Plot the electric field and curl
        arrows!(ax, vec(permutedims(points, [3, 2, 1])), vec(E_display), lengthscale=0.05, arrowsize=0.01, color=vec(E_colors))

        # Plot the magnetic field and curl
        arrows!(ax, vec(permutedims(points, [3, 2, 1])), vec(B_display), lengthscale=0.05, arrowsize=0.01, color=vec(B_colors))

        # Plot the point charges with colors based on their charge
        for pc in charge_list
            scatter!(ax, pc.x(t)..., color=RGBf(0, 0.5-pc.q, 0.5+pc.q))
        end

        display(fig)
        sleep(0.001)
    end
end

animate()