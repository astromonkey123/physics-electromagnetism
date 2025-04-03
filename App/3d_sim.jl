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

astromonkey123 12/2024
"""

# ------------------------------ Definitions ------------------------------ #

using GLMakie
using CalculusWithJulia
using BenchmarkTools

# Physics parameters
const ϵ₀ = 8.854187817e-12 # Permittivity of space
const μ₀ = 1.25663706127e-6 # Permeability of space
const c = 299792458 # Speed of light

# Simulation parameters
const dt = 7.5e-10 # Time step 
const steps = 250 # Number of steps
const xs = [i for i in -5:0.625:5]; const ys = [i for i in -5:0.625:5]; const zs = [i for i in -5:0.625:5] # Spatial bounds

# Point charge object with position and charge
mutable struct PointCharge
    x::Function # Position as a function of time
    q::Float64 # Charge
end

# Vector field object containing the field with its divergence and curl at a specific time
Base.@kwdef mutable struct VectorField
    time::Float64 # Timestamp
    field::Array = [zeros(3) for _ in xs, _ in ys, _ in zs] # Vector Field
    div::Array = [0.0 for _ in xs, _ in ys, z in zs] # Divergence
    curl::Array = [zeros(3) for _ in xs, _ in ys, _ in zs] # Curl
end


# ------------------------------ Point Charges ------------------------------ #

# Position functions
function x1(t)
    t < 0 ? t = 0 : t = t
    [3*cos(0.25*c*t), 3*sin(0.25*c*t), 0] # Donut
    # [0, 0, 3*sin(0.25*c*t)] # Hourglass
    # [3*cos(0.05*c*t), 2*sin(0.3*c*t), 0] # Snake
    # R = 3; a = 20; s = 1/(steps * dt); [sqrt(R^2 - (R*s*t)^2) * cos(a*s*t), sqrt(R^2 - (R*s*t)^2) * sin(a*s*t), R*s*t] # Half-sphere
end

# ------------------------------ Electric Field ------------------------------ #

# Generate the electric field throughout space at a given time
function generate_electric_field(t, Eₙ₋₁, Bₙ₋₁, charge_list)
    vector_field::Vector{Vector{Float64}} = []
    div_field::Vector{Float64} = []
    curl_field::Vector{Vector{Float64}} = []

    # Ampere's Law
    # ∇×B = μ₀ ( J + ϵ₀ ( ∂E / ∂t ) )
    # ∂E = ∂t( ( ( ∇×B / μ₀ ) - J ) / ϵ₀ )

    # ∂E = dt * (((Bₙ₋₁.curl / μ₀) - 0) / ϵ₀)
    # ∂E = Bₙ₋₁.curl / μ₀

    # Iterate through all positions
    for x in xs, y in ys, z in zs
        vector = zeros(3)
        div = 0.0
        curl = zeros(3)

        # Iterate through all point charges (pc)
        for pc in charge_list
            # Find the contribution to the electric field from pc at a location
            # Must be (x, y, z) for the curl operator to work, but it's always evaluated at z = 0
            function electric_field(x, y, z)
                # TODO: Think about using hypot((pos - pc.x(t))...) to consolidate for r and r_delay
                r = hypot(x - pc.x(t)[1], y - pc.x(t)[2], z - pc.x(t)[3]) # Distance from location to the pc
                delay = r / c # Delay from the speed of light
                r_delay = hypot(x - pc.x(t - delay)[1], y - pc.x(t - delay)[2], z - pc.x(t - delay)[3]) # Distance from location to the pc's previous position
    
                # Electric field strength from that pc by Coulomb's Law 
                r = (1/4*π*ϵ₀) * pc.q * (1/r_delay^2)
                û = (x - pc.x(t - delay)[1]) / r_delay
                v̂ = (y - pc.x(t - delay)[2]) / r_delay
                ŵ = (z - pc.x(t - delay)[3]) / r_delay

                # TODO: Ampere's Law
    
                return r * [û, v̂, ŵ]
            end
            E(v) = electric_field(v...)
    
            # Add the contributions to the field and its divergence and curl from that pc
            vector += electric_field(x, y, z)
            div += (∇⋅E)(x, y, z)
            curl += (∇×E)(x, y, z)

        end

        # Append the field and its divergence and curl to the previously calculated locations
        append!(vector_field, [vector])
        append!(div_field, div)
        append!(curl_field, [curl])
    end

    field = reshape(vector_field, length(xs), length(ys), length(zs))
    div = div_field
    curl = reshape(curl_field, length(xs), length(ys), length(zs))

    # Reshape the electric field and curl into 2D matrices and return as an ElectricField instance
    return VectorField(t, field, div, curl)
end

# ------------------------------ Magnetic Field ------------------------------ #

# Generate the magnetic field throughout space at a given time
function generate_magnetic_field(t, Eₙ₋₁, Bₙ₋₁, charge_list)
    vector_field::Vector{Vector{Float64}} = []
    div_field::Vector{Float64} = []
    curl_field::Vector{Vector{Float64}} = []

    # Faraday's Law
    # ∇×E = -( ∂B / ∂t )
    # ∂B = -∂t( ∇×E )

    ∂B = -dt * Eₙ₋₁.curl
    magnetic_field_vec = vec(Bₙ₋₁.field)
    ∂B_vec = vec(∂B)

    # Iterate through all positions
    for x in xs, y in ys, z in zs
        vector = zeros(3)
        div = 0.0
        curl = zeros(3)

        # Find the magnetic field at a location
        function magnetic_field(x, y, z)
            x_step = xs[2] - xs[1]
            y_step = ys[2] - ys[1]
            z_step = zs[2] - zs[1]

            idx = trunc(Int, (length(ys) * length(zs) * ((x - xs[1]) / x_step)) + (length(zs) * ((y - ys[1]) / y_step)) + ((z - zs[1]) / z_step) + 1)

            return magnetic_field_vec[idx] + ∂B_vec[idx]
        end
        B(v) = magnetic_field(v...)

        # Add the contributions to the field and its divergence and curl from that pc
        vector = magnetic_field(x, y, z)
        div = (∇⋅B)(x, y, z)
        curl = (∇×B)(x, y, z)

        # Append the field and its divergence and curl to the previously calculated locations
        append!(vector_field, [vector])
        append!(div_field, div)
        append!(curl_field, [curl])
    end

    field = reshape(vector_field, length(xs), length(ys), length(zs))
    div = div_field
    curl = reshape(curl_field, length(xs), length(ys), length(zs))

    # TODO: Helmholtz decomposition (make B divergence-free)

    # Reshape the magnetic field and curl into 2D matrices and return as an VectorField instance
    return VectorField(t, field, div, curl)
end

# ------------------------------ Plotting ------------------------------ #

function animate()
    # Graphing parameters
    set_theme!(theme_black())
    fig = Figure(size=(650, 650))
    ax = Axis3(fig[1, 1], limits=(xs[begin], xs[end], ys[begin], ys[end], zs[begin], zs[end])) # Axis with the spatial bounds
    points = [Point3f(x, y, z) for x in xs, y in ys, z in zs] # Points for the vectors

    # Point charge instances
    pc_1 = PointCharge(x1, 1.0)

    charge_list = [pc_1]

    # Electric and magnetic fields
    Eₙ = VectorField(time=0)
    Bₙ = VectorField(time=0)
    Eₙ₋₁ = VectorField(time=0)
    Bₙ₋₁ = VectorField(time=0)

    for t in 0:dt:(steps * dt)
        # Store the previous electric and magnetic fields
        Eₙ₋₁ = Eₙ
        Bₙ₋₁ = Bₙ

        # Generate the new fields based on the ones from the previous time step
        Eₙ = generate_electric_field(t, Eₙ₋₁, Bₙ₋₁, charge_list) # Generate the electric field
        Bₙ = generate_magnetic_field(t, Eₙ₋₁, Bₙ₋₁, charge_list) # Generate the magnetic field

        # Colors for E
        E_lengths = [hypot(vector...) for vector in Eₙ.field]
        avg_E_length = median(E_lengths)
        avg_E_length = 2e-13 # TODO: Magic number
        E_lengths = [length / avg_E_length for length in E_lengths]
        E_colors = [RGBAf(l/50, 0, 0, l - 15) for l in E_lengths]
        Eₙ_field_display = [Vec3f(vector/hypot(vector...)) for vector in Eₙ.field]

        # Colors for ∇×E
        # E_curl_lengths = [hypot(vector...) for vector in Eₙ.curl]
        # E_curl_colors = [RGBAf(0, l * 1e12, 0, l * 5e11) for l in E_curl_lengths]
        # Eₙ_curl_display = [Vec3f(vector) for vector in Eₙ.curl]

        # Colors for B
        B_lengths = [hypot(vector...) for vector in Bₙ.field]
        avg_B_length = median(B_lengths)
        avg_B_length = 3e-22 # TODO: Magic number
        B_lengths = [length / avg_B_length for length in B_lengths]
        B_colors = [RGBAf(0, l/500, l/100, l - 50) for l in B_lengths]
        Bₙ_field_display = [Vec3f(vector/hypot(vector...)) for vector in Bₙ.field]

        # Colors for ∇×B
        # B_curl_lengths = [hypot(vector...) for vector in Bₙ.curl]
        # B_curl_colors = [RGBAf(0, l * 1e12, 0, l * 5e11) for l in B_curl_lengths]
        # Bₙ_curl_display = [Vec3f(vector) for vector in Bₙ.curl]

        empty!(ax) # Clear the axes each animation cycle
        
        # Plot the electric field and curl
        arrows!(ax, vec(permutedims(points, [3, 2, 1])), vec(Eₙ_field_display), lengthscale=0.1, arrowsize=0.02, color=vec(E_colors))
        # arrows!(ax, vec(permutedims(points, [3, 2, 1])), vec(Eₙ_curl_display), lengthscale=1e10, arrowsize=0.02, color=vec(E_curl_colors))

        # Plot the magnetic field and curl
        arrows!(ax, vec(permutedims(points, [3, 2, 1])), vec(Bₙ_field_display), lengthscale=0.1, arrowsize=0.02, color=vec(B_colors))
        # arrows!(ax, vec(permutedims(points, [3, 2, 1])), vec(Bₙ_curl_display), lengthscale=1, arrowsize=1, color=:green)

        # Plot the point charges with colors based on their charge
        for pc in charge_list
            scatter!(ax, pc.x(t)[1], pc.x(t)[2], pc.x(t)[3], color=RGBf(0, 0.5-pc.q, 0.5+pc.q))
        end

        display(fig)
        sleep(0.001)
    end
end

animate()