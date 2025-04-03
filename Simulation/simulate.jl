# Step 1: Find E based on the location of point charges
# Step 2: Find ∇⋅E and ∇×E
# Step 3: Update B based on ∇×E (Faraday)
# Step 4: Find ∇⋅B and ∇×B
# Step 5: Update E based on ∇×B (Ampere)

# Physics parameters
const c = 299792458 # Speed of light
const ϵ₀ = 8.854187817e-12 # Permittivity of space
const μ₀ = 1.25663706127e-6 # Permeability of space

# Simulation parameters
const dt = 1e-9 # Simulation time step
const steps = 250 # Simulation step number
const xs = [i for i in -10:10]; const ys = [i for i in -10:10]; const zs = [i for i in -10:10] # Simulation spatial bounds

using CalculusWithJulia
using GLMakie
using JLD2

mutable struct PointCharge
    x::Function
    q::Float64
end

# ---------------------------------------------------------------------------------------------------- #

function find_electric_field(x, y, z, t)
    # Parameters:
    #   - x::Float64 | X coordinate
    #   - y::Float64 | Y coordinate
    #   - z::Float64 | Z coordinate
    #   - t::Float64 | Time coordinate
    # Returns:
    #   - Vector{Float64} | Electric field vector
    #   - Float64 | Electric field divergence
    #   - Vector{Float64} | Electric field curl
    electric_vector = zeros(3)
    electric_div = 0
    electric_curl = zeros(3)

    for point_charge in charge_list
        q = point_charge.q
        pos = point_charge.x

        function E(x₀, y₀, z₀)
            radius = sqrt((x₀ - pos(t)[1])^2 + (y₀ - pos(t)[2])^2 + (z₀ - pos(t)[3])^2)

            t_delay = t - (radius / c)

            radius_delay = sqrt((x₀ - pos(t_delay)[1])^2 + (y₀ - pos(t_delay)[2])^2 + (z₀ - pos(t_delay)[3])^2)

            u = q * (1/radius) * (1/radius_delay) * (x₀ - pos(t_delay)[1])
            v = q * (1/radius) * (1/radius_delay) * (y₀ - pos(t_delay)[2])
            w = q * (1/radius) * (1/radius_delay) * (z₀ - pos(t_delay)[3])
            # u = q * (1/radius_delay) * (x₀ - pos(t_delay)[1])
            # v = q * (1/radius_delay) * (y₀ - pos(t_delay)[2])
            # w = q * (1/radius_delay) * (z₀ - pos(t_delay)[3])
            return [u, v, w]
        end
        E_field(v) = E(v...)

        electric_vector += E(x, y, z)
        electric_div += (∇⋅E_field)(x, y, z)
        electric_curl += (∇×E_field)(x, y, z)
    end

    return electric_vector, electric_div, electric_curl
end

function find_magnetic_field(x, y, z, t, magnetic_field, i)
    # Parameters:
    #   - x::Float64 | X coordinate
    #   - y::Float64 | Y coordinate
    #   - z::Float64 | Z coordinate
    #   - t::Float64 | Time coordinate
    #   - magnetic_field::Vector{Vector{Float64}} | The magnetic field
    #   - i::Int64 | Index of current coordinates in the field
    # Returns:
    #   - Vector{Float64} | Magnetic field vector
    #   - Float64 | Magnetic field divergence
    #   - Vector{Float64} | Magnetic field curl
    magnetic_vector = zeros(3)
    magnetic_div = 0
    magnetic_curl = zeros(3)

    function B(x₀, y₀, z₀)
        index = trunc(Int64, (21*21)*(x₀ + 10) + 21(y₀ + 10) + (z₀ + 11))

        # Get ∇×E at (x, y, z)
        _, _, electric_curl = find_electric_field(x₀, y₀, z₀, t)

        # the x y z of the new magnetic field is the x y z of the old magnetic field minus the curl times dt
        # This is Faraday's Law
        u = magnetic_field[index][1] - (electric_curl[1] * dt)
        v = magnetic_field[index][2] - (electric_curl[2] * dt)
        w = magnetic_field[index][3] - (electric_curl[3] * dt)
        return [u, v, w]
    end
    B_field(v) = B(v...)

    magnetic_vector = B(x, y, z)
    magnetic_div += 0
    magnetic_curl += (∇×B_field)(x, y, z)

    return magnetic_vector, magnetic_div, magnetic_curl
end

function generate_electric_field(t, electric_field, electric_div, electric_curl)
    coords = [[x, y, z] for x in xs for y in ys for z in zs]
    new_electric_field = [Vec3f(0, 0, 0) for x in xs for y in ys for z in zs] # Vector field

    for i in eachindex(coords)
        field_vector, div_scalar, curl_vector = find_electric_field(coords[i][1], coords[i][2], coords[i][3], t)
        new_electric_field[i] = field_vector
        electric_div[i] = div_scalar
        electric_curl[i] = curl_vector
    end

    return new_electric_field, electric_div, electric_curl
end

function generate_magnetic_field(t, magnetic_field, magnetic_div, magnetic_curl)
    coords = [[x, y, z] for x in xs for y in ys for z in zs]
    new_magnetic_field = [Vec3f(0, 0, 0) for x in xs for y in ys for z in zs] # Vector field

    for i in eachindex(coords)
        field_vector, div_scalar, curl_vector = find_magnetic_field(coords[i][1], coords[i][2], coords[i][3], t, magnetic_field, i)
        new_magnetic_field[i] = field_vector
        magnetic_div[i] = div_scalar
        magnetic_curl[i] = curl_vector
    end

    return new_magnetic_field, magnetic_div, magnetic_curl
end

# ---------------------------------------------------------------------------------------------------- #

function x1(t)
    if t < 0
        t = 0
    end
    [0, 0, 2sin(0.25*c*t)]
end

pc_1 = PointCharge(x1, -1.0)

charge_list = [pc_1]

magnetic_field = [Vec3f(0, 0, 0) for x in xs for y in ys for z in zs] # Vector field
magnetic_div = [0.0 for x in xs for y in ys for z in zs] # Scalar field
magnetic_curl = [Vec3f(0, 0, 0) for x in xs for y in ys for z in zs] # Vector field

electric_field = [Vec3f(0, 0, 0) for x in xs for y in ys for z in zs] # Vector field
electric_div = [0.0 for x in xs for y in ys for z in zs] # Scalar field
electric_curl = [Vec3f(0, 0, 0) for x in xs for y in ys for z in zs] # Vector field

for t in 0:dt:(steps*dt)
    println("Simulated $(round((100 * t / (steps * dt)), digits=2))%")

    # Create the electric field, as well as its divergence and curl
    electric_field, electric_div, electric_curl = generate_electric_field(t, electric_field, electric_div, electric_curl)

    magnetic_field, magnetic_div, magnetic_curl = generate_magnetic_field(t, magnetic_field, magnetic_div, magnetic_curl)

    # Farady's law is already accoounted for (it is used in the generation of the magnetic field)

    # Ampere's law
    electric_field += magnetic_curl * dt * (1/ϵ₀) * (1/μ₀)

    electric_field_lengths = norm.(electric_field)
    electric_field_colors = [RGBAf(1, 0.2, 0.2, l - 0.25) for l in electric_field_lengths]

    magnetic_field_lengths = norm.(magnetic_field)
    magnetic_field_colors = [RGBAf(0.2, (l * 1e9) - 0.25, 1, (l * 1e9) - 0.25) for l in magnetic_field_lengths]

    save("cache/electric_field$(trunc(Int, (t/dt))).jld2", "E$(trunc(Int, (t/dt)))", electric_field)
    save("cache/magnetic_field$(trunc(Int, (t/dt))).jld2", "B$(trunc(Int, (t/dt)))", magnetic_field)
end