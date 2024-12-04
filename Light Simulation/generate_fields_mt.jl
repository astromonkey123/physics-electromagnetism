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
const extent_x = 10 # Simualtion size
const extent_y = 10
const extent_z = 10
const xs = [i for i in -extent_x:extent_x]; const ys = [i for i in -extent_y:extent_y]; const zs = [i for i in -extent_z:extent_z] # Simulation spatial bounds

using CalculusWithJulia
using GLMakie
using Base.Threads

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

            return [u, v, w]
        end
        E_field(v) = E(v...)

        electric_vector += E(x, y, z)
        electric_div += (∇⋅E_field)(x, y, z)
        electric_curl += (∇×E_field)(x, y, z)
    end

    return electric_vector, electric_div, electric_curl
end

function find_magnetic_field(x, y, z, t, magnetic_field)
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
        index = trunc(Int64, ((2 * extent_y + 1) * (2 * extent_z + 1))*(x₀ + extent_x) + (2 * extent_z + 1)*(y₀ + extent_y) + (z₀ + extent_z + 1))

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

    num_splits = 8 # 8 for now
    len = length(coords)
    sections = []
    field = [[], [], [], [], [], [], [], []]
    div = [[], [], [], [], [], [], [], []]
    curl = [[], [], [], [], [], [], [], []]

    # Split coords into sections
    for i in 1:num_splits
        start_index = trunc(Int, round(((i - 1) / num_splits) * len) + 1)
        end_index = trunc(Int, round((i / num_splits) * len))
        append!(sections, [coords[start_index:end_index]])
    end

    @threads for i in eachindex(sections)
        for j in eachindex(sections[i])
            field_vector, div_scalar, curl_vector = find_electric_field(sections[i][j][1], sections[i][j][2], sections[i][j][3], t)
            append!(field[i], [Vec3f(field_vector[1], field_vector[2], field_vector[3])])
            append!(div[i], [div_scalar])
            append!(curl[i], [curl_vector])
        end
    end

    new_electric_field = reduce(vcat, field)
    electric_div = reduce(vcat, div)
    electric_curl = reduce(vcat, curl)

    return new_electric_field, electric_div, electric_curl
end

function generate_magnetic_field(t, magnetic_field, magnetic_div, magnetic_curl)
    coords = [[x, y, z] for x in xs for y in ys for z in zs]

    num_splits = 8 # 8 for now
    len = length(coords)
    sections = []
    field = [[], [], [], [], [], [], [], []]
    div = [[], [], [], [], [], [], [], []]
    curl = [[], [], [], [], [], [], [], []]

    # Split coords into sections
    for i in 1:num_splits
        start_index = trunc(Int, round(((i - 1) / num_splits) * len) + 1)
        end_index = trunc(Int, round((i / num_splits) * len))
        append!(sections, [coords[start_index:end_index]])
    end

    @threads for i in eachindex(sections)
        for j in eachindex(sections[i])
            field_vector, div_scalar, curl_vector = find_magnetic_field(sections[i][j][1], sections[i][j][2], sections[i][j][3], t, magnetic_field)
            append!(field[i], [Vec3f(field_vector[1], field_vector[2], field_vector[3])])
            append!(div[i], [div_scalar])
            append!(curl[i], [curl_vector])
        end
    end

    new_magnetic_field = reduce(vcat, field)
    magnetic_div = reduce(vcat, div)
    magnetic_curl = reduce(vcat, curl)

    return new_magnetic_field, magnetic_div, magnetic_curl
end

# ---------------------------------------------------------------------------------------------------- #

fig = Figure(size=(700, 700))   
ax = Axis3(fig[1, 1], limits=(-extent_x, extent_x, -extent_y, extent_y, -extent_z, extent_z))

points = [Point3f(x, y, z) for x in xs for y in ys for z in zs]

function x1(t)
    if t < 0
        t = 0
    end
    [0, 0, 2*sin(0.25*c*t)]
    # [2*cos(0.3*c*t), 2*sin(0.3*c*t), 0]
end

pc_1 = PointCharge(x1, -1.0)

charge_list = [pc_1]

magnetic_field = [Vec3f(0, 0, 0) for x in xs for y in ys for z in zs] # Vector field
magnetic_div = [0.0 for x in xs for y in ys for z in zs] # Scalar field
magnetic_curl = [Vec3f(0, 0, 0) for x in xs for y in ys for z in zs] # Vector field

electric_field = [Vec3f(0, 0, 0) for x in xs for y in ys for z in zs] # Vector field
electric_div = [0.0 for x in xs for y in ys for z in zs] # Scalar field
electric_curl = [Vec3f(0, 0, 0) for x in xs for y in ys for z in zs] # Vector field

set_theme!(theme_dark())

for t in 0:dt:(steps*dt)

    # Create the electric field, as well as its divergence and curl
    electric_field, electric_div, electric_curl = generate_electric_field(t, electric_field, electric_div, electric_curl)

    magnetic_field, magnetic_div, magnetic_curl = generate_magnetic_field(t, magnetic_field, magnetic_div, magnetic_curl)

    # Farady's law is already accoounted for (it is used in the generation of the magnetic field)

    # Ampere's law
    electric_field += magnetic_curl * dt * (1/ϵ₀) * (1/μ₀)

    # electric_field_display = [Vec3f(0, 0, 0) for x in xs for y in ys for z in zs] # Vector field
    # for i in eachindex(electric_field)
    #     if i % (2 * extent_z + 1) == extent_z + 1
    #         electric_field_display[i] = Vec3f(0, 0, electric_field[i][3])
    #     end
    # end

    # magnetic_field_display = [Vec3f(0, 0, 0) for x in xs for y in ys for z in zs] # Vector field
    # for i in eachindex(magnetic_field)
    #     if i % (2 * extent_z + 1) == extent_z + 1
    #         magnetic_field_display[i] = Vec3f(magnetic_field[i][1], magnetic_field[i][2], 0)
    #     end
    # end

    # electric_field_lengths = norm.(electric_field_display)
    # electric_field_colors = [RGBAf(1, 0.2, 0.2, 1) for l in electric_field_lengths]

    # magnetic_field_lengths = norm.(magnetic_field_display)
    # magnetic_field_colors = [RGBAf(0.2, (l * 1e10) - 0.75, 1, 1) for l in magnetic_field_lengths]

    # for i in eachindex(electric_field_display)
    #     electric_field_display[i] *= 1/electric_field_lengths[i]
    # end

    # for i in eachindex(magnetic_field_display)
    #     magnetic_field_display[i] *= 1/magnetic_field_lengths[i]
    # end

    electric_field_lengths = norm.(electric_field)
    electric_field_colors = [RGBAf(1, 0.2, 0.2, l - 0.25) for l in electric_field_lengths]

    magnetic_field_lengths = norm.(magnetic_field)
    magnetic_field_colors = [RGBAf(0.2, (l * 1e9) - 0.25, 1, (l * 1e9) - 0.25) for l in magnetic_field_lengths]

    empty!(ax)

    # arrows!(ax, points, electric_field_display, color=electric_field_colors, lengthscale=500, linewidth=0.05, arrowsize = Vec3f(0.05, 0.05, 0.05))
    # arrows!(ax, points, magnetic_field_display, color=magnetic_field_colors, lengthscale=5e10, linewidth=0.05, arrowsize = Vec3f(0.05, 0.05, 0.05))

    arrows!(ax, points, electric_field, color=electric_field_colors, lengthscale=0.5, linewidth=0.05, arrowsize = Vec3f(0.1, 0.1, 0.1))
    arrows!(ax, points, magnetic_field, color=magnetic_field_colors, lengthscale=5e8, linewidth=0.05, arrowsize = Vec3f(0.1, 0.1, 0.1))

    for point_charge in charge_list
        scatter!(ax, point_charge.x(t)[1], point_charge.x(t)[2], point_charge.x(t)[3], color=:white)
    end

    # electric_curl_lengths = norm.(electric_curl)
    # electric_curl_colors = [RGBAf(1, 0.8, 0.8, (l * 1e2) - 0.5) for l in electric_curl_lengths]
    # arrows!(ax, points, electric_curl, color=electric_curl_colors, lengthscale=10, linewidth=0.015, arrowsize = Vec3f(0.075, 0.075, 0.075))

    # magnetic_curl_lengths = norm.(magnetic_curl)
    # magnetic_curl_colors = [RGBAf(0, 1, 0, (10*l) - 1) for l in magnetic_curl_lengths]
    # arrows!(ax, points, magnetic_curl, color=magnetic_curl_colors, lengthscale=2, linewidth=0.015, arrowsize = Vec3f(0.075, 0.075, 0.075))

    display(fig)
    sleep(dt)
end