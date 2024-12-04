using GLMakie; gl = GLMakie
using CalculusWithJulia
using Base.Threads

# Physics constants
# const ϵ₀ = 8.85418782e-12 # Permittivity of free space
# const μ₀ = 1.25663706e-6 # Permeability of free space
# const c = 1/√(ϵ₀⋅μ₀) # Speed of light
const ϵ₀ = 1 # Permittivity of free space (testing)
const μ₀ = 1 # Permeability of free space (testing)
const c = 2 # Speed of light (testing)

# Simulation constants
const dt = 0.1 # Simulation time step
const duration = 20 # Length of the simulation

# Size and resolution
const extent = (x_lower=-5, x_upper=5, y_lower=-5, y_upper=5, z_lower=-5, z_upper=5)
const field_resolution = 15
const operator_resolution = 15

# Plot settings
const dark_latexfonts = merge(theme_dark(), theme_latexfonts())
const fig = Figure(size=(700, 700))
ax = Axis3(fig[1,1], limits=(extent.x_lower, extent.x_upper, extent.y_lower, extent.y_upper, extent.z_lower, extent.z_upper))

# Bounds for the EM field
const xs = LinRange(extent.x_lower, extent.x_upper, field_resolution)
const ys = LinRange(extent.y_lower, extent.y_upper, field_resolution)
const zs = LinRange(extent.z_lower, extent.z_upper, field_resolution)

# Bounds for the curl and divergence fields
const map_xs = LinRange(extent.x_lower, extent.x_upper, operator_resolution)
const map_ys = LinRange(extent.y_lower, extent.y_upper, operator_resolution)
const map_zs = LinRange(extent.z_lower, extent.z_upper, operator_resolution)

# Fields YAY
magnetic_field = [Vec3f(0, 0, 0) for x in map_xs for y in map_ys for z in map_zs]

mutable struct point_charge
    x
    q
end

function E(x, y, z, pc::point_charge, t)
    q = pc.q
    pos = pc.x

    radius = sqrt((x - pos(t)[1])^2 + (y - pos(t)[2])^2 + (z - pos(t)[3])^2)

    t_delay = t - (radius / c)

    radius_delay = sqrt((x - pos(t_delay)[1])^2 + (y - pos(t_delay)[2])^2 + (z - pos(t_delay)[3])^2)

    u = q * (1/radius) * (1/radius_delay) * (x - pos(t_delay)[1])
    v = q * (1/radius) * (1/radius_delay) * (y - pos(t_delay)[2])
    w = q * (1/radius) * (1/radius_delay) * (z - pos(t_delay)[3])

    return [u, v, w]
end

function E_curl(x, y, z, pc::point_charge, t)
    q = pc.q
    pos = pc.x

    function field(x, y, z)
        radius = sqrt((x - pos(t)[1])^2 + (y - pos(t)[2])^2 + (z - pos(t)[3])^2)

        t_delay = t - (radius / c)

        radius_delay = sqrt((x - pos(t_delay)[1])^2 + (y - pos(t_delay)[2])^2 + (z - pos(t_delay)[3])^2)

        u = q * (1/radius_delay) * (x - pos(t_delay)[1])
        v = q * (1/radius_delay) * (y - pos(t_delay)[2])
        w = q * (1/radius_delay) * (z - pos(t_delay)[3])
        return [u, v, w]
    end
    field(v) = field(v...)

    return (∇ × field)(x, y, z)
end

function E_div(x, y, z, pc::point_charge, t)
    q = pc.q
    pos = pc.x

    function field(x, y, z)
        radius = sqrt((x - pos(t)[1])^2 + (y - pos(t)[2])^2 + (z - pos(t)[3])^2)

        t_delay = t - (radius / c)

        radius_delay = sqrt((x - pos(t_delay)[1])^2 + (y - pos(t_delay)[2])^2 + (z - pos(t_delay)[3])^2)

        u = q * (1/radius_delay) * (x - pos(t_delay)[1])
        v = q * (1/radius_delay) * (y - pos(t_delay)[2])
        w = q * (1/radius_delay) * (z - pos(t_delay)[3])
        return [u, v, w]
    end
    field(v) = field(v...)

    return (∇ ⋅ field)(x, y, z)
end

# function B(x, y, z, t)
#     u = 0
#     v = 0
#     w = 0

#     # Take the integral from 0 to t of the curl of the electric field at a point (x, y, z)
#     for t in 0:dt:t
#         for point_charge in charge_list
#             E_dx = E_curl(x, y, z, point_charge, t)[1]
#             E_dy = E_curl(x, y, z, point_charge, t)[2]
#             E_dz = E_curl(x, y, z, point_charge, t)[3]

#             u += E_dx * dt
#             v += E_dy * dt
#             w += E_dz * dt
#         end
#     end

#     return [u, v, w]
# end

function B_curl(x, y, z, pc::point_charge, t, B_field)
    function B(x, y, z)
        x_index = trunc(Int, (14x/10)+(14/2)+1)
        y_index = trunc(Int, (14y/10)+(14/2)+1)
        z_index = trunc(Int, (14z/10)+(14/2)+1)

        # Get ∇×E at (x, y, z)
        electric_curl = E_curl(x, y, z, pc, t)

        # The (u, v, w, t) vector of the magnetic field at (x, y, z, t) is -1 times ∇×E at (x, y, z, t) times dt plus (u, v, w, t-dt)
        # (u, v, w) = (u, v, w) - ∇×E * dt
        u, v, w = reshape(B_field, field_resolution, field_resolution, field_resolution)[x_index, y_index, z_index] - electric_curl * dt
        # global magnetic_field = reshape(updated_B_field, field_resolution^3)

        return [u, v, w]
    end
    B(v) = B(v...)

    return (∇ × B)(x, y, z)
end

function x1(t)
    if t < 0
        t = 0
    end
    [0, 0, 3sin(t)]
end

pc_1 = point_charge(x1, -1.0)
charge_list = [pc_1]

for t in 0:dt:duration
    set_theme!(dark_latexfonts)

    empty!(ax)

    points = [Point3f(x, y, z) for x in xs for y in ys for z in zs]

    electric_field = [Vec3f(0, 0, 0) for x in xs for y in ys for z in zs]
    electric_field_curl = [Vec3f(0, 0, 0) for x in map_xs for y in map_ys for z in map_zs]
    electric_field_div = zeros(operator_resolution, operator_resolution, operator_resolution)

    magnetic_field_curl = [Vec3f(0, 0, 0) for x in map_xs for y in map_ys for z in map_zs]

    xf = [0, 0, 0]
    xi = [0, 0, 0]
    for point_charge in charge_list
        electric_field += [E(x, y, z, point_charge, t) for x in xs for y in ys for z in zs]
        electric_field_curl += [E_curl(x, y, z, point_charge, t) for x in map_xs for y in map_ys for z in map_zs]
        electric_field_div += [E_div(x, y, z, point_charge, t) for x in map_xs, y in map_ys, z in map_zs]

        magnetic_field_curl += [B_curl(x, y, z, pc_1, t, magnetic_field) for x in map_xs for y in map_ys for z in map_zs]
        xf += [point_charge.x(t)[1], point_charge.x(t)[2], point_charge.x(t)[3]]
        xi += [point_charge.x(t-dt)[1], point_charge.x(t-dt)[2], point_charge.x(t-dt)[3]]
    end
    v = (xf - xi) / dt

    # Maxwell-Faraday
    magnetic_field -= electric_field_curl * dt

    electric_field_div_flat = reshape(electric_field_div, operator_resolution^3)
    
    # Current density
    J = [div * v for div in electric_field_div_flat]

    # Maxwell-Ampere
    # electric_field += (((magnetic_field_curl / μ₀) - J) / ϵ₀) * dt
    electric_field += (magnetic_field_curl - J) * dt

    # electric_field += (((magnetic_field_curl / μ₀) - J) / ϵ₀) * dt
    # J = ρv, Current density is the product of charge density and the velocity of the charge
    # We know charge density is related to the divergence of the electric field
    # ∇⋅E = ρ/ϵ₀
    # At a given point, ρ = (∇⋅E) * ϵ₀
    # J = ((∇⋅E) * ϵ₀) * v, where (∇⋅E) * ϵ₀ is some scalar, v is (∂x/∂t, ∂y/∂t, ∂z/∂t)
    # Adding this all into the original equation,
    # electric_field += (((magnetic_field_curl / μ₀) - ((∇⋅E) * ϵ₀) * v) / ϵ₀) * dt
    # C urrently ϵ₀ and μ₀ are 1 so that reduces to
    # electric_field += (magnetic_field_curl - (electric_field_divergence * v)) * dt, for v = (vx, vy, vz)
    # TODO: We need to find magnetic_field_curl and v and then we are done!

    display_electric_field = electric_field
    display_magnetic_field = magnetic_field
    # display_electric_field = [Vec3f(0, 0, vec[3]) for vec in electric_field]
    # display_magnetic_field = [Vec3f(vec[1], vec[2], 0) for vec in magnetic_field]

    electric_field_lengths = norm.(electric_field)
    electric_field_colors = [RGBAf(1, 0, 0, √l - 0.5) for l in electric_field_lengths]

    magnetic_field_lengths = norm.(magnetic_field)
    magnetic_field_colors = [RGBAf(0, √l - 1, 1, √l - 1) for l in magnetic_field_lengths]

    a = (magnetic_field_curl - J) * dt
    a_lengths = norm.(a)
    a_colors = [RGBAf(l - 0.5, 1, 0, l - 0.5) for l in a_lengths]
    println(typeof(magnetic_field_curl))

    arrows!(ax, points, display_electric_field, color=electric_field_colors, lengthscale=0.2, linewidth=0.015, arrowsize = Vec3f(0.1, 0.1, 0.1))
    arrows!(ax, points, display_magnetic_field, color=magnetic_field_colors, lengthscale=0.05, linewidth=0.015, arrowsize = Vec3f(0.1, 0.1, 0.1))
    arrows!(ax, points, magnetic_field_curl, color=a_colors, lengthscale=0.1, linewidth=0.015, arrowsize = Vec3f(0.1, 0.1, 0.1))
    # arrows!(ax, points, (magnetic_field_curl - J) * dt, color=a_colors, lengthscale=0.1, linewidth=0.015, arrowsize = Vec3f(0.1, 0.1, 0.1))
    # arrows!(ax, points, electric_field_div_vec, color=electric_field_div_colors, lengthscale=0.1, linewidth=0.015, arrowsize = Vec3f(0.1, 0.1, 0.1))
    for point_charge in charge_list
        scatter!(ax, point_charge.x(t)[1], point_charge.x(t)[2], point_charge.x(t)[3], color=:white)
    end
    display(fig)

    sleep(dt/10)
end