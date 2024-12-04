using GLMakie; gl = GLMakie
using CalculusWithJulia; calcjl = CalculusWithJulia

const dt = 0.1
const c = 0.5

const extent = (x_lower=-5, x_upper=5, y_lower=-5, y_upper=5, z_lower=-5, z_upper=5)
const field_resolution = 15
const operator_resolution = 15

const xs = LinRange(extent.x_lower, extent.x_upper, field_resolution)
const ys = LinRange(extent.y_lower, extent.y_upper, field_resolution)
const zs = LinRange(extent.z_lower, extent.z_upper, field_resolution)

const map_xs = LinRange(extent.x_lower, extent.x_upper, operator_resolution)
const map_ys = LinRange(extent.y_lower, extent.y_upper, operator_resolution)
const map_zs = LinRange(extent.z_lower, extent.z_upper, operator_resolution)

mutable struct point_charge
    x
    m
end

function position1(t)
    if t < 0
        [0, 0, 0]
    else
        [0, 5sin(t/5), 0]
    end
end

function position2(t)
    if t < 0
        [(5-t)cos(π+t*0.2), (5-t)sin(π+t*0.2), 0]
    else
        [(5-t)cos(π+t*0.2), (5-t)sin(π+t*0.2), 0]
    end
end

mass1 = point_charge(position1, 1.0)
mass2 = point_charge(position2, 1.0)
charge_list = [mass1]

function E(x, y, z, pm::point_charge, t)
    m = pm.m
    pos = pm.x

    radius = sqrt((x - pos(t)[1])^2 + (y - pos(t)[2])^2 + (z - pos(t)[3])^2)

    t_delay = (t - radius) / c

    radius_delay = sqrt((x - pos(t_delay)[1])^2 + (y - pos(t_delay)[2])^2 + (z - pos(t_delay)[3])^2)

    u = m * (1/radius) * (1/radius_delay) * (x - pos(t_delay)[1])
    v = m * (1/radius) * (1/radius_delay) * (y - pos(t_delay)[2])
    w = m * (1/radius) * (1/radius_delay) * (z - pos(t_delay)[3])

    return [u, v, w]
end

function B(x, y, z, E, E_curl, t)
    # ∇⋅B=0
    # ∇×B=μ(J+∂E/∂t)
    # And ∇×E=-∂B/∂t
    # If we know B(t=0) then the rate of change will tell us how it evolves.
    # What is B(t=0)?

    u = E_curl[1] * t
    v = E_curl[2] * t
    w = E_curl[3] * t

    return [u, v, w]
end

function E_div(x, y, z, pm::point_charge, t)
    m = pm.m
    pos = pm.x

    function field(x, y, z)
        radius = sqrt((x - pos(t)[1])^2 + (y - pos(t)[2])^2 + (z - pos(t)[3])^2)

        t_delay = (t - radius) / c

        radius_delay = sqrt((x - pos(t_delay)[1])^2 + (y - pos(t_delay)[2])^2 + (z - pos(t_delay)[3])^2)
    
        # u = m * (1/radius) * (1/radius_delay) * (x - pos(t_delay)[1])
        # v = m * (1/radius) * (1/radius_delay) * (y - pos(t_delay)[2])
        # w = m * (1/radius) * (1/radius_delay) * (z - pos(t_delay)[3])

        u = m * (1/radius_delay) * (x - pos(t_delay)[1])
        v = m * (1/radius_delay) * (y - pos(t_delay)[2])
        w = m * (1/radius_delay) * (z - pos(t_delay)[3])
        return [u, v, w]
    end
    field(v) = field(v...)

    return (∇ ⋅ field)(x, y, z)
end

function E_curl(x, y, z, pm::point_charge, t)
    m = pm.m
    pos = pm.x

    function field(x, y, z)
        radius = sqrt((x - pos(t)[1])^2 + (y - pos(t)[2])^2 + (z - pos(t)[3])^2)

        t_delay = (t - radius) / c

        radius_delay = sqrt((x - pos(t_delay)[1])^2 + (y - pos(t_delay)[2])^2 + (z - pos(t_delay)[3])^2)
    
        # u = m * (1/radius) * (1/radius_delay) * (x - pos(t_delay)[1])
        # v = m * (1/radius) * (1/radius_delay) * (y - pos(t_delay)[2])
        # w = m * (1/radius) * (1/radius_delay) * (z - pos(t_delay)[3])

        u = m * (1/radius_delay) * (x - pos(t_delay)[1])
        v = m * (1/radius_delay) * (y - pos(t_delay)[2])
        w = m * (1/radius_delay) * (z - pos(t_delay)[3])
        return [u, v, w]
    end
    field(v) = field(v...)    

    return (∇ × field)(x, y, z)
end

dark_latexfonts = merge(theme_dark(), theme_latexfonts())

fig = Figure(size=(700, 700))
ax = Axis3(fig[1,1], limits=(extent.x_lower, extent.x_upper, extent.y_lower, extent.y_upper, extent.z_lower, extent.z_upper))

magnetic_field = [Vec3f(0, 0, 0) for x in map_xs for y in map_ys for z in map_zs]

for t in 0:dt:10
    set_theme!(dark_latexfonts)

    empty!(ax)

    points = [Point3f(x, y, z) for x in xs for y in ys for z in zs]

    electric_field = [Vec3f(0, 0, 0) for x in xs for y in ys for z in zs]
    electric_field_div = zeros(operator_resolution, operator_resolution, operator_resolution)
    electric_field_curl = [Vec3f(0, 0, 0) for x in map_xs for y in map_ys for z in map_zs]

    for point_charge in charge_list
        electric_field += [E(x, y, z, point_charge, t) for x in xs for y in ys for z in zs]
        electric_field_div += [E_div(x, y, z, point_charge, t) for x in map_xs, y in map_ys, z in map_zs]
        electric_field_curl += [E_curl(x, y, z, point_charge, t) for x in map_xs for y in map_ys for z in map_zs]
    end

    magnetic_field -= electric_field_curl * dt

    electric_field_lengths = norm.(electric_field)
    electric_field_colors = [RGBAf(1, 0, 0, l^2 - 0.2) for l in electric_field_lengths]

    electric_field_curl_lengths = norm.(electric_field_curl)
    electric_field_curl_colors = [RGBAf(0, 1, 0, l^2 - 0.2) for l in electric_field_curl_lengths]

    magnetic_field_lengths = norm.(magnetic_field)
    magnetic_field_colors = [RGBAf(0, 0, 1, l^2 - 0.2) for l in magnetic_field_lengths]

    arrows!(ax, points, electric_field, colormap=:viridis, color=electric_field_colors, lengthscale=0.2, linewidth=0.035, arrowsize = Vec3f(0.15, 0.15, 0.15))
    arrows!(ax, points, electric_field_curl, colormap=:viridis, color=electric_field_curl_colors, lengthscale=0.5, linewidth=0.025, arrowsize = Vec3f(0.1, 0.1, 0.1))
    # arrows!(ax, points, magnetic_field, colormap=:curl, color=magnetic_field_colors, lengthscale=0.5, linewidth=0.025, arrowsize = Vec3f(0.15, 0.15, 0.15))
    for point_charge in charge_list
        scatter!(ax, point_charge.x(t)[1], point_charge.x(t)[2], point_charge.x(t)[3], color=:blue)
    end
    display(fig)

    sleep(0.01)
end