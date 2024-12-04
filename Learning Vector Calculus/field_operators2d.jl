using GLMakie; gl = GLMakie
using CalculusWithJulia; calcjl = CalculusWithJulia

const dt = 0.1
const c = 0.5

const xs = LinRange(-10, 10, 50)
const ys = LinRange(-10, 10, 50)

const map_xs = LinRange(-10, 10, 50)
const map_ys = LinRange(-10, 10, 50)

mutable struct point_charge
    x
    m
end

function position1(t)
    if t < 0
        [0, 0, 0]
    else
        [0, sin(t/2), 0]
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

    radius = hypot(x - pos(t)[1], y - pos(t)[2], z - pos(t)[3])

    t_delay = t - (radius / c)

    radius_delay = hypot(x - pos(t_delay)[1], y - pos(t_delay)[2], z - pos(t_delay)[3])

    u = m * (1/radius) * (1/radius_delay) * (x - pos(t_delay)[1])
    v = m * (1/radius) * (1/radius_delay) * (y - pos(t_delay)[2])
    w = m * (1/radius) * (1/radius_delay) * (z - pos(t_delay)[3])

    return [u, v, w]
end

function get_divergence(x, y, z, pm::point_charge, t)
    m = pm.m
    pos = pm.x

    function field(x, y, z)
        radius = hypot(x - pos(t)[1], y - pos(t)[2], z - pos(t)[3])

        t_delay = t - (radius / c)

        radius_delay = hypot(x - pos(t_delay)[1], y - pos(t_delay)[2], z - pos(t_delay)[3])

        u = m * (1/radius_delay) * (x - pos(t_delay)[1])
        v = m * (1/radius_delay) * (y - pos(t_delay)[2])
        w = m * (1/radius_delay) * (z - pos(t_delay)[3])
        return [u, v, w]
    end
    field(v) = field(v...)

    return (∇ ⋅ field)(x, y, z)
end

function get_curl(x, y, z, pm::point_charge, t)
    m = pm.m
    pos = pm.x

    function field(x, y, z)
        radius = hypot(x - pos(t)[1], y - pos(t)[2], z - pos(t)[3])

        t_delay = t - (radius / c)

        radius_delay = hypot(x - pos(t_delay)[1], y - pos(t_delay)[2], z - pos(t_delay)[3])

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
ax = Axis(fig[1,1], limits=(-10, 10, -10, 10))

for t in 0:dt:30
    set_theme!(dark_latexfonts)

    empty!(ax)

    us = zeros(50, 50)
    vs = zeros(50, 50)
    divergence_map = zeros(50, 50)
    curl_map = zeros(50, 50)

    for point_charge in charge_list
        us += [E(x, y, 0, point_charge, t)[1] for x in xs, y in ys]
        vs += [E(x, y, 0, point_charge, t)[2] for x in xs, y in ys]
        divergence_map += [get_divergence(x, y, 0, point_charge, t) for x in map_xs, y in map_ys]
        curl_map += [get_curl(x, y, 0, point_charge, t)[3] for x in map_xs, y in map_ys]
    end

    strength = vec(sqrt.(us.^2 + vs.^2))

    image!(ax, (-10, 10), (-10, 10), divergence_map, colormap=:curl, colorrange=(-5, 5))
    # image!(ax, (-10, 10), (-10, 10), curl_map, colormap=:curl, colorrange=(-2.5, 2.5))
    arrows!(ax, xs, ys, us, vs, colormap=:viridis, arrowcolor=strength, linecolor=strength)
    for point_charge in charge_list
        scatter!(ax, point_charge.x(t)[1], point_charge.x(t)[2], color=:green)
    end
    display(fig)

    sleep(0.01)
end