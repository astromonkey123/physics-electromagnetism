using GLMakie; gl = GLMakie
using CalculusWithJulia; calcjl = CalculusWithJulia
using LinearAlgebra; la = LinearAlgebra

const c = 1

const dt = 0.1

function normalize(u, v, norm)
    u_norm = zeros(50, 50)
    v_norm = zeros(50, 50)

    for i in 1:length(u)
        for j in 1:length(u[i])
            u_index = u[i][j]
            v_index = v[i][j]
            radius = hypot(u_index, v_index)
            scale = norm / radius

            u_norm[i][j] = scale * u[i][j]
            v_norm[i][j] = scale * v[i][j]
        end
    end
    return u_norm, v_norm
end

mutable struct point_charge
    q
    pos
end

function x₁(t)
    if t < 0
        t = 0
    end
    position = [0, sin(t)]
    # Check to ensure the displacement is less than the speed of light
    return position
end

charge_1 = point_charge(1.0, x₁)

charge_1 = Observable(charge_1)

charge_list = [charge_1]

function E(x, y, charge::point_charge, t)
    pos = charge.pos
    q = charge.q

    radius = √((x - pos(t)[1])^2 + (y - pos(t)[2])^2)

    t_delayed = t - (radius / c)

    radius_delayed = √((x - pos(t_delayed)[1])^2 + (y - pos(t_delayed)[2])^2)

    u = -q * (1 / radius) * (1 / radius_delayed) * (x - pos(t_delayed)[1])
    v = -q * (1 / radius) * (1 / radius_delayed) * (y - pos(t_delayed)[2])

    return [u, v]
end

xs = LinRange(-25, 25, 50)
ys = LinRange(-25, 25, 50)

dark_latexfonts = merge(theme_dark(), theme_latexfonts())

fig = Figure(size = (700, 700))
ax = Axis(fig[1, 1], limits = (-25, 25, -25, 25))

for t in 1:dt:50
    set_theme!(dark_latexfonts)
    empty!(ax)

    us = zeros(50, 50)
    vs = zeros(50, 50)
    us_delayed = zeros(50, 50)
    vs_delayed = zeros(50, 50)

    for charge in charge_list
        us += [E(x, y, charge[], t)[1] for x in xs, y in ys]
        vs += [E(x, y, charge[], t)[2] for x in xs, y in ys]
        scatter!(ax, charge[].pos(t)[1], charge[].pos(t)[2], color = :white)
    end

    us_norm = (1 ./ sqrt.(us .^ 2 .+ vs .^ 2)) .* us
    vs_norm = (1 ./ sqrt.(us .^ 2 .+ vs .^ 2)) .* vs

    # Encode the strength of the vector with color
    strength = vec(sqrt.(vs .^ 2))
    # strength = vec(log10.(vs .^ 2))
    # strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))

    # Normalize a direction
    # us = us_norm
    # vs = vs_norm

    # Set a direction to 0
    # us = [0 for _ in us]
    # vs = [0 for _ in vs]

    arrows!(ax, xs, ys, us, vs, lengthscale = 0.7, arrowsize = 5, arrowcolor = strength, linecolor = strength)

    display(fig)

    sleep(0.01)
end