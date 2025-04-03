using GLMakie
using CalculusWithJulia
using JLD

# Physics parameters
const c = 299792458 # Speed of light
const ϵ₀ = 8.854187817e-12 # Permittivity of space
const μ₀ = 1.25663706127e-6 # Permeability of space

# Simulation parameters
const dt = 1e-9 # Simulation time step
const steps = 250 # Simulation step number
const xs = [i for i in -10:10]; const ys = [i for i in -10:10]; const zs = [i for i in -10:10] # Simulation spatial bounds

mutable struct PointCharge
    x::Function
    q::Float64
end

fig = Figure(size=(700, 700))
ax = Axis3(fig[1, 1], limits=(-10, 10, -10, 10, -10, 10))

points = [Point3f(x, y, z) for x in xs for y in ys for z in zs]

function x1(t)
    if t < 0
        t = 0
    end
    [0, 0, 2sin(0.25*c*t)]
end

pc_1 = PointCharge(x1, -1.0)

charge_list = [pc_1]

electric_field_dict = load("cache/electric_field0.jld2")
magnetic_field_dict = load("cache/magnetic_field0.jld2")

for t in 0:dt:(steps*dt)
    println("Loaded $(round((100 * t / (steps * dt)), digits=2))%")
    merge!(electric_field_dict, load("cache/electric_field$(trunc(Int, (t/dt))).jld2"))
    merge!(magnetic_field_dict, load("cache/magnetic_field$(trunc(Int, (t/dt))).jld2"))
end

for t in 0:dt:(steps*dt)
    set_theme!(theme_dark())
    empty!(ax)

    electric_field = [Vec3f(E) for E in electric_field_dict["E$(trunc(Int, (t/dt)))"]]
    magnetic_field = [Vec3f(B) for B in magnetic_field_dict["B$(trunc(Int, (t/dt)))"]]

    electric_field_lengths = norm.(electric_field)
    electric_field_colors = [RGBAf(1, 0.2, 0.2, l - 0.25) for l in electric_field_lengths]

    magnetic_field_lengths = norm.(magnetic_field)
    magnetic_field_colors = [RGBAf(0.2, (l * 1e9) - 0.25, 1, (l * 1e9) - 0.25) for l in magnetic_field_lengths]

    arrows!(ax, points, electric_field, color=electric_field_colors, lengthscale=0.5, linewidth=0.05, arrowsize = Vec3f(0.1, 0.1, 0.1))
    arrows!(ax, points, magnetic_field, color=magnetic_field_colors, lengthscale=5e8, linewidth=0.05, arrowsize = Vec3f(0.1, 0.1, 0.1))

    for point_charge in charge_list
        scatter!(ax, point_charge.x(t)[1], point_charge.x(t)[2], point_charge.x(t)[3], color=:white)
    end

    display(fig)
    sleep(dt)
end

# for t in 0:dt:duration
#     empty!(ax)
#     set_theme!(dark_latexfonts)
    
#     electric_field = [Vec3f(E) for E in electric_field_dict["E$t"]]
#     magnetic_field = [Vec3f(B) for B in magnetic_field_dict["B$t"]]

#     electric_field_lengths = norm.(electric_field)
#     electric_field_colors = [RGBAf(1, 0, 0, l - 0.2) for l in electric_field_lengths]
    
#     # Only display vertical components
#     # for i in eachindex(electric_field)
#     #     electric_field[i] = Vec3f(0.0, 0.0, electric_field[i][3]*5)
#     # end

#     magnetic_field_lengths = norm.(magnetic_field)
#     magnetic_field_colors = [RGBAf(0, √l - 1, 1, l^2 - 2) for l in magnetic_field_lengths]

#     arrows!(ax, points, electric_field, color=electric_field_colors, lengthscale=0.2, linewidth=0.025, arrowsize = Vec3f(0.15, 0.15, 0.15))
#     arrows!(ax, points, magnetic_field, color=magnetic_field_colors, lengthscale=0.1, linewidth=0.015, arrowsize = Vec3f(0.1, 0.1, 0.1))
#     for point_charge in charge_list
#         scatter!(ax, point_charge.x(t)[1], point_charge.x(t)[2], point_charge.x(t)[3], color=:pink)
#     end
#     display(fig)

#     sleep(dt/10)
# end