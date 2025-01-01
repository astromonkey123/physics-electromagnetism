using CalculusWithJulia
using GLMakie

xs = collect(range(-2π, 2π, 17))
ys = collect(range(-2π, 2π, 17))

# Scalar valued function f:R² -> R
f(x, y) = 3x + sin(y)
F(v) = f(v...)

# Vector valued function g:R² -> R²
g(x, y) = [y * cos(x), -x * sin(y)]
G(v) = g(v...)

h(x, y) = Point2(g(x, y)[1], g(x, y)[2])

points = [Point2f(x, y) for x in xs for y in ys]
vectors = [g(x, y) for x in xs for y in ys]

div = [(∇⋅G)(x, y) for x in xs, y in ys]
curl = [(∇×G)(x, y) for x in xs, y in ys]

fig = Figure()
ax = Axis(fig[1, 1], limits=(xs[begin], xs[end], ys[begin], ys[end]))

# div_map = heatmap!(ax, xs, ys, transpose(div), colormap=:magma, interpolate=true)
curl_map = heatmap!(ax, xs, ys, transpose(curl), colormap=:magma, interpolate=true)

arrows!(ax, points, vectors, lengthscale=0.2, color=:white)

# streamplot!(ax, h, xs, ys, colormap=:magma)

Colorbar(fig[2, 1], div_map, vertical=false)

fig