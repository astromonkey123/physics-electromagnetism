using GLMakie

c = 1
x(t) = [t - 5, 1, 0]
P = [0, 0, 0]

function get_delayed_distance(P, x, t)
    r(t) = hypot( (P - x(t))...)
    t_guess = t
    error(guess) = r(guess) - hypot( (P - x(t - r(guess)/c))...)
    # Newton's method to find the root (minimizing error)
    Δt = 1e-10
    for _ in 1:10
        δt = (error(t_guess + Δt) - error(t_guess)) / Δt
        t_guess = t_guess - ( error(t_guess) / δt)
    end
    return t_guess
end

fig = Figure()
ax = Axis(fig[1, 1], ylabel="Distance (m)", xlabel="Time (s)")
times = []
true_dists = []
observer_dists = []

for t in 0:0.05:10
    true_dist = hypot(( P - x(t))...)
    t_delayed = get_delayed_distance(P, x, t)
    observer_dist = hypot(( P - x(t_delayed))...)

    append!(times, [t])
    append!(true_dists, [true_dist])
    append!(observer_dists, [observer_dist])
end

lines!(ax, times, true_dists, color=:black)
lines!(ax, times, observer_dists, color=:blue)
display(fig)
