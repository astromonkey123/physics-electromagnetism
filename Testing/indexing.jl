const xs = [i for i in -5:0.5:5]; const ys = [i for i in -5:0.5:5]; const zs = [i for i in -5:0.5:5] # Spatial bounds
using Base.Threads

points = [[x, y, z] for x in xs, y in ys, z in zs]

function pos_to_idx(x, x_min, x_max, x_step)
    if x < x_min || x > x_max; error("Position out of bounds"); end
    index = (x - x_min)/x_step + 1
    return trunc(Int, index)
end

println(pos_to_idx(0, -5, 5, 0.5))

vectors = vec([zeros(3) for x in xs, y in ys, z in zs])

num_threads = 8
len = length(vec(vector_field))
split_len = ceil(len / num_threads)

index_splits = append!([trunc(Int, n) for n in 1:split_len:len], len)

positions = [[0, 1, 0], [0, -1, 0], [1, 0, 0], [-1, 0, 0]]

@threads for i in 1:num_threads
    start_idx = index_splits[i]
    end_idx = index_splits[i + 1]
    working_points = vec(points)[start_idx:end_idx]
    for pos in positions
        for j in eachindex(working_points)
            point = working_points[j]
            vectors[j + start_idx - 1] = pos - point / hypot((pos - point)...)
        end
    end
end
