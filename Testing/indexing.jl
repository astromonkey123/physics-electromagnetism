const xs = [i for i in -5:0.5:5]; const ys = [i for i in -5:0.5:5]; const zs = [i for i in -5:0.5:5] # Spatial bounds

vector_field = [[x, y, z] for x in xs, y in ys, z in zs]

function index(x, y, z)
    x_pos = (x - xs[begin]) / (xs[end] - xs[begin])
    y_pos = (y - ys[begin]) / (ys[end] - ys[begin])
    z_pos = (z - zs[begin]) / (zs[end] - zs[begin])
    idx_x = trunc(Int, x_pos * (length(xs) - 1) + 1)
    idx_y = trunc(Int, y_pos * (length(xs) - 1) + 1)
    idx_z = trunc(Int, z_pos * (length(xs) - 1) + 1)
    println(vector_field[idx_x, idx_y, idx_z])
    # return findall(i->i==[x, y, z], vector_field)
end

index(-5, -1.5, 3.5)

# println(vector_field[1, end, end])