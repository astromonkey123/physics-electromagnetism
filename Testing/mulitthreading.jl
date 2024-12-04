using Base.Threads

xs = collect(-2:2)
ys = collect(-2:2)
zs = collect(-2:2)

coords = [[x, y, z] for x in xs for y in ys for z in zs]

function square()
    operation = [[0, 0, 0] for _ in eachindex(coords)]

    for i in eachindex(coords)
        operation[i] = coords[i] .^ 2
    end
    return operation
end

function square_mt()
    num_splits = 8

    len = length(coords)
    
    sections = []
    
    # Split coords into sections
    for i in 1:num_splits
        start_index = trunc(Int, round(((i - 1) / num_splits) * len) + 1)
        end_index = trunc(Int, round((i / num_splits) * len))
        append!(sections, [coords[start_index:end_index]])
    end
    
    operation = [[], [], [], [], [], [], [], []]
    
    @threads for i in eachindex(sections)
        for j in eachindex(sections[i])
            append!(operation[i], [sections[i][j] .^ 2])
        end
    end
    return reduce(vcat, operation)
end

a = square()
b = square_mt()