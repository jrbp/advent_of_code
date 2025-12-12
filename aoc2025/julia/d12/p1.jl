module AOC25_11
macro comment(e...) end

function parse_region(ln)
    s_dims, s_nums = split(ln, ":")
    x, y = parse.(Int, split(s_dims, "x"))
    reqs = parse.(Int, split(s_nums))
    @assert length(reqs) == 6
    (;x, y, reqs)
end

function readinput(filename)
    line_iter = eachline(filename)
    ln, state = iterate(line_iter)
    shapes = BitMatrix[]
    # always 6 shapes on 3x3 grid
    while length(shapes) < 6
        @assert ln == "$(length(shapes)):"
        push!(shapes, mapfoldl(vcat, 1:3) do _
                  ln, state = iterate(line_iter, state)
                  map(==('#'), collect(ln))'
              end)
        ln, state = iterate(line_iter, state)
        @assert ln == ""
        ln, state = iterate(line_iter, state)
    end
    #regions = @NamedTuple{x::Int, y::Int, reqs::Vector{Int}}[]
    regions = [parse_region(ln),]
    while !Base.isdone(line_iter, state)
        ln, state = iterate(line_iter, state)
        push!(regions, parse_region(ln))
    end
    shapes, regions
end
function canfit(shapes, region)
    tot_region_size = region.x * region.y
    tot_shape_size = mapreduce(+, shapes, region.reqs) do s, n
        n * sum(s)
    end
    tot_shape_size > tot_region_size && return false
    ((region.x ÷ 3) * (region.y ÷ 3)) >= sum(region.reqs) && return true
    error("math is hard")
end
canfit(shapes) = r -> canfit(shapes, r)

function main_1(filename="d12/input")
    shapes, regions = readinput(filename)
    sum(canfit(shapes), regions)
end
@assert main_1("d12/input") |> isequal(476) # lol

end

# Local Variables:
# julia-snail-repl-buffer: "*julia-aoc*"
# julia-snail-port: 10012
# End:
