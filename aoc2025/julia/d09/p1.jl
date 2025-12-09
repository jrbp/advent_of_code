module AOC25_09
using Transducers

readpos(filename) = eachline(filename) |>
    Map(x->split(x, ","; limit=2)) |>
    Map(p->(;x=first(p), y=last(p))) |> # to 2-tuple
    Map(x->map(xi->parse(Int, xi), x)) |>
    collect


function main_1(filename="d09/sample")
    positions = readpos(filename)
    Iterators.product(eachindex(positions), eachindex(positions)) |>
        Filter(splat(<)) |>
        Map(x->getindex.(Ref(positions), x)) |>
        Map() do (p1, p2)
            prod(zip(p1, p2)) do (c1, c2)
                abs(c2 - c1) + 1
            end
        end |>
            foldxl(max; init=0)
end
@assert main_1("d09/sample") |> isequal(50)
@assert main_1("d09/input") |> isequal(4763932976)

function main_2(filename="d09/sample")
    positions = readpos(filename)
    miny, maxy = extrema(last.(positions))
    #println(miny)
    vlines = sort(positions; by=first) |> PartitionBy(first) |> Map(copy) |> Map() do pts
        first(pts).x => (pts |> Map(last) |> foldxl(TeeRF(min, max)) |> splat(range))
        end |> collect
    #println(vlines)
    allowed_x = miny:maxy |> Map() do yy
        bounding_lines = vlines |> Filter(r->yy ∈ last(r)) |> collect
        isempty(bounding_lines) ? (1:0) : bounding_lines |>
            Map(first) |> extrema |> splat(range)
    end |> collect
    allowedpt(p) = p.x ∈ allowed_x[p.y - miny + 1]
    # allowedpt_lz(p) = let blines = filter(r -> p.y ∈ last(r), vlines)
    #     isempty(blines) ? false : p.x ∈ range(extrema(first, blines)...)
    # end
    Iterators.product(eachindex(positions), eachindex(positions)) |>
        Filter(splat(<)) |>
        Map(x->getindex.(Ref(positions), x)) |>
        # Filter() do (p1, p2)
        #     # this is still wrong
        #     xr = min(p1.x, p2.x):max(p1.x, p2.x)
        #     yr = min(p1.y, p2.y):max(p1.y, p2.y)
        #     _bl = vlines |> Filter(r->issubset(yr, last(r))) |> Map(first) |> collect
        #     # println(p1, ", ", p2)
        #     # println("    ", xr, ", ", yr, ", ", extrema(_bl))
        #     # println("      ", issubset(xr, range(extrema(_bl)...)))
        #     !isempty(_bl) && issubset(xr, range(extrema(_bl)...))
        # end |>
        # Filter() do (p1, p2)
        #     # ~5x slower version of below, doesn't bother with early termination from corners
        #     @view(allowed_x[min(p1.y, p2.y):max(p1.y, p2.y) .- miny .+ 1]) |> Map() do xr
        #         issubset(min(p1.x, p2.x):max(p1.x, p2.x), xr)
        #     end  |> ReduceIf(!identity) |> foldxl(&; init=true)
        # end |>
        Filter() do (p1, p2)
            s = sign(mapreduce(-, *, p2, p1))
            iszero(s) || begin
                if s > 0
                    tr = (;x=max(p1.x, p2.x), y=min(p1.y, p2.y))
                    bl = (;x=min(p1.x, p2.x), y=max(p1.y, p2.y))
                    allowedpt(tr) && allowedpt(bl) && begin
                        @view(allowed_x[(tr.y):(bl.y) .- miny .+ 1]) |> Map() do xr
                            issubset((bl.x):(tr.x), xr)
                        end |> foldxl(&; init=true)
                    end
                else
                    br = (;x=max(p1.x, p2.x), y=max(p1.y, p2.y))
                    tl = (;x=min(p1.x, p2.x), y=min(p1.y, p2.y))
                    allowedpt(br) && allowedpt(tl) && begin
                        @view(allowed_x[(tl.y):(br.y) .- miny .+ 1]) |> Map() do xr
                            issubset((tl.x):(br.x), xr)
                        end |> foldxl(&; init=true)
                    end
                end
            end
        end |>
        Map() do (p1, p2)
            prod(zip(p1, p2)) do (c1, c2)
                abs(c2 - c1) + 1
            end
        end |>
            foldxl(max; init=0)

end
@assert main_2("d09/sample") |> isequal(24)
# @assert main_2("d09/input") |> isequal(1501292304)
#@time main_2("d09/input")

end

# Local Variables:
# julia-snail-repl-buffer: "*julia-aoc*"
# julia-snail-port: 10012
# End:
