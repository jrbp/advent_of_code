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
        Filter() do (p1, p2)
            delta = map(-, p2, p1)
            s = sign(prod(delta))
            @assert allowedpt(p1)
            @assert allowedpt(p2)
            iszero(s) || begin
                if s > 0
                    tl, br = p1.x < p2.x ? (p1, p2) : (p2, p1)
                    tr = (;x=br.x, y=tl.y)
                    bl = (;x=tl.x, y=br.y)
                    @assert tr == (;x=max(p1.x, p2.x), y=min(p1.y, p2.y))
                    @assert bl == (;x=min(p1.x, p2.x), y=max(p1.y, p2.y))
                    allowedpt(tr) && allowedpt(bl) && begin
                        allowed_x[(tr.y):(bl.y) .- miny .+ 1] |> Map() do xr
                            issubset((bl.x):(tr.x), xr)
                        end |> foldxl(&; init=true)
                    end
                else
                    tr, bl = p1.x > p2.x ? (p1, p2) : (p2, p1)
                    br = (;x=tr.x, y=bl.y)
                    tl = (;x=bl.x, y=tr.y)
                    @assert br == (;x=max(p1.x, p2.x), y=max(p1.y, p2.y)) "s=$(s), p1=$(p1), p2=$(p2), br=$(br)"
                    @assert tl == (;x=min(p1.x, p2.x), y=min(p1.y, p2.y))
                    allowedpt(br) && allowedpt(tl) && begin
                        allowed_x[(tl.y):(br.y) .- miny .+ 1] |> Map() do xr
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
#@assert main_2("d09/input") |> !isequal(4613357005)
#main_2("d09/input") # -> 1501292304

end

# Local Variables:
# julia-snail-repl-buffer: "*julia-aoc*"
# julia-snail-port: 10012
# End:
