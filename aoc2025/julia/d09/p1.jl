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
function merge_nondisjoint(l, r)
    lstart, lstop = extrema(l)
    rstart, rstop = extrema(r)
    range(min(lstart, rstart), max(lstop, rstop))
end
function merge_ranges(rin)
    rs = sort(rin; by=first)
    maxi = length(rs)
    res = UnitRange{Int}[]
    i = 1
    while i <= maxi
        thisr = rs[i]
        nexti = i + 1
        while ((nexti <= maxi) && !isdisjoint(thisr, rs[nexti]))
            thisr = merge_nondisjoint(thisr, rs[nexti])
            nexti += 1
        end
        i = nexti
        push!(res, thisr)
    end
    res
end

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
    Iterators.product(eachindex(positions), eachindex(positions)) |>
        Filter(splat(<)) |>
        Map(x->getindex.(Ref(positions), x)) |>
        # Filter() do (p1, p2)
        #     # currently as slow as  version below
        #     # could be best if we could avoid the collect + call
        #     # to merge_ranges and have vlines fold to left_ok/right_ok
        #     xr = min(p1.x, p2.x):max(p1.x, p2.x)
        #     yr = min(p1.y, p2.y):max(p1.y, p2.y)
        #     left_vl = vlines |> Filter(v->first(v) <= minimum(xr)) |> Map(last) |> collect
        #     isempty(left_vl) && return false
        #     left_ok = merge_ranges(left_vl) |> Map(x->issubset(yr, x)) |> ReduceIf(identity) |> foldxl(|)
        #     left_ok || return false
        #     right_vl = vlines |> Filter(v->first(v) >= maximum(xr)) |> Map(last) |> collect
        #     isempty(right_vl) && return false
        #     right_ok = merge_ranges(right_vl) |> Map(x->issubset(yr, x)) |> ReduceIf(identity) |> foldxl(|)
        #     return right_ok
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
@assert main_2("d09/input") |> isequal(1501292304)
# main_2("d09/sample") |> isequal(24) |> println
# tt = Threads.@spawn @time main_2("d09/input")
# fetch(tt)
#println(main_2("d09/input") == 1501292304)
@time main_2("d09/input")

end

# Local Variables:
# julia-snail-repl-buffer: "*julia-aoc*"
# julia-snail-port: 10012
# End:
