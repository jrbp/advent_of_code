module AOC25_05
using Transducers

_fl(x) = first(x), last(x)
function ranges_avail(filename)
    read_range = Ref(true)
    ranges = UnitRange{Int}[]
    avail = Int[]
    for l in eachline(filename)
        iszero(length(l)) ? read_range[] = false :
        !(read_range[]) ? push!(avail, parse(Int, l)) :
        push!(ranges, range(map(Base.Fix1(parse, Int), _fl(split(l, "-"; limit=2)))...))
    end
    ranges, avail
end

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

function main_1(filename)
    _ranges, avail = ranges_avail(filename)
    ranges = merge_ranges(_ranges)
    avail |> Map() do i
        ranges |>
        Map(r-> i ∈ r) |>
        ReduceIf(identity) |>
        foldxl(|; init=false)
    end |> foldxl(+)
end
@assert main_1("d05/sample") |> isequal(3)
@assert main_1("d05/input") |> isequal(896)


function main_2(filename)
    ranges_avail(filename) |>
    first |>
    merge_ranges |>
    Map(length) |>
    foldxl(+)
end

# using LinearAlgebra: I
# let uniq_rangs = ranges_avail("d05/input") |> first |> merge_ranges
#     @assert map(splat(!isdisjoint), Iterators.product(uniq_rangs, uniq_rangs)) == I
# end

@assert main_2("d05/sample") |> isequal(14)
@assert main_2("d05/input") |> isequal(346240317247002)
end
