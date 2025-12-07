module AOC25_04
using Transducers

function parse_chararray(f, filename)
    eachline(filename) |>
    Map(collect) |>
    Map(x->f.(x)) |>
    Map(adjoint) |>
    foldxl(vcat)
end

function getindex_withdefault(arr, ind::CartesianIndex, def=0)
    ind ∈ CartesianIndices(arr) ? arr[ind] : def
end

const _stencil = Iterators.product(-1:1, -1:1) |> Filter(!=((0,0))) |> Map(CartesianIndex) |> collect
removable(floor) = 
    CartesianIndices(floor) |>
    Filter(i->floor[i]) |> 
    Map() do i
        _stencil |>
        Map(si->getindex_withdefault(floor, i+si)) |>
        foldxl(+)
    end |>
    Map(x->x < 4) |>
    foldxl(+)

function main_1(filename="d04/sample")
    removable(parse_chararray(==('@'), filename))
end
@assert AOC25_04.main_1("d04/sample") == 13
@assert main_1("d04/input") |> isequal(1445)

removable_inds(floor) = 
    CartesianIndices(floor) |>
    Filter(i->floor[i]) |> 
    Filter() do i
        (_stencil |>
        Map(si->getindex_withdefault(floor, i+si)) |>
        foldxl(+)) < 4
    end

function remove_step!(floor)
    removable_inds(floor) |> collect |>
    Map() do i
        floor[i] = false
        1
    end |> foldxl(+; init=0)
end

function main_2(filename="d04/sample")
    floor = parse_chararray(==('@'), filename)
    totremoved = 0
    removed_some = true
    while removed_some
        rm_thisstep = remove_step!(floor)
        removed_some = !iszero(rm_thisstep)
        totremoved += rm_thisstep
    end
    totremoved
end
@assert main_2("d04/sample") |> isequal(43)
@assert main_2("d04/input") |> isequal(8317)
end

