module AOC25_08
using Transducers
using StaticArrays

function readpositions(filename)
    eachline(filename) |>
        Map(l->split(l, ",")) |>
        Map(x->parse.(Int, x)) |>
        Map(splat(SVector{3})) |>
        collect
end

function to_circuits(connections)
    foldxl(connections; init=Set{Int}[]) do sofar, c
        sepfrom = isdisjoint.(sofar, Ref(c))
        all(sepfrom) ? [sofar..., Set(c)] : begin
            [union(Set(c), sofar[(!).(sepfrom)]...), sofar[sepfrom]...]
        end
    end
end

function main_1(filename="d08/sample", nconnect=10)
    positions = readpositions(filename)
    Iterators.product(eachindex(positions), eachindex(positions)) |> Filter() do (i1, i2)
        i2 > i1
    end |> Map() do (i1, i2)
        (i1, i2), sum(x->x^2, positions[i1] .- positions[i2])
    end |> collect |> (x->sort(x; by=last)) |> Map(first) |> Take(nconnect) |>
        to_circuits |> Map(length) |> collect |> sort |> TakeLast(3) |> foldxl(*)
end

@assert main_1("d08/sample", 10) |> isequal(40)
@assert main_1("d08/input", 1000) |> isequal(127551)

function to_circuits_earlyterm(connections, stopsize)
    foldxl(connections; init=Set{Int}[]) do sofar, c
        sepfrom = isdisjoint.(sofar, Ref(c))
        all(sepfrom) ? [sofar..., Set(c)] : begin
            next = [union(Set(c), sofar[(!).(sepfrom)]...), sofar[sepfrom]...]
            if (length(next) == 1) && (length(first(next)) == stopsize)
                reduced([Set(c),]) # wrapped in [Set] for type stability
            else
                next
            end
        end
    end
end

function main_2(filename="d08/sample")
    positions = readpositions(filename)
    fi1, fi2 = Iterators.product(eachindex(positions), eachindex(positions)) |>
        Filter() do (i1, i2)
            i2 > i1
        end |> Map() do (i1, i2)
            (i1, i2), sum(x->x^2, positions[i1] .- positions[i2])
        end |> collect |> (x->sort(x; by=last)) |> Map(first) |>
            (x->to_circuits_earlyterm(x, length(positions))) |>
            first
    first(positions[fi1]) * first(positions[fi2])

end
@assert main_2("d08/sample") |> isequal(25272)
main_2("d08/input") # -> 2347225200
end

# Local Variables:
# julia-snail-repl-buffer: "*julia-aoc*"
# julia-snail-port: 10012
# End:
