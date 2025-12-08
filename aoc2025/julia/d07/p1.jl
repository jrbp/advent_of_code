module AOC25_07
using Transducers

function main_1(filename="d07/sample")
    sourceline, otherlines = Iterators.peel(eachline(filename))
    initbeams = BitVector(map(==('S'), collect(sourceline)))
    println(String(map(x->x ? '|' : ' ', initbeams)))
    splits, _ = otherlines |> Map(collect) |> Map(x -> BitVector(map(==('^'), x))) |>
        foldxl(; init=(0, initbeams)) do (count, beams), r
        splitterhits = r .& beams
        newcount = count + sum(splitterhits)
        newbeams = (beams .& xor.(beams, r)) .| <<(splitterhits, 1) .| >>(splitterhits, 1) # allocates
        #println(newbeams)
        println(String(map(x->x ? '|' : ' ', newbeams)))
        newcount, newbeams
    end
    splits
end
#@assert main_1("d07/sample") |> isequal(21)
#main_1("d07/input") # -> 1640

function main_2(filename="d07/sample")
    ASCII0 = Int(first("0")) # 48
    sourceline, otherlines = Iterators.peel(eachline(filename))
    initbeams = map(Int ∘ ==('S'), collect(sourceline))
    #println(initbeams)
    splits = otherlines |> Map(collect) |> Map(x -> BitVector(map(==('^'), x))) |>
        foldxl(; init=initbeams) do beammul, splitters
            splitbeams = splitters .* beammul
            newbeams = (beammul .- splitbeams) .+ circshift(splitbeams, 1) .+ circshift(splitbeams, -1)
            #println(newbeams)
            newbeams
    end |> sum
end
@assert main_2("d07/sample") |> isequal(40)
main_2("d07/input") # -> 40999072541589
end
#@btime AOC25_07.main_1("d07/input")

# Local Variables:
# julia-snail-repl-buffer: "*julia-aoc*"
# julia-snail-port: 10012
# End:
