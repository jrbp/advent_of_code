module AOC25_06
using Transducers

function readmo(filename)
    ops = eachline(filename) |> last |> split .|> ==("*")
    nums = eachline(filename) |>
           Map(split) |>
           Filter(isdigit ∘ first ∘ collect ∘ first) |>
           Map(x -> (parse.(Int, x))') |>
           foldxl(vcat)
    nums, ops
end

function main_1(filename="d06/sample")
    mat, ops = readmo(filename)
    zip(eachcol(mat), ops) |> Map() do (v, o)
        reduce(o ? (*) : (+), v)
    end |>
    foldxl(+)
end
@assert main_1("d06/sample") |> isequal(4277556)
main_1("d06/input") # -> 3968933219902


function cephread(filename)
    ops = eachline(filename) |> last |> split .|> ==("*")
    numchars = eachline(filename) |>
               DropLast(1) |>
               Map(collect) |>
               Map(permutedims) |>
               foldxl(vcat)
    sepinds = eachcol(numchars) |> Enumerate() |> Filter(((i, x),)->all(x .== ' ')) |> Map(first) |> collect
    blockcols = UnitRange{Int}[]
    let ii = 1
        for si in sepinds
            push!(blockcols, ii:(si-1))
            ii = si+1
        end
        push!(blockcols, ii:(size(numchars, 2)))
    end
    nums = blockcols |> Map(cs -> view(numchars, :, cs)) |> Map() do blk
        eachcol(blk) |> Map() do x
            foldl(*, x)
        end |> Map(x->parse(Int, x)) |> collect
    end |> collect
    zip(nums, ops) |> Map() do (ns, o)
        reduce(o ? (*) : (+), ns)
    end |> foldxl(+)
end

function main_2(filename="d06/sample")
    cephread(filename)
end
@assert main_2("d06/sample") |> isequal(3263827)
@assert main_2("d06/input") |> isequal(6019576291014)

end
