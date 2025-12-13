module AOC25_10
using Transducers
using SparseArrays
using LinearAlgebra
using Combinatorics: multiexponents

macro comment(x...) end

function tosparsemat(btns, nlights)
    is, js = Int[], Int[]
    for b in eachindex(btns)
        for c in eachindex(btns[b])
            push!(is, btns[b][c] + 1)
            push!(js, b)
        end
    end
    nconnect, nbtns = sum(length, btns), length(btns)
    sparse(is, js, ones(Bool, nconnect), nlights, nbtns)
end

function parse_line(ln)
    parts = map(x -> x[begin+1:end-1], split(ln))
    goal = map(==('#'), collect(first(parts)))
    others = map(x -> parse.(Int, split(x, ',')), parts[begin+1:end])
    joltreq = others[end]
    @assert length(goal) == length(joltreq)
    schem = tosparsemat(others[begin:end-1], length(goal))
    goal, schem, joltreq
end

function goalafterpress(goal::AbstractVector{Bool}, schematic, npresses)
    (!iszero).(mod.(schematic * npresses, 2)) == goal
end

# function soln_xf(g::AbstractVector{Bool}, s)
#     nbtns = size(s, 2)
#     MapCat(npress->multiexponents(nbtns, npress)) ⨟
#      Filter(x -> all(y-> y ∈ (0,1) , x)) ⨟
#      Filter(x -> goalafterpress(g, s, x))
# end
function soln_itr(g::AbstractVector{Bool}, s)
    nbtns = size(s, 2)
    0:nbtns |> MapCat(npress->multiexponents(nbtns, npress)) |>
     Filter(x -> all(y-> y ∈ (0,1) , x)) |>
     Filter(x -> goalafterpress(g, s, x))
end

function minpress(goal, schematic)
    res = soln_itr(goal, schematic) |>
          Take(1) |>
          collect
    isempty(res) ? error("takes more than") : first(res)
end

function main_1(filename="d10/sample")
    eachline(filename) |> collect |>
    Map(parse_line) |>
    Map(((g, s, _),) -> minpress(g, s)) |>
    Map(sum) |>
    foldxt(+)
end
@assert main_1("d10/sample") |> isequal(7)
@assert main_1("d10/input") |> isequal(459)

# from brute force attempt
# function goalafterpress(goal::AbstractVector{<:Int}, schematic, npresses)
#     schematic * npresses == goal
# end

function nminpress(goal::AbstractVector{Int}, schematic,
    cache=Dict((zeros(Int, length(goal))=>0,))
)
    get!(cache, goal) do
        soln_itr((@. Bool(mod(goal, 2))), schematic) |>
        Map(xev->(sum(xev), (goal - (schematic * xev)) .÷ 2)) |>
        Filter(xg->all(i->i>=0, last(xg))) |>
        Map() do (sxev, newg)
            2 * nminpress(newg, schematic, cache) + sxev
        end |>
        foldxl(min; init=2^60)
    end
end

function main_2(filename="d10/sample")
    eachline(filename) |> collect |>
    Map(parse_line) |>
    Map(((_, s, g),) -> nminpress(g,s)) |>
    foldxt(+; init=0)
end

@assert main_2("d10/sample") |> isequal(33)
@assert main_2("d10/input") |> isequal(18687) # slow

end

# Local Variables:
# julia-snail-repl-buffer: "*julia-aoc*"
# julia-snail-port: 10012
# End:
