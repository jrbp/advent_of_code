module AOC25_10_2
using Transducers
using SparseArrays
using LinearAlgebra
using NormalForms: hnfr
using Combinatorics: multiexponents
using LLLplus: integerfeasibility

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

function goalafterpress(goal::AbstractVector{<:Bool}, schematic, npresses)
    (!iszero).(mod.(schematic * npresses, 2)) == goal
end

function goalafterpress(goal::AbstractVector{<:Int}, schematic, npresses)
    res = (schematic * npresses) == goal
    res && println("found soln at: ", npresses)
    res
end

# onehot(n, l) = map(==(n), Base.OneTo(l))

# nwc(k, n) = binomial(n+k-1, k-1)
# #sum(k->nwc(4, k), 1:50) # -> 316250
# nwc(2,1)
# iterator over ways to put n balls in k bins
# these are the 'weak compositions' of n of fixed size k
# tot number is binomial(n+k-1, k-1)
# maps to all subsets of size k-1 formed from set of size n+k-1
# 1,n -> [[n,]]
# k,1 -> map(x->onehot(x, k), 1:k)
# 2,2 -> [[2,0], [1,1], [0,2]]
# 3,2 -> [[2,0,0], [1,1,0], [1,0,1], [0,1,1], [0,0,2]]
#          **||     *|*|     *||*     |*|*     ||**
# just using Combinatorics.multiexponents

function minpress(goal, schematic; maxcheck=20)
    nbtns = size(schematic, 2)
    for npress in Base.OneTo(maxcheck)
        (multiexponents(nbtns, npress) |>
         Map(x -> goalafterpress(goal, schematic, x)) |>
         ReduceIf(identity) |> foldxl(|; init=false)) && return npress
    end
    error("takes more than $maxcheck")
end
function minpress_pat(goal, schematic; maxcheck=30)
    nbtns = size(schematic, 2)
    for npress in Base.OneTo(maxcheck)
        pat = (multiexponents(nbtns, npress) |>
         Map(x -> schematic * x) |>
         ReduceIf(==(goal)) |> foldxl(right))
        pat == goal && return pat
    end
    error("takes more than $maxcheck")
end

function main_1(filename="d10/sample")
    eachline(filename) |> Map(parse_line) |>
    Map(((g, s, _),) -> minpress(g, s)) |> foldxl(+)
end
@assert main_1("d10/sample") |> isequal(7)
@assert main_1("d10/input") |> isequal(459)

nzerorows(mat) = isempty(mat) ? 0 : begin
    Iterators.reverse(eachrow(mat)) |> Map(x->all(iszero, x)) |> ReduceIf(!identity) |> foldxl(+)
end

function nozerorows(G, H)
    zrows = BitVector(map(x->all(iszero, x), eachrow(H)))
    nzrows = map(x->!x, zrows)
    @assert all(iszero, @view H[zrows, :])
    @assert all(iszero, @view G[zrows])
    Hr = @view H[nzrows, :]
    Gr = @view G[nzrows]
    Gr, Hr, zrows
end

isnatural(x; tol=1e-8) = (x >= 0) && (abs(round(x) - x) < tol)
isnatural(x::AbstractArray) = all(isnatural, x)
naturalize(x) = isnatural(x) ? Int.(x) : error("x is not natural: x=$x")

function soln_xf(G::AbstractVector{Int}, H)
    let G=G, H=H
        nr, nc = size(H)
        @assert nc >= nr
        Hl, Hr = (@view(H[:, begin:nr]), @view(H[:, (nr+1):end]))
        Hlinv = inv(Hl)
        HlinvG = Hlinv * G
        HlinvHr = Hlinv * Hr
        (MapCat(n->multiexponents(size(Hr, 2), n)) ⨟
         Map(xr->vcat(HlinvG - HlinvHr * xr, xr)) ⨟
         Filter(isnatural) #⨟ Map(x->Int.(x))
        )
    end
end

function smallernat(l, r)
    !isnatural(l) ? r : begin
        sl, sr = sum(l), sum(r)
        sl < sr ? l : r
    end
end

function wide_hnf_solve_int_mateqn(G::AbstractVector{Int},H)
    nr, nc = size(H)
    @assert nc >= nr
    Hl, Hr = (@view(H[:, begin:nr]), @view(H[:, (nr+1):end]))
    nzhl = nzerorows(Hl)
    if !iszero(nzhl)
        # block form
        #   Htl  Htr
        #   0    Hbr
        trange = 1:(nr-nzhl)
        brange = (nr-nzhl + 1):nr
        @assert all(iszero, @view(Hl[brange, :]))
        Htl = @view(Hl[trange, :])
        Htr = @view(Hr[trange, :])
        Hbr = @view(Hr[brange, :])
        Gt = @view(G[trange])
        Gb = @view(G[brange])
        #xr0 = hnf_solve_int_mateqn(Gb, Hbr)
        #Main.@infiltrate !isnatural(xr0)
        # nGl = Gt .- (Htr * xr0)
        # xl = wide_hnf_solve_int_mateqn(nGl, Htl)
        # x0 = vcat(xl, xr0)
        #ns = Int(sum(xr0))
        #ns = Int(sum(xr0))
        xres = 0:(2 * sum(abs, G)) |> soln_xf(Gb, Hbr) |> Map() do xr
            (xr, Gt .- (Htr * xr))
        end |>
                Filter(x->all(c->abs(round(c) - c)<1e-8, last(x))) |> # filter noninteger G
                Map() do (x, G)
                    x, Int.(G)
                end |>
                Map() do (xr, nG)
                    xl = wide_hnf_solve_int_mateqn(nG, Htl)
                    vcat(xl, xr)
                end |>
                Filter(isnatural) |>
                foldxl(smallernat; init=-ones(Int, size(H, 2)))
        isnatural(xres) || error("failed to find positive integer solution for xr")
        xres
        # res = foldxl(soln_xf(G, H)'(smallernat),
        #     0:maxcheck; init=-ones(Int, size(H, 2)))
    else
        maxcheck = 2 * sum(abs, G)
        #Main.@infiltrate !isnatural(G)
        res = foldxl(soln_xf(G, H)'(smallernat),
            0:maxcheck; init=-ones(Int, size(H, 2)))
        isnatural(res) || error("search failed to find positive integer solution with $maxcheck, G=$G")
        res
        #Main.@infiltrate
        #-ones(Int, size(H, 2))
    end
end

function hnf_solve_int_mateqn(G, H)
    nz = nzerorows(H)
    if iszero(nz) && ==(size(H)...)
        inv(H) * G
        #isnatural(res) ? Int.(res) : throw(DomainError("G, H does not have natural solution"))
    elseif !iszero(nz)
        Gr, Hr, zrows = nozerorows(G, H)
        hnf_solve_int_mateqn(Gr, Hr)
    else
        wide_hnf_solve_int_mateqn(G, H)
    end
end

function solve_int_mateqn(g, a)
    # want integer v which solves
    # g = a * v
    hform = hnfr(a)
    # obtain U which moves problem to form
    # U * g = G = H * v
    # where H = U * a, and H is upper triangular
    G, H = hform.U * g, hform.H 
    @assert !iszero(H[1,1]) && all(iszero, @view H[(begin + 1):end, begin]) "pivot is not main diag"
    soln = hnf_solve_int_mateqn(G, H) |> naturalize
    #soln = isnatural(_soln) ? Int.(soln) : error("could not find natural solution")
    #Main.@infiltrate !(g == a * soln)
    @assert g == a * soln
    soln
end

function minpress_jolt(goal::Vector{Int}, schematic)
    soln = solve_int_mateqn(goal, Matrix{Int}(schematic))
    sum(soln)
end

function main_2(filename="d10/sample")
    res = eachline(filename) |> Map(parse_line) |>
          Map(((_, s, g),) -> minpress_jolt(g, s)) |> foldxl(+)
    res

end

macro comment(x...) end
#@assert main_2("d10/sample") |> isequal(33)
#main_2("d10/input")
#println(main_2("d10/input"))

@comment begin

let machines = parse_line.((eachline("d10/input")))
    for (ii, (lg, as, g)) in enumerate(machines)
        res = solve_int_mateqn(g, Matrix{Int}(as))
        println(ii, ", ", res)
    end
end
#allinps = parse_line.(eachline(("d10/sample")))
# #allinps = parse_line.(eachline(("d10/input")))
# # # #(lg, as, g) = allinps[6] # square mat, but with zero bottom row
# # # (lg, as, g) = allinps[7] # 2D surface of solutions
# # #(lg, as, g) = allinps[9]
# #(lg, as, g) = allinps[12]
# (lg, as, g) = parse_line.((eachline("d10/sample")))[3]
# a = Matrix{Int}(as)
# hform = hnfr(a)
# h = hform.H
# x = integerfeasibility(h[1:4, :], g[1:4])
# hform.U
# hform.U * [x..., 0,0]
# g
#solve_int_mateqn(g, a)
# minpress(g, a)
# solve_int_mateqn(g, a)
end

end

# Local Variables:
# julia-snail-repl-buffer: "*julia-aoc*"
# julia-snail-port: 10012
# End:
