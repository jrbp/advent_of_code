module AOC25_10_2
using Transducers
using SparseArrays
using LinearAlgebra
using NormalForms: hnfr
using Combinatorics: multiexponents

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

nzerorows(mat) = isempty(mat) ? 0 : sum(x -> all(iszero, x), eachrow(mat))

# function reduce_singular(G, H)
#     nsq = min(size(H)...)
#     nzr = nzerorows(H)
#     iszero(nzr) ? (@view(G[begin:nsq]), @view(H[begin:nsq, begin:nsq])) :
#     # @assert all(iszero, @view H[(end-nzr+1):end, :])
#     # @assert all(iszero, G[(end-nzr+1):end])
#     Hr = @view H[begin:(end-nzr), begin:(end-nzr)]
#     Gr = @view G[begin:(end-nzr)]
#     nzerorows(Hr) > 0 ? reduce_singular(Gr, Hr) : (Gr, Hr)
# end
# function nozerorows(G, H, Hb)
#     nzr = nzerorows(H)
#     #nc = size(H, 1)
#     # Hr = @view H[begin:(nc-nzr), begin:(nc-nzr)]
#     # Gr = @view G[begin:(end-nzr)]
#     Gr = @view G[begin:(end-nzr)]
#     Hr = @view H[begin:(end-nzr), :]
#     Hb2 = @view Hb[begin:(end-nzr), :]
#     Gr, Hr, Hb2
# end
function nozerorows(G, H)
    nzr = nzerorows(H)
    @assert all(iszero, @view H[(end-nzr+1):end, :])
    @assert all(iszero, @view G[(end-nzr+1):end])
    Hr = @view H[begin:(end-nzr), :]
    Gr = @view G[begin:(end-nzr)]
    Gr, Hr
end

# function relG_Hinv(G, H, Hb)
#     nr, nc = size(H)
#     nz = nzerorows(H)
#     if iszero(nz) && (nr == nc)
#         G, H, Hb
#     elseif !iszero(nz)
#         Gr, Hr, Hbb = nozerorows(G, H, Hb)
#         relG_Hinv(Gr, Hr, Hbb)
#     else
#         Ht, Hbb = @view(H[:, begin:nr]), hcat(@view(H[:, (nr+1):end]), Hb)
#         relG_Hinv(G, Ht, Hbb)
#     end
# end

# function _search_int_mateqn(G, Ht, Hb; maxcheck=50)
#     # search over vr = N^(size(Hb, 2))
#     # for set of integer vl
#     # which solve Gr = Ht * vl + Hb * vr
#     #Gr, Htt, Hbb = relG_Hinv(G, Ht, Hb)
#     Hbb = Hb
#     Gr = G
#     u, s, v = svd(Ht)
#     Htinv = v * diagm(map(x->abs(x)<1e-8 ? 1 : 1/x, s)) * u'
#     for n in 0:maxcheck
#         # multiexponents(sizeof(Hb, 2), n) |>
#         #     Map(vr-> Hinv * (Gr - Hb * vr))
#         for vr in multiexponents(size(Hbb, 2), n)
#             vl = Htinv * (Gr - Hbb * vr)
#             all(x->(x >= 0) && (abs(round(x) - x)<1e-8), vl) && return round.(Int,[vl..., vr...])
#         end
#     end
#     display(hcat(Ht, Hb))
#     display(G)
#     abs.(G)
# end

function wide_hnf_solve_int_mateqn(G,H; maxcheck=20)
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
        xr = hnf_solve_int_mateqn(Gb, Hbr)
        xl = wide_hnf_solve_int_mateqn(Gt .- (Htr * xr), Htl) # missing options? # a mistake?
        vcat(xl, xr)
    else
        Hlinv = inv(Hl)
        ## a = Hlinv * G
        # B = Hlinv * Hr
        ugh = Vector{Int}[]
        for n in 0:maxcheck
            for xr in multiexponents(size(Hr, 2), n)
                xl = hnf_solve_int_mateqn(Hlinv * (G - Hr * xr), Hl)
                @assert all(x-> x<1e-8, abs.(H * vcat(xl, xr) .- G))
                # all(x-> x<1e-8, abs.(H * vcat(xl, xr) .- G)) || begin
                #     println()
                #     println("wat?")
                #     println("G")
                #     display(G)
                #     display(Hl * xl + Hr * xr)
                #     _soln = vcat(xl, xr)
                #     println("x")
                #     display(_soln)
                #     println("H * x")
                #     display(H * _soln)
                # end
                #all(x->(x >= 0) && (abs(round(x) - x)<1e-8), xl) && return round.(Int,vcat(xl, xr))
                all(x->(x >= 0) && (abs(round(x) - x)<1e-8), xl) && push!(ugh, round.(Int,vcat(xl, xr)))
            end
        end
        if !isempty(ugh) 
            nps = sum.(ugh)
            return ugh[argmin(nps)]
        end
        println("could not find solution in $(maxcheck)")
        -ones(Int, size(H, 2))
    end
end

function hnf_solve_int_mateqn(G, H)
    nz = nzerorows(H)
    if iszero(nz) && ==(size(H)...)
        Int.(inv(H) * G)
    elseif !iszero(nz)
        Gr, Hr = nozerorows(G, H)
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
    soln = hnf_solve_int_mateqn(G, H)
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
#
@assert main_2("d10/sample") |> isequal(33)
 main_2("d10/input")
#println(main_2("d10/input"))

@comment begin

let machines = parse_line.((eachline("d10/sample")))
    (lg, as, g) = parse_line(first(eachline("d10/sample")))
    println()
    (lg, as, g) = machines[1]
    println(sum(minpress(g, as)))
    #res = solve_int_mateqn(g, Matrix{Int}(as))
    # println(res)
    # println(sum(res))
    println()
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
