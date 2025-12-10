module AOC25_10
using Transducers
using SparseArrays
using LinearAlgebra
using NormalForms: hnfr

#The Combinations iterator # from Combinatorics.jl
struct Combinations
    n::Int
    t::Int
end

@inline function Base.iterate(c::Combinations, s = [min(c.t - 1, i) for i in 1:c.t])
    if c.t == 0 # special case to generate 1 result for t==0
        isempty(s) && return (s, [1])
        return
    end
    for i in c.t:-1:1
        s[i] += 1
        if s[i] > (c.n - (c.t - i))
            continue
        end
        for j in i+1:c.t
            s[j] = s[j-1] + 1
        end
        break
    end
    s[1] > c.n - c.t + 1 && return
    (s, s)
end

Base.length(c::Combinations) = binomial(c.n, c.t)

Base.eltype(::Type{Combinations}) = Vector{Int}

"""
    combinations(a, n)

Generate all combinations of `n` elements from an indexable object `a`. Because the number
of combinations can be very large, this function returns an iterator object.
Use `collect(combinations(a, n))` to get an array of all combinations.
"""
function combinations(a, t::Integer)
    if t < 0
        # generate 0 combinations for negative argument
        t = length(a) + 1
    end
    reorder(c) = [a[ci] for ci in c]
    (reorder(c) for c in Combinations(length(a), t))
end


"""
    combinations(a)

Generate combinations of the elements of `a` of all orders. Chaining of order iterators
is eager, but the sequence at each order is lazy.
"""
combinations(a) = Iterators.flatten([combinations(a, k) for k = 0:length(a)])

struct MultiExponents{T}
    c::T
    nterms::Int
end

# Standard stars and bars:
# https://en.wikipedia.org/wiki/Stars_and_bars_(combinatorics)
function Base.iterate(m::MultiExponents, s = nothing)
    next = s === nothing ? iterate(m.c) : iterate(m.c, s)
    next === nothing && return
    stars, ss = next

    # stars minus their consecutive
    # position becomes their index
    result = zeros(Int, m.nterms)
    for (i, s) in enumerate(stars)
        result[s-i+1] += 1
    end

    result, ss
end

Base.length(m::MultiExponents) = length(m.c)
Base.eltype(::Type{MultiExponents{T}}) where {T} = Vector{Int}

"""
    multiexponents(m, n)

Returns the exponents in the multinomial expansion (x₁ + x₂ + ... + xₘ)ⁿ.

For example, the expansion (x₁ + x₂ + x₃)² = x₁² + x₁x₂ + x₁x₃ + ...
has the exponents:

```julia-repl
julia> collect(multiexponents(3, 2))
6-element Vector{Vector{Int64}}:
 [2, 0, 0]
 [1, 1, 0]
 [1, 0, 1]
 [0, 2, 0]
 [0, 1, 1]
 [0, 0, 2]
```
"""
function multiexponents(m, n)
    # number of stars and bars = m+n-1
    c = combinations(1:m+n-1, n)

    MultiExponents(c, m)
end

function tosparsemat(btns, nlights)
    is, js = Int[], Int[]
    for b in eachindex(btns)
        for c in eachindex(btns[b])
            push!(is, btns[b][c]+1)
            push!(js, b)
        end
    end
    nconnect, nbtns = sum(length, btns), length(btns)
    sparse(is, js, ones(Bool, nconnect), nlights, nbtns)
end

function parse_line(ln)
    parts = map(x->x[begin+1:end-1], split(ln))
    goal = map(==('#'), collect(first(parts)))
    others = map(x->parse.(Int, split(x, ',')), parts[begin+1:end])
    joltreq = others[end]
    @assert length(goal) == length(joltreq)
    schem = tosparsemat(others[begin:end-1], length(goal))
    goal, schem, joltreq
end

function goalafterpress(goal::Vector{Bool}, schematic, npresses)
    (!iszero).(mod.(schematic * npresses, 2)) == goal
end

function goalafterpress(goal::Vector{Int}, schematic, npresses)
    res = (schematic * npresses) == goal
    res && println("found soln at: ", npresses)
    res
end

# onehot(n, l) = map(==(n), Base.OneTo(l))

# nwc(k, n) = binomial(n+k-1, k-1)
# #sum(k->nwc(4, k), 1:50) # -> 316250
# nwc(2,1)

function weak_compositions(k, n)
    multiexponents(k, n)
    # iterator over ways to put n balls in k bins
    # these are the 'weak compositions' of n of fixed size k
    # tot number is binomial(n+k-1, k-1)
    # maps to all subsets of size k-1 formed from set of size n+k-1
    # 1,n -> [[n,]]
    # k,1 -> map(x->onehot(x, k), 1:k)
    # 2,2 -> [[2,0], [1,1], [0,2]]
    # 3,2 -> [[2,0,0], [1,1,0], [1,0,1], [0,1,1], [0,0,2]]
    #          **||     *|*|     *||*     |*|*     ||**
end

function minpress(goal, schematic; maxcheck = 20)
    nbtns = size(schematic, 2)
    for npress in Base.OneTo(maxcheck)
        (weak_compositions(nbtns, npress) |>
        Map(x-> goalafterpress(goal, schematic, x)) |>
        ReduceIf(identity) |> foldxl(|; init=false)) && return npress
    end
    error("takes more than $maxcheck")
end
# function minpress_maybe(goal::Vector{Int}, schematic; maxcheck = 100)
#     nbtns = size(schematic, 2)
#     u, s, v = svd(schematic)
# end

function main_1(filename="d10/sample")
    eachline(filename) |> Map(parse_line) |>
        Map(((g, s, _),) -> minpress(g, s)) |> foldxl(+)
end
@assert main_1("d10/sample") |> isequal(7)
@assert main_1("d10/input") |> isequal(459)

nzerorows(mat) = sum(x->all(iszero, x), eachrow(mat))

function _search_int_mateqn(Gr, Ht, Hb; maxtry=10)
    # search over vr = N^(size(Hb, 2))
    # for set of integer vl
    # which solve Gr = Ht * vl + Hb * vr
    Htinv = inv(Ht)
    for n in 0:maxtry
        # weak_compositions(sizeof(Hb, 2), n) |>
        #     Map(vr-> Hinv * (Gr - Hb * vr))
        for vr in weak_compositions(sizeof(Hb, 2), n)
            vl = Htinv * (Gr - Hb * vr)
            # do we need to valiedate it?
            error("not finished")
        end
    end
end

function _nz_solve_int_mateqn(Gr, Hr)
    nr, nc = size(Hr)
    @assert nc >= nr
    if nr == nc
        Int.(inv(Hr) * Gr)
    else
        Ht, Hb = @view(Hr[:, begin:nr]), @view(Hr[:, (nr+1):end])
        #_search_int_mateqn(Gr, Ht, Hb)
        [Int.(inv(Ht) * Gr)..., zeros(Int, size(Hb, 2))...]
    end
end

function solve_int_mateqn(g, a)
    # want integer v which solves
    # g = a * v
    hform = hnfr(a)
    # obtain U which moves problem to form
    # U * g = G = H * v
    # where H = U * a, and H is upper triangular
    H = hform.H
    G = hform.U * g
    nzr = nzerorows(H)
    @assert all(iszero, @view H[(end-nzr+1):end, :])
    @assert all(iszero, G[(end-nzr+1):end])
    Hr = @view H[begin:(end-nzr), :]
    Gr = @view G[begin:(end-nzr)]
    _soln = _nz_solve_int_mateqn(Gr, Hr)
    topad = size(a, 2) - length(_soln)
    soln = [_soln..., zeros(Int, topad)...]
    @assert g == a * soln
    soln
end

function main_2(filename="d10/sample")
    res = eachline(filename) |> Map(parse_line) |>
        Map(((_, s, g),) -> minpress(g, s; maxcheck=58)) |> foldxl(+)
    println("finished with with res = ", res)
    res

end
@assert main_2("d10/sample") |> isequal(33)

# allinps = parse_line.(eachline(("d10/input")))
# #(lg, as, g) = allinps[6] # square mat, but with zero bottom row
# # (lg, as, g) = allinps[7] # 2D surface of solutions
# #(lg, as, g) = allinps[9]
# (lg, as, g) = allinps[9]
# a = Matrix{Int}(as)
# solve_int_mateqn(g, a)

end

# Local Variables:
# julia-snail-repl-buffer: "*julia-aoc*"
# julia-snail-port: 10012
# End:
