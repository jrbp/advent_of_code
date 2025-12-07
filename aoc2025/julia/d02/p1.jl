module AOC25_02

using Transducers
using Primes

struct CommaSeperatedData
    filepath::String
end
Base.IteratorSize(::Type{<:CommaSeperatedData}) = Base.SizeUnknown()
Base.isdone(csd::CommaSeperatedData, state=open(csd.filepath)) = eof(state)
function Base.iterate(csd::CommaSeperatedData, state=open(csd.filepath))
    !eof(state) ? (chomp(readuntil(state, ',')), state) : begin
        close(state)
        nothing
    end
end

function isvalid_id(id::Int, sdigits=Vector{Int}(undef, ndigits(id)))
    n = ndigits(id)
    nd2 = n ÷ 2
    digits!(sdigits, id)
    mapreduce(!==, |, # does this return early on a true?
        @view(sdigits[begin:nd2]), @view(sdigits[(nd2 + 1):n]))
end

function split_by_numdigits(x::UnitRange{Int})
    l, h = extrema(x) 
    nl, nh = ndigits(l), ndigits(h)
    nl:nh |> Map() do _n
    #map(nl:nh) do _n
        _l = _n == nl ? l : 10 ^ (_n-1)
        _r = _n == nh ? h : (10 ^ (_n) - 1)
        _l:_r
    end
end

function allvalid(x::UnitRange)
    nl = ndigits(minimum(x))
    isodd(nl) && nl == ndigits(maximum(x)) 
end

function main_1(filename)
    CommaSeperatedData(filename) |> collect |>
    Map(x -> split(x, '-'; limit=2)) |> 
    Map(x -> parse.(Int, x)) |> 
    Map(splat(range)) |> 
    MapCat(split_by_numdigits) |> # e.g. 20:1500 -> [20:99, 100:999, 1000:1500]
    Filter(!allvalid) |>          # e.g. skip 100:999 since all odd ndigits are valid
    Cat() |>
    Filter(!isvalid_id) |>
    #collect
    foldxt(+)
end

function isinvalid_id_two(id::Int, sdigits=Vector{Int}(undef, ndigits(id)))
    digits!(sdigits, id)
    n = ndigits(id)
    1:(n ÷ 2) |> Filter(sl->iszero(rem(n, sl))) |> Map() do seqlength
        seq = view(sdigits, 1:seqlength)
        repeats = n ÷ seqlength
        2:repeats |> Map() do r
            view(sdigits, (1 + (r-1)*seqlength):(r*seqlength))
        end |> Map(==(seq)) |> foldxl(&)
    end |> foldxl(|; init=false)
end

function main_2(filename)
    CommaSeperatedData(filename) |> collect |>
    Map(x -> split(x, '-'; limit=2)) |> 
    Map(x -> parse.(Int, x)) |> 
    Map(splat(range)) |> 
    #MapCat(split_by_numdigits) |> # e.g. 20:1500 -> [20:99, 100:999, 1000:1500]
    Cat() |>
    Filter(isinvalid_id_two) |>
    foldxt(+)
end
@assert main_2("d02/sample") |> isequal(4174379265)
@assert main_2("d02/input") |> isequal(19058204438)

# function ismultipleofany(n, fs)
#     foldxl(|, Map(x->n % x == 0) ⨟ ReduceIf(identity), fs; init = false)
# end
# function sieve(sofar, x)
#     ismultipleofany(x, sofar) ? (0, sofar) : (x, push!!(sofar, x))
# end
# function primesupto(n)
#     prime_xf = ScanEmit(sieve, Int[]) ⨟ Filter(!iszero)
#     foldxl(push!!, prime_xf, 2:n; init=Int[])
# end

@inline function mobius(n)
    eachfactor(n) |>
    Map(last) |>
    ReduceIf(p->p>1) |>
    Map() do x
        x > 1 ? 0 : -1
    end |>
    foldxl(*; init=1)
end

divisorof(a) = x->iszero(mod(a, x))
function _sum_invalid_ids_of_length(p)
    p == 1 ? 0 :
    2:p |> Filter(divisorof(p)) #|> Filter(!iszero ∘ mobius) |>
    Map() do r
        q = p ÷ r
        pat2n = (10^p - 1) ÷ (1 - 10^q)
        patterns = (10^(q-1)):(10^q - 1)
        mobius(r) * sum(patterns) * pat2n
    end |> foldxl(+)
end

@inline firstndigits(x, n) = x ÷ 10^(ndigits(x) - n) 
@inline lastndigits(x, n) = x % 10^(n) 

function sum_invalid_ids(rr::UnitRange)
    p = ndigits(first(rr)) # must equal ndigits(last(r))
    minrr = minimum(rr)
    maxrr = maximum(rr)
    p == 1 ? 0 : 
    2:p |> Filter(divisorof(p)) |> Filter(!iszero ∘ mobius) |>
    Map() do r
        q = p ÷ r
        pat2n = (10^p - 1) ÷ (10^q - 1)
        patmin = let _s = firstndigits(minrr, q)
            _s + ((_s * pat2n) < minrr)
        end
        patmax = let _s = firstndigits(maxrr, q)
            _s - ((_s * pat2n) > maxrr)
        end
        patterns = patmin:patmax
        -1 * mobius(r) * sum(patterns) * pat2n
    end |> foldxl(+)
end
_fl(x) = first(x), last(x)
function main_2fast(filename)
    CommaSeperatedData(filename) |>
    Map(x -> _fl(split(x, '-'; limit=2))) |> 
    Map(x -> map(Base.Fix1(parse, Int), x)) |> 
    collect |>
    Map(splat(range)) |> 
    MapCat(split_by_numdigits) |> # e.g. 20:1500 -> [20:99, 100:999, 1000:1500]
    Map(sum_invalid_ids) |>
    foldxl(+)
end
@assert main_2fast("d02/sample") |> isequal(4174379265)
@assert main_2fast("d02/input") |> isequal(19058204438)
end

#@btime AOC25_02.main_2("d02/input")
@btime AOC25_02.main_2fast("d02/input")
