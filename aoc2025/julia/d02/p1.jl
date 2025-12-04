module AOC25_02

using Transducers

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
    #collect
    foldxt(+)
end

end
