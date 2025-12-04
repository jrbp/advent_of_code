const safeknob_length = 100
parsetoshift(l::String) = let (d, s) = iterate(l)
    clicks = parse(Int, Base.rest(l, s))
    #shift = d == 'L' ? -clicks : d == 'R' ? clicks : error("unknown direction $(d)")
    d == 'L' ? -clicks : clicks
end

struct SafeKnobBV
    pos::BitVector
    function SafeKnobBV()
        new(BitVector(map(==(51),
            Base.OneTo(safeknob_length))))
    end
end
atzero(sk::SafeKnobBV) = first(sk.pos)
turn!(sk::SafeKnobBV, l::String) = circshift!(sk.pos, parsetoshift(l))
function mainBV(file)
    sk = SafeKnobBV()
    res = 0
    for l in eachline(file)
        turn!(sk, l)
        res += atzero(sk)
    end
    res
end


struct SafeKnobInt
    pos::Int
    function SafeKnobInt(pos = 50)
        new(mod(pos, safeknob_length))
    end
end
turn(sk::SafeKnobInt, clicks::Int) = SafeKnobInt(sk.pos + clicks)
atzero(sk::SafeKnobInt) = iszero(sk.pos)

function mainInt(file)
    # foldl(eachline(file); init = (SafeKnobInt(), 0)) do (sk, count), l
    #     newsk = turn(sk, parsetoshift(l))
    #     (newsk, count + atzero(newsk))
    # end |> last
    sk = SafeKnobInt()
    res = 0
    for l in eachline(file)
        sk = turn(sk, parsetoshift(l))
        res += atzero(sk)
    end
    res
end
