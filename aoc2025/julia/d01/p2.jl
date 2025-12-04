const safeknob_length = 100
parsetoshift(l::String) = let (d, s) = iterate(l)
    d == '#' ? 0 : let clicks = parse(Int, Base.rest(l, s))
        d == 'L' ? -clicks : 
        d == 'R' ? clicks :
        error("unknown direction $(d)")
    end
    #shift = d == 'L' ? -clicks : d == 'R' ? clicks : error("unknown direction $(d)")
end

struct SafeKnobInt
    pos::Int
    function SafeKnobInt(pos = 50)
        new(mod(pos, safeknob_length))
    end
end
function turn_count0s(sk::SafeKnobInt, clicks::Int)
    startszero = iszero(sk.pos)
    naivesum = sk.pos + clicks
    newsk = SafeKnobInt(naivesum)
    wind = ((clicks < 0) * (iszero(newsk.pos) - iszero(sk.pos))
            + abs(newsk.pos - naivesum) ÷ safeknob_length)
    #println(clicks, ", ", newsk.pos, ", ", wind)
    newsk, wind
end

function mainInt(file)
    sk = SafeKnobInt()
    res = 0
    for l in eachline(file)
        sk, n0s  = turn_count0s(sk, parsetoshift(l))
        res += n0s
    end
    res
end
