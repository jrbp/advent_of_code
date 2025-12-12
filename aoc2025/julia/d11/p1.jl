module AOC25_11

macro comment(e...) end

function parse_connections(filename)
    Dict(map(eachline(filename)) do l
        s_dev, s_outputs = split(l, ':')
        dev = Symbol(s_dev)
        outputs = map(Symbol, split(s_outputs))
        dev => outputs
    end)
end

function npaths_to_out(node, devouts, cache=Dict{Symbol, Int}((:out=>1,)))
    get!(cache, node) do
        sum(n->npaths_to_out(n, devouts, cache), devouts[node]; init=0)
    end
end

function main_1(filename="d11/sample")
    devouts = parse_connections(filename)
    npaths_to_out(:you, devouts)
end
@assert main_1("d11/sample") |> isequal(5)
@assert main_1("d11/input") |> isequal(791)
@comment begin
    Main.@btime main_1("d11/input") #773.099 μs (6411 allocations: 527.28 KiB)
end

function npaths_between(devouts, fnode, lnode)
    # kinda hacky naming now, but w/e
    npaths_to_out(fnode, devouts, Dict{Symbol, Int}((lnode=>1,)))
end

function main_2(filename="d11/sample2")
    devouts = parse_connections(filename)
    devouts[:out] = Symbol[]
    firstnode, lastnode = (:svr, :out)
    midnodes = (:dac, :fft)
    secondnode, thirdnode = if iszero(npaths_between(devouts, midnodes...))
        last(midnodes), first(midnodes)
    else
        midnodes
    end
    path = (firstnode, secondnode, thirdnode, lastnode)
    prod(ns->npaths_between(devouts, ns...), zip(path[1:(end-1)], path[2:end]))
end
@assert main_2("d11/sample2") |> isequal(2)
@assert main_2("d11/input") |> isequal(520476725037672)
@comment begin
    devouts = parse_connections("d11/input")
    npaths_to_out(:svr, devouts, Dict{Symbol, Int}((:out => 1))) # -> 355562200179074951
    # okay we're not iterating over those, lol
    #npaths_to_out(:svr, devouts, Dict{Symbol, Int}((:dac => 1))) # -> 1406220880404
    npaths_to_out(:dac, devouts, Dict{Symbol, Int}((:fft => 1))) # -> 0
    # ah, no loops (that would be infinite)
    # must go svr->fft->dac->out
    npaths_to_out(:svr, devouts, Dict{Symbol, Int}((:fft => 1))) # -> 13844
    npaths_to_out(:fft, devouts, Dict{Symbol, Int}((:dac => 1))) # -> 4504653
    npaths_to_out(:dac, devouts, Dict{Symbol, Int}((:out => 1))) # -> 8346
end

end

# Local Variables:
# julia-snail-repl-buffer: "*julia-aoc*"
# julia-snail-port: 10012
# End:
