module AOC25_03

using Transducers

function maxjoltage_2bat(bank)
    b1 = argmax(@view(bank[begin:(end-1)]))
    vb1 = bank[b1] * 10
    vb2 = maximum(@view(bank[(b1+1):end]))
    vb1 + vb2
end

function maxjoltage_12bat(bank)
    inds = zeros(Int, 12)
    last_ind = 0
    for ii in 1:12
        _avail = @view(bank[(last_ind+1):(end+1-(12-(ii-1)))])
        inds[ii] = last_ind + argmax(_avail)
        last_ind = inds[ii]
    end
    evalpoly(10, reverse(view(bank, inds)))
end

function main_1(filename)
    eachline(filename) |> Map() do l
        l |> Map(c->parse(Int, c)) |> collect
    end |> Map(maxjoltage_2bat) |>
    foldxl(+)
    #collect
end

function main_2(filename)
    eachline(filename) |> Map() do l
        l |> Map(c->parse(Int, c)) |> collect
    end |> Map(maxjoltage_12bat) |>
    foldxl(+)
    #collect
end

end

# println(AOC25_03.main_1("sample") == 357)
# println(AOC25_03.main_1("input"))

#println(AOC25_03.main_2("sample") == 3121910778619)
