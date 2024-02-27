include("trgutils.jl")

module Potts

function weight(β; q)
  Main.Diagonal([exp(β) - 1 for _ in 1:q]) + ones(q, q)
end

function bulk(β; q)
  Main.bulk(weight(β; q))
end

function horizontalboundary(β; q)
  Main.horizontalboundary(weight(β; q))
end

function criticaltemperature(;q)
  1 / log(1 + √q)
end

function centralcharge(;q)
  if q == 3
    return 0.8
  elseif q == 4
    return 1.0
  end
end

function quantumdimension(;q)
  if q == 3
    return sqrt(3 + 6 / √5)
  elseif q == 4
    return (3 + 2√2) / 2
  end
end

end