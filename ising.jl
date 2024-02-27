include("trgutils.jl")
using QuadGK

module Ising

function weight(β)
  [exp(β) exp(-β); exp(-β) exp(β)]
end

function bulk(β)
  Main.bulk(weight(β))
end

function horizontalboundary(β)
  Main.horizontalboundary(weight(β))
end

function criticaltemperature()
  2 / log(1 + √2)
end

function freeenergy(β)
  (Main.quadgk(k1 ->
  Main.quadgk(k2 ->
    log(cosh(2β)^2 - sinh(2β) * (cos(k1) + cos(k2))) / 8pi^2,
  0, 2pi)[1],
  0, 2pi)[1] + log(2)) / -β
end

function centralcharge()
  0.5
end

function quantumdimension()
  1 + 1 / √2
end

end