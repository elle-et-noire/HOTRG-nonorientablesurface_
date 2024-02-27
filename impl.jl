include("hotrg.jl")
include("ising.jl")
include("potts.jl")
include("utils.jl")

default(
  fontfamily="Computer Modern",
  guidefontsize=12,
  tickfontsize=10,
  legendfontsize=8,
  margin=5Plots.mm,
  size=(600, 400),
  grid = false,
  foreground_color_legend = nothing,
  background_color_legend = colorant"rgba(255,255,255,0.6)"
)

function impl(;maxdim, stepnum, model, eigvalnum = 3)
  if model == :ising
    β = inv(Ising.criticaltemperature())
    A, iA = Ising.bulk(β)
  end
  if String(model)[1:5] == "potts"
    q = parse(Int, String(model)[end])
    β = inv(Potts.criticaltemperature(; q))
    A, iA = Potts.bulk(β; q)
  end

  norms, cftdata = hotrg(A, iA; maxdim, stepnum, eigvalnum)

  items = Any[("norm", norms)]
  for (way, d) in cftdata
    for (key, val) in d
      push!(items, ("$key ($way)", val))
    end
  end

  memo("memo/$model/$maxdim.txt", sort(items))
end

markershapes = Dict()
function markershape(label)
  if !haskey(markershapes, label)
    for s in [:circle, :star5, :utriangle, :diamond, :star4, :hexagon]
      if !(s in values(markershapes))
        global markershapes[label] = s
      end
    end
  end
  markershapes[label]
end

markercolors = Dict()
function markercolor(label)
  if !haskey(markercolors, label)
    for c in [:crimson, :blue, :green, :orange, :violet, :gray]
      if !(c in values(markercolors))
        global markercolors[label] = c
      end
    end
  end
  markercolors[label]
end


function plotfr(;model, maxdim, fitrange, plotrange = nothing)
  data = read("memo/$model/$maxdim.txt")

  c = 0.
  if model == :ising
    c = Ising.centralcharge()
  elseif String(model)[1:5] == "potts"
    q = parse(Int, String(model)[end])
    c = Potts.centralcharge(;q)
  end

  frplt = plot()
  frerrplt = plot()
  for (key, val) in data
    if length(key) < 5 || key[1:5] != "<R|i>"
      continue
    end
    label = string(split(key)[2:end]...)

    fr = log.(abs.(val[1, :]))
    N = collect(eachindex(fr))
    β = 2 .^ N

    f_test(x, p) = @. (c / 4) * log(x) + p[1]
    fit = curve_fit(f_test, β[fitrange], fr[fitrange], [0.])
    b = fit.param[1]
    println("b = $b")

    if isnothing(plotrange)
      plotrange = eachindex(fr)
    end

    scatter!(frplt, β[plotrange], fr[plotrange], xlabel = "β", ylabel = "F_R",
      title = "$model, D = $maxdim, fit in $fitrange, $method"; label, xscale = :log10,
      legend = :topleft, markershape = markershape(label), color = markercolor(label)
    )
    plot!(frplt, β[plotrange], f_test(β[plotrange], fit.param), label = "y = (c / 4) ln β + $(f2s(b)), c = $c")

    err = abs.(fr - f_test(β, fit.param))
    plot!(frerrplt, β[plotrange], err[plotrange], xlabel = "β", ylabel = "error of F_R",
      title = "$model, D = $maxdim, fit in $fitrange, $method"; label, xscale = :log10,
      markershape = markershape(label), color = markercolor(label)
    )
  end

  savefig!(frplt, "plot/$method/$model/fr-maxdim$maxdim.png")
  savefig!(frerrplt, "plot/$method/$model/frerr-maxdim$maxdim.png")
end

function plotfc(;model, maxdim, method = nonorientablehotrg, plotrange, plotlevels = [1])
  data = read("memo/$method/$model/$maxdim.txt")

  g = 0.
  if model == :ising
    g = Ising.quantumdimension()
  elseif String(model)[1:5] == "potts"
    q = parse(Int, String(model)[end])
    g = Potts.quantumdimension(;q)
  end

  plot()
  for (key, val) in data
    if length(key) <= 5 || key[1:5] != "<C|i>"
      continue
    end
    label = string(split(key)[2:end]...)

    gamma2 = val .^ 2
    N = collect(eachindex(gamma2))
    β = 2 .^ N

    plot!(β[plotrange], [gamma2[j, plotrange] for j in plotlevels], xlabel = "β", ylabel = "Γ^2",
      title = "$model, D = $maxdim, $method"; label, xscale = :log10,
      markershape = markershape(label), color = markercolor(label)
    )
  end
  vals_theory = [g]
  if model == :ising
    append!(vals_theory, [0, 1 - 1 / √2])
  end
  label = vals_theory .|> x -> @sprintf("%.5f", x)
  for i in eachindex(vals_theory)
    if i in plotlevels
      hline!([vals_theory[i]], label = label[i])
    end
  end

  savefig!("plot/$method/$model/gamma2-maxdim$maxdim.png")
end

function plotall(;model, maxdim, method = nonorientablehotrg, fcplotrange, fcplotlevels = [1], fitrange)
  plotfc(;model, maxdim, method, plotrange = fcplotrange, plotlevels = fcplotlevels)
  plotfr(;model, maxdim, method, fitrange)
end

if length(ARGS) < 4
  println("usage: julia impl.jl %model% %maxdim% %stepnum% %method%")
else
  impl(;maxdim = parse(Int64, ARGS[2]), stepnum = parse(Int64, ARGS[3]), method = getfield(Main, Symbol(ARGS[4])), model = Symbol(ARGS[1]))
end