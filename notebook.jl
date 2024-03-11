### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ ea189394-df64-11ee-2bb8-ef75ba1c2bf8
begin
	using Pkg
	Pkg.activate("HOTRG_nonorientablesurface/")
	using HOTRG_nonorientablesurface
	include("hotrg.jl")
	include("ising.jl")
	include("potts.jl")
	using LsqFit, Plots, Printf, LaTeXStrings
	using Logging
	Logging.disable_logging(Logging.Warn)
end;

# ╔═╡ 0e2327e1-4236-4c6f-aca8-811c10e5f1f9
begin
	theme(:dao)
	default(
	  fontfamily = "Times Roman",
	  guidefontsize = 12,
	  tickfontsize = 10,
	  legendfontsize = 10,
	  markerstrokewidth = 2,
	  margin = 5Plots.mm,
	  size = (600, 400),
	  grid = false,
	  foreground_color_legend = nothing,
	  background_color_legend = colorant"rgba(255,255,255,0.6)",
	  linewidth = 1
	)
end;

# ╔═╡ a0efd6cc-1153-4452-b64f-78f2c7439492
@__MODULE__

# ╔═╡ c1de7b4e-ead9-4a06-96eb-9a20bbecb54e
md"""
## Calculation of crosscap/rainbow free energy term for the Ising model and three-state Potts model
"""

# ╔═╡ 9f26526b-f9a9-4206-8944-f3119d5caab7
begin
	χ = 16
	stepnum = 7
	eigvalnum = 3;
end;

# ╔═╡ fd03016d-ae64-42b5-806e-032748972e6c
import Ising: bulk

# ╔═╡ a62dc594-440d-4ac5-8a37-31abf9d1485b
begin
	_, cftdata_ising = hotrg(Ising.bulk(inv(Ising.criticaltemperature()))...; maxdim = χ, stepnum, eigvalnum);
end

# ╔═╡ 47163114-deae-4769-8b9c-0b60886f85ff


# ╔═╡ Cell order:
# ╠═ea189394-df64-11ee-2bb8-ef75ba1c2bf8
# ╠═0e2327e1-4236-4c6f-aca8-811c10e5f1f9
# ╠═a0efd6cc-1153-4452-b64f-78f2c7439492
# ╟─c1de7b4e-ead9-4a06-96eb-9a20bbecb54e
# ╠═9f26526b-f9a9-4206-8944-f3119d5caab7
# ╠═fd03016d-ae64-42b5-806e-032748972e6c
# ╠═a62dc594-440d-4ac5-8a37-31abf9d1485b
# ╠═47163114-deae-4769-8b9c-0b60886f85ff
