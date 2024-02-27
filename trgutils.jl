using ITensors, LinearAlgebra, Parameters

@with_kw mutable struct Indices
  up = missing
  left = missing
  down = missing
  right = missing
  proj = missing
end

function i2v(iA::Indices)
  collect(skipmissing(propertynames(iA) .|> x -> getproperty(iA, x)))
end

function altind(i)
  Index(dim(i))
end

function altinds(iA::Indices, default...)
  iB = Indices()
  for (property, val) in default
    setproperty!(iB, property, val)
  end
  for property in propertynames(iA)
    if ismissing(getproperty(iA, property)) || !ismissing(getproperty(iB, property))
      continue
    end
    setproperty!(iB, property, altind(getproperty(iA, property)))
  end
  iB
end

function replaceindices(A::ITensor, iA::Indices, iB::Indices)
  replaceinds(A, i2v(iA), i2v(iB))
end

function replaceindices(A::ITensor, iA::Indices, rep_inds::Pair{Symbol, <:Index}...)
  for (prop, ind) in rep_inds
    replaceinds(A, getproperty(iA, prop) => ind)
    setproperty!(iA, prop, ind)
  end
  A, iA
end

function replaceindices!(A::ITensor, iA::Indices, iB::Indices)
  replaceinds!(A, i2v(iA), i2v(iB))
end

function copybelow_align(A::ITensor, iA::Indices)
  iB = altinds(iA, :up => iA.down)
  replaceindices(A, iA, iB), iB
end

function copybelow_reflect(A::ITensor, iA::Indices)
  iB = altinds(iA, :up => iA.down)
  replaceinds(A, iA.up => iB.down, iA.left => iB.left, iA.right => iB.right), iB
end

function copyright_align(A::ITensor, iA::Indices)
  iB = altinds(iA, :left => iA.right)
  replaceindices(A, iA, iB), iB
end

function copyright_reflect(A::ITensor, iA::Indices)
  iB = altinds(iA, :left => iA.right)
  replaceinds(A, iA.up => iB.up, iA.down => iB.down, iA.left => iB.right), iB
end

"""
A is symmetryc with iA.left <-> iA.right (same for B).
"""
function horizontalprojector(A::ITensor, iA::Indices, B = missing, iB = missing; maxdim)
  if ismissing(B)
    B, iB = copybelow_align(A, iA)
  end
  @assert iA.down == iB.up

  AA = A * prime(conj(A), iA.left, iA.down)
  BB = B * prime(conj(B), iB.left, iB.up)
  U, S, _ = svd(AA * BB, (iA.left, iB.left); maxdim)
  iU = Indices(up = iA.left, down = iB.left, proj = commonind(U, S))
  U, iU
end

"""
A is symmetryc with iA.up <-> iA.down (same for B).
"""
function verticalprojector(A::ITensor, iA::Indices, B = missing, iB = missing; maxdim)
  if ismissing(B)
    B, iB = copyright_align(A, iA)
  end
  @assert iA.right == iB.left

  AA = A * prime(conj(A), iA.up, iA.right)
  BB = B * prime(conj(B), iB.up, iB.left)
  U, S, _ = svd(AA * BB, (iA.up, iB.up); maxdim)
  iU = Indices(left = iA.up, right = iB.up, proj = commonind(U, S))
  U, iU
end

function contract_horizontalprojector(A::ITensor, iA::Indices, U::ITensor, iU::Indices, B = missing, iB = missing)
  if ismissing(B)
    B, iB = copybelow_align(A, iA)
  end
  @assert iA.down == iB.up
  @assert !ismissing(iU.up) && !ismissing(iU.down) && !ismissing(iU.proj)

  A *= B
  A *= replaceinds(U, iU.up => iA.left, iU.down => iB.left)
  iUproj = altind(iU.proj)
  A *= replaceinds(U, iU.up => iA.right, iU.down => iB.right, iU.proj => iUproj)
  A, Indices(up = iA.up, left = iU.proj, right = iUproj, down = iB.down)
end

function contract_horizontalprojector(AB::ITensor, iA::Indices, iB::Indices, U::ITensor, iU::Indices)
  @assert iA.down == iB.up
  @assert hassameinds((iA.up, iA.left, iA.right, iB.left, iB.right, iB.down), AB)
  @assert !ismissing(iU.up) && !ismissing(iU.down) && !ismissing(iU.proj)

  AB *= replaceinds(U, iU.up => iA.left, iU.down => iB.left)
  iUproj = altind(iU.proj)
  AB *= replaceinds(U, iU.up => iA.right, iU.down => iB.right, iU.proj => iUproj)
  AB, Indices(up = iA.up, left = iU.proj, right = iUproj, down = iB.down)
end

function contract_verticalprojector(A::ITensor, iA::Indices, U::ITensor, iU::Indices, B = missing, iB = missing)
  if ismissing(B)
    B, iB = copyright_align(A, iA)
  end
  @assert iA.right == iB.left
  @assert !ismissing(iU.left) && !ismissing(iU.right) && !ismissing(iU.proj)

  A *= B
  A *= replaceinds(U, iU.left => iA.up, iU.right => iB.up)
  iUproj = altind(iU.proj)
  A *= replaceinds(U, iU.left => iA.down, iU.right => iB.down, iU.proj => iUproj)
  A, Indices(up = iU.proj, left = iA.left, right = iB.right, down = iUproj)
end

function contract_verticalprojector_twister(O::ITensor, iO::Indices, U::ITensor, iU::Indices)
  @assert !ismissing(iU.left) && !ismissing(iU.right) && !ismissing(iU.proj)
  @assert !ismissing(iO.up) && !ismissing(iO.down)

  _iO = altinds(iO)
  _O = replaceindices(O, iO, _iO)
  O *= _O * replaceinds(U, iU.left => iO.up, iU.right => _iO.up)
  iUproj = altind(iU.proj)
  O *= replaceinds(U, iU.left => _iO.down, iU.right => iO.down, iU.proj => iUproj)
  iO = Indices(up = iU.proj, down = iUproj)
  O, iO
end

function contract_horizontalprojector_twister(O::ITensor, iO::Indices, U::ITensor, iU::Indices)
  @assert !ismissing(iU.up) && !ismissing(iU.down) && !ismissing(iU.proj)
  @assert !ismissing(iO.right) && !ismissing(iO.left)

  _iO = altinds(iO)
  _O = replaceindices(O, iO, _iO)
  O *= _O * replaceinds(U, iU.up => iO.left, iU.down => _iO.left)
  iUproj = altind(iU.proj)
  O *= replaceinds(U, iU.up => _iO.right, iU.down => iO.right, iU.proj => iUproj)
  iO = Indices(left = iU.proj, right = iUproj)
  O, iO
end

function trg(A...; maxdim, stepnum, step, normalize, observer = [])
  norms = []
  observable = [[] for _ in eachindex(observer)]

  for i in 1:stepnum
    A = step(A...; maxdim)
    n, A = normalize(A...)
    push!(norms, n)
    for i in eachindex(observer)
      push!(observable[i], observer[i](A...))
    end
  end
  norms, observable..., A...
end

function logpartfunc(norms; sitenum_per_step, initial_sitenum = 1)
  stepnum = length(norms)
  sum = 0.
  sitenum = initial_sitenum
  lnz = Vector{Float64}(undef, stepnum)

  for i in 1:stepnum
    sitenum *= sitenum_per_step
    sum += log(norms[i]) / sitenum
    lnz[i] = sum
  end
  lnz
end

function scalingval(lnz, vals = ones(length(lnz)); sitenum_per_step, initial_sitenum = 1)
  -[lnz[i] * sitenum_per_step^i * initial_sitenum + log(vals[i]) for i in eachindex(lnz)]
end

function bulk(vweight, hweight = vweight)
  dv = size(vweight)[1]
  dh = size(hweight)[1]
  A = δ(Index.([dv, dh, dv, dh]))

  V = ITensor(vweight, inds(A)[1], inds(A)[3])
  U1, S, V1 = svd(V, inds(V)[1])
  iU1 = commonind(S, V1)
  iV1 = commonind(S, U1)
  U1 *= sqrt.(S)
  V1 *= sqrt.(S)

  H = ITensor(hweight, inds(A)[2], inds(A)[4])
  U2, S, V2 = svd(H, inds(H)[1])
  iU2 = commonind(S, V2)
  iV2 = commonind(S, U2)
  U2 *= sqrt.(S)
  V2 *= sqrt.(S)

  A *= U1
  A *= V1
  A *= U2
  A *= V2

  iA = Indices(up = iU1, down = iV1, left = iU2, right = iV2)
  A, iA
end

function horizontalboundary(vweight, hweight = vweight)
  dv = size(vweight)[1]
  dh = size(hweight)[1]
  A = δ(Index.([dh, dv, dh]))

  V = ITensor(vweight, inds(A)[2], inds(A)[2]')
  U1, S, V1 = svd(V, inds(V)[1])
  iU1 = commonind(S, V1)
  U1 *= sqrt.(S)

  H = ITensor(hweight, inds(A)[1], inds(A)[3])
  U2, S, V2 = svd(H, inds(H)[1])
  iU2 = commonind(S, V2)
  iV2 = commonind(S, U2)
  U2 *= sqrt.(S)
  V2 *= sqrt.(S)

  A *= U1
  A *= U2
  A *= V2

  iA = Indices(down = iU1, left = iU2, right = iV2)
  A, iA
end

function phasetransdetector(M)
  tr(M)^2 / tr(M^2)
end

function trace_torus(A::ITensor, iA::Indices)
  tr(transfermatrix_torus(A, iA))
end

function normalize_torus(A::ITensor, iA::Indices)
  f = trace_torus(A, iA)
  A /= f
  f, (A, iA)
end

function transfermatrix_torus(A::ITensor, iA::Indices; direction = :horizontal)
  if direction == :horizontal
    matrix(A * δ(iA.up, iA.down))
  elseif direction == :vertical
    matrix(A * δ(iA.left, iA.right))
  end
end