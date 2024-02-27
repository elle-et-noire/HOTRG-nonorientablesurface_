include("trgutils.jl")

function hotrg(A::ITensor, iA::Indices; maxdim, stepnum, eigvalnum, copybelow = copybelow_align, copyright = copyright_align)
  norms = zeros(stepnum)
  cftval = Dict(zip(["<C|i>", "<R|i>", "eigval"], [zeros(eigvalnum, stepnum) for _ in 1:3]))

  # init spatial reflection operator
  iO = Indices(up = altind(iA.up), down = altind(iA.down))
  O = δ(iO.up, iO.down)
  refl(i, j) = replaceinds(O, iO.up => i, iO.down => j)

  for i in 1:stepnum
    B, iB = copybelow(A, iA)
    AB = A * B

    # measure cft data of crosscap and rainbow boundary state
    M = AB * δ(iA.up, iB.down); iM1 = iA.left; iM2 = iB.left
    U1, S, _ = svd(M, (iM1, iM2); maxdim = eigvalnum)
    S = storage(S)
    norms[i] = S[1]; S /= norms[i]; A /= norms[i]
    D = length(S)
    cftval["eigval"][1:D, i] = S

    # if the copy is reflected, correspondence between crosscap/rainbow and contracting δ_ij/O_ij flips
    if copybelow == copybelow_align
      cftval["<C|i>"][1:D, i] = storage(U1 * δ(iM1, iM2))
      cftval["<R|i>"][1:D, i] = storage(U1 * refl(iM1, iM2))
    elseif copybelow == copybelow_reflect
      cftval["<R|i>"][1:D, i] = storage(U1 * δ(iM1, iM2))
      cftval["<C|i>"][1:D, i] = storage(U1 * refl(iM1, iM2))
    end

    U, iU = horizontalprojector(A, iA, B, iB; maxdim)
    A, iA = contract_horizontalprojector(AB, iA, iB, U, iU)

    B, iB = copyright(A, iA)
    V, iV = verticalprojector(A, iA, B, iB; maxdim)
    A, iA = contract_verticalprojector(A, iA, V, iV, B, iB)

    O, iO = contract_horizontalprojector_twister(O, iO, U, iU)
  end

  norms, cftval
end
