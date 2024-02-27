using LsqFit, Plots, Printf

function read(path)
  dic = Dict()
  vecs = []
  label = ""
  open(path, "r") do fp
    for line in eachline(fp)
      if isempty(line)
        continue
      end
      if line[1] == '#'
        if !isempty(label)
          dic[label] = vcat(vecs...)
          vecs = []
        end
        label = line[3:findlast(':', line) - 1]
        if isnumeric(line[end])
          dic[label] = parse(Float64, split(line)[end])
          label = ""
        end
      else
        push!(vecs, parse.(Float64, split(line))')
      end
    end
    dic[label] = vcat(vecs...)
  end
  dic
end

function getfolder(path)
  i = findlast('/', path)
  isnothing(i) ? "" : path[1:(i - 1)]
end

function memo(path, dic)
  mkpath(getfolder(path))
  open(path, "w") do fp
    for (key, value) in dic
      if isempty(size(value))
        println(fp, "# $key: $value")
      else
        println(fp, "# $key:")
        Base.print_array(fp, value)
        println(fp)
      end
    end
  end
end

function savefig!(path)
  mkpath(getfolder(path))
  savefig(path)
end

function savefig!(plt, path)
  mkpath(getfolder(path))
  savefig(plt, path)
end

function f2s(x)
  @sprintf("%.3f", x)
end

function linfit(x, y)
  f(x, p) = @. p[1] * x + p[2]
  fit = curve_fit(f, x, y, [0., 0])
  x -> f(x, fit.param), fit.param...
end