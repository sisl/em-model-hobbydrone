using Pkg

if !in("BayesNets", keys(Pkg.installed()))
	Pkg.add("BayesNets")
end

if !in("LightGraphs", keys(Pkg.installed()))
	Pkg.add("LightGraphs")
end

if !in("DataFrames", keys(Pkg.installed()))
	Pkg.add("DataFrames")
end

if !in("CSV", keys(Pkg.installed()))
	Pkg.add("CSV")
end

if !in("JLD2", keys(Pkg.installed()))
	Pkg.add("JLD2")
end