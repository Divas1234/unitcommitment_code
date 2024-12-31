using Pkg
Pkg.activate("/.pkg/")
# Pkg.add("Revise")
# Pkg.add("JuMP")
# Pkg.add("Gurobi")
# Pkg.add("Test")
# Pkg.add("DelimitedFiles")
# Pkg.add("PlotlyJS")
# Pkg.add("LaTeXStrings")
# Pkg.add("Plots")
# Pkg.add("JLD")
# Pkg.add("DataFrames")
# Pkg.add("Clustering")
# Pkg.add("StatsPlots")
# Pkg.add("Distributions")
# Pkg.add("XLSX")
using Revise, JuMP, Gurobi, Test, DelimitedFiles, PlotlyJS, LaTeXStrings, Plots, JLD, DataFrames,Distributions, XLSX
using Clustering, StatsPlots

# using JuliaFormatter
# plotlyjs()
gr()
using Random
Random.seed!(1234)
include("./formatteddata.jl")
include("./renewableenergysimulation.jl")
include("./showboundrycase.jl")
include("./readdatafromexcel.jl")
include("./SUCuccommitmentmodel.jl")
include("./FCUCuccommitmentmodel.jl")
include("./casesploting.jl")
include("./creatfrequencyconstraints.jl")
include("./saveresult.jl")
include("./BFLib_consideringFRlimit.jl")
include("./enhance_FCUCuccommitmentmodel.jl")
include("./generatefittingparameters.jl")

