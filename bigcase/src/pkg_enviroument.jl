using Pkg
Pkg.activate("/Users/yuanyiping/Documents/GitHub/unit_commitment_code/bigcase/.pkg/")
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
include("/Users/yuanyiping/Documents/GitHub/unit_commitment_code/bigcase/src/formatteddata.jl")
include("/Users/yuanyiping/Documents/GitHub/unit_commitment_code/bigcase/src/renewableenergysimulation.jl")
include("/Users/yuanyiping/Documents/GitHub/unit_commitment_code/bigcase/src/showboundrycase.jl")
include("/Users/yuanyiping/Documents/GitHub/unit_commitment_code/bigcase/src/readdatafromexcel.jl")
include("/Users/yuanyiping/Documents/GitHub/unit_commitment_code/bigcase/src/SUCuccommitmentmodel.jl")
include("/Users/yuanyiping/Documents/GitHub/unit_commitment_code/bigcase/src/FCUCuccommitmentmodel.jl")
include("/Users/yuanyiping/Documents/GitHub/unit_commitment_code/bigcase/src/casesploting.jl")
include("/Users/yuanyiping/Documents/GitHub/unit_commitment_code/bigcase/src/creatfrequencyconstraints.jl")
include("/Users/yuanyiping/Documents/GitHub/unit_commitment_code/bigcase/src/saveresult.jl")
include("/Users/yuanyiping/Documents/GitHub/unit_commitment_code/bigcase/src/BFLib_consideringFRlimit.jl")
include("/Users/yuanyiping/Documents/GitHub/unit_commitment_code/bigcase/src/enhance_FCUCuccommitmentmodel.jl")
include("/Users/yuanyiping/Documents/GitHub/unit_commitment_code/bigcase/src/generatefittingparameters.jl")

