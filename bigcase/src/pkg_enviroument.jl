using Pkg

# Activate the project environment
Pkg.activate(joinpath(@__DIR__, "../.pkg/"))

# Import necessary packages
using Revise, JuMP, Gurobi, Test, DelimitedFiles, CSV, LaTeXStrings, Plots, JLD,
	  DataFrames, Distributions, XLSX, Clustering, StatsPlots

# Setup plotting
gr()

# Set random seed for reproducibility
using Random
Random.seed!(1234)

if Sys.isapple()
	# Include necessary source files
	base_path = joinpath(@__DIR__, "../src")
	include(joinpath(base_path, "formatteddata.jl"))
	include(joinpath(base_path, "renewableenergysimulation.jl"))
	include(joinpath(base_path, "showboundrycase.jl"))
	include(joinpath(base_path, "readdatafromexcel.jl"))
	include(joinpath(base_path, "SUCuccommitmentmodel.jl"))
	include(joinpath(base_path, "FCUCuccommitmentmodel.jl"))
	include(joinpath(base_path, "casesploting.jl"))
	include(joinpath(base_path, "creatfrequencyconstraints.jl"))
	include(joinpath(base_path, "saveresult.jl"))
	include(joinpath(base_path, "BFLib_consideringFRlimit.jl"))
	include(joinpath(base_path, "enhance_FCUCuccommitmentmodel.jl"))
	include(joinpath(base_path, "generatefittingparameters.jl"))
elseif Sys.iswindows()
	include("C:/Users/yyp_uestc/Downloads/unitcommitment_code-master/bigcase/src/formatteddata.jl")
	include("C:/Users/yyp_uestc/Downloads/unitcommitment_code-master/bigcase/src/renewableenergysimulation.jl")
	include("C:/Users/yyp_uestc/Downloads/unitcommitment_code-master/bigcase/src/showboundrycase.jl")
	include("C:/Users/yyp_uestc/Downloads/unitcommitment_code-master/bigcase/src/readdatafromexcel.jl")
	include("C:/Users/yyp_uestc/Downloads/unitcommitment_code-master/bigcase/src/SUCuccommitmentmodel.jl")
	include("C:/Users/yyp_uestc/Downloads/unitcommitment_code-master/bigcase/src/FCUCuccommitmentmodel.jl")
	include("C:/Users/yyp_uestc/Downloads/unitcommitment_code-master/bigcase/src/casesploting.jl")
	include("C:/Users/yyp_uestc/Downloads/unitcommitment_code-master/bigcase/src/creatfrequencyconstraints.jl")
	include("C:/Users/yyp_uestc/Downloads/unitcommitment_code-master/bigcase/src/saveresult.jl")
	include("C:/Users/yyp_uestc/Downloads/unitcommitment_code-master/bigcase/src/BFLib_consideringFRlimit.jl")
	include("C:/Users/yyp_uestc/Downloads/unitcommitment_code-master/bigcase/src/enhance_FCUCuccommitmentmodel.jl")
	include("C:/Users/yyp_uestc/Downloads/unitcommitment_code-master/bigcase/src/generatefittingparameters.jl")
end
