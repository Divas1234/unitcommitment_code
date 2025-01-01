include("src/pkg_enviroument.jl")

# using Debugger
UnitsFreqParam, WindsFreqParam, StrogeData, DataGen, GenCost, DataBranch, LoadCurve, DataLoad = readxlssheet()
config_param, units, lines, loads, stroges, NB, NG, NL, ND, NT, NC = forminputdata(
	DataGen, DataBranch, DataLoad, LoadCurve, GenCost, UnitsFreqParam, StrogeData)

winds, NW = genscenario(WindsFreqParam, 1)

#* calculating different value along with different residual settings.
residual_scenarios_num = 5

# e_x₀, e_p₀, e_pᵨ, e_pᵩ, e_seq_sr⁺, e_seq_sr⁻, e_pss_charge_p⁺, e_pss_charge_p⁻, e_su_cost, e_sd_cost, e_prod_cost, e_cost_sr⁺, e_cost_sr⁻ =
# 	enhance_FCUC_scucmodel(NT, NB, NG, ND, NC, units, loads, winds, lines, config_param)

# #! NOTE basesline UC and FCUC
# x₀, p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, pss_charge_p⁺, pss_charge_p⁻, su_cost, sd_cost, prod_cost, cost_sr⁺, cost_sr⁻ = FCUC_scucmodel(NT, NB, NG, ND, NC, units, loads, winds, lines, config_param)

big_M = 1e6
f_base = 50.0
RoCoF_max = 1.5
f_nadir = 49.5
f_qss = 49.5
Δp = maximum(units.p_max[:, 1]) * 1.0
# ------
x = ones(NG, 1)
x[1:30] .= 0

sum(winds.Mw[:, 1] .* winds.Fcmode[:, 1] .* winds.p_max[:, 1]) +
2 * sum(x .* units.Hg[:, 1] .* units.p_max[:, 1])
Δp * f_base / RoCoF_max * (sum(units.p_max[:, 1]) + sum(winds.Fcmode .* winds.p_max)) / 50

flag_method_type = 1
NN = 1 # sampled scenarios for chance constraints
μ, σ = 0, 0.5e-3
fittingparameter_vector, whitenoise_parameter, whitenoise_parameter_probability = generatefreq_fittingparameters(
	units, winds, NG, NW, NN, flag_method_type, μ, σ)

adjustablewindsVSCpower = winds.Fcmode .* winds.p_max
inverse_winds_Rw = zeros(NW, 1)
for i in 1:NW
	if winds.Fcmode[i, 1] == 0
		inverse_winds_Rw[i, 1] = 1 / winds.Rw[i, 1]
	end
end
current_Kw = 1.0
current_Dw = sum(winds.Dw .* adjustablewindsVSCpower) / sum(winds.p_max) # Dw
current_Mw = sum(winds.Mw .* adjustablewindsVSCpower) / sum(winds.p_max) # Mw
current_Hw = current_Mw / 2
current_Rw = 1 / (sum(winds.Kw .* inverse_winds_Rw .* (ones(NW, 1) - winds.Fcmode) .*
				  winds.p_max) /
			  sum((ones(NW, 1) .* winds.p_max)))

#  powers for intia frequency response
localapparentpower = (sum(units.p_max[:, 1]) + sum(winds.p_max .* winds.Fcmode))
sumapparentpower = (localapparentpower - sum(winds.p_max .* winds.Fcmode) +
					sum(winds.p_max))
fittingparameter = fittingparameter_vector[1, :] * (-1)

fittingparameter[1] / sumapparentpower *
(sum(x .* units.Hg .* units.p_max) + sum(current_Mw .* adjustablewindsVSCpower)) +
fittingparameter[2] / sum(units.p_max) *
(sum(x .* units.Kg .* units.Fg ./ units.Rg .* units.p_max)) +
fittingparameter[3] / sum(units.p_max) *
(sum(x .* units.Kg ./ units.Rg .* units.p_max)) +
fittingparameter[4]
(f_base - f_nadir) * 2.0

# -----
# FIXME redebug
# Sampling_Statue[1, 1] = 1
# creatingfrequencyresponsesamplingdata(units, winds, NW, NG, Sampling_Statue)

# normalized winds parameters through COI
vsmFC_number = sum(winds.Fcmode[:, 1])
doopFC_number = length(winds.Fcmode[:, 1]) - vsmFC_number
adjustablewindsVSCpower = winds.Fcmode .* winds.p_max
inverse_winds_Rw = zeros(NW, 1)
for i in 1:NW
	if winds.Fcmode[i, 1] == 0
		inverse_winds_Rw[i, 1] = 1 / winds.Rw[i, 1]
	end
end
current_Kw = 1.0
current_Dw = sum(winds.Dw .* adjustablewindsVSCpower) / sum(adjustablewindsVSCpower) # Dw
current_Mw = sum(winds.Mw .* adjustablewindsVSCpower) / sum(adjustablewindsVSCpower) # Mw
current_Hw = current_Mw / 2
current_Rw = 1 / (sum(winds.Kw .* inverse_winds_Rw .* (ones(NW, 1) - winds.Fcmode) .*
				  winds.p_max) / sum((ones(NW, 1) - winds.Fcmode) .* winds.p_max))

# units parameters
Sampling_Statue = ones(NG, 1)
adjustabletheramlpower = units.p_max .* Sampling_Statue

indices_to_PSSunits = zeros(NG, 1)
indices_to_PSSunits[findall(x -> x < 100, units.Rg)] .= 1

current_Kg = 1.0
current_Tg = mean(units.Tg)
current_Fg_div_Rg = sum(indices_to_PSSunits .* units.Kg .* units.Fg ./ units.Rg .*
						adjustabletheramlpower) /
					sum(adjustabletheramlpower)
inverse_current_Rg = 1 / (sum(indices_to_PSSunits .* units.Kg ./ units.Rg .*
						  adjustabletheramlpower) / sum(adjustabletheramlpower)) # Kg
current_Rg = inverse_current_Rg
current_Fg = current_Fg_div_Rg * current_Rg
current_Dg = sum(units.Dg .* adjustabletheramlpower) / sum(adjustabletheramlpower)
current_Hg = sum(units.Hg .* adjustabletheramlpower) / sum(adjustabletheramlpower)
current_Mg = current_Hg * 2

#  powers for intia frequency response
localapparentpower = (sum(units.p_max[:, 1]) + sum(winds.p_max .* winds.Fcmode))
sumapparentpower   = (localapparentpower - sum(winds.p_max .* winds.Fcmode) + sum(winds.p_max))
p_step             = maximum(units.p_max) * 1.0

# sumD and sumH
adjustabletheramlpower
current_sumD = (sum(units.Dg .* adjustabletheramlpower) +
				sum(winds.Dw .* winds.p_max .* winds.Fcmode +
					winds.Kw .* winds.p_max .* (ones(NW, 1) - winds.Fcmode)) * 20) /
			   (sum(winds.p_max) + sum(adjustabletheramlpower))
current_sumH = (sum(current_Mg .* adjustabletheramlpower) +
				sum(current_Mw .* adjustablewindsVSCpower)) /
			   (sum(adjustablewindsVSCpower) + sum(adjustabletheramlpower)) / 2

D, H, F, R, T, K, δp = current_sumD,
current_sumH, current_Fg, current_Rg, current_Tg, current_Kg, p_step
# H = H * 4
# T = 7
# R = 1/R

R = collect(0:0.05:10)
F = collect(0:0.05:10)
ll = length(R)
res = zeros(ll, ll)
res1 = zeros(ll, ll)
H = 4.96
D = 1.2
for i in 1:ll
	for j in 1:ll
		ωₙ = sqrt((D * R[i] + 1) / (2 * H * R[i] * T))
		ς = (D * R[i] * T + 2 * H * R[i] + F[j] * T) / (2 * (D * R[i] + 1)) * ωₙ
		# ζ = (2 * H + (F * (K / R) + D) * T) / (2 * sqrt(2 * H * T) * sqrt(D + K / R))
		if ς <= 1
			res[i,j] = ς
			ωᵣ = ωₙ * sqrt(1 - ς * ς)
			tem = 2 * T * ς * ωₙ + (T * ωᵣ) * (T * ωᵣ)
			# res1[i, j] = tem
			if tem < 1
				# α = sqrt((1 - tem) / (1 - ς * ς))
				res1[i, j] = tem
			end
		end
	end
end
Plots.plot(res, legend=false)
Plots.plot(res1, legend=false)
@show res
@show res1











# FIXME - sfr -----
# @show [D, H, F, T, R, δp]

D, M, H, Tᵣ, Fₕ, R, K = 1.2, 4.96 * 2, 4.96, 10.8, 0.24, 1 / 20.8, 1.00
# D, M, H, T, F, R, K = Dg, Mg, Hg, Tg, Fg, Rg, K
D = 4.9567
H = 8
F = 0.2898
R = 0.02065

f_nadir, t_nadir = calculation_frequencynadirdata(H, R, D, F, T, δp)

f_nadir = f_nadir * (-1)

ωᵣ = ωₙ * sqrt(1 - ς * ς)
α = sqrt((1 - 2 * T * ς * ωₙ + (T * ωᵣ) * (T * ωᵣ)) / (1 - ς * ς))

(D * R * T + 2 * H * R + F * T)
(2 * (D * R + 1))

D / H
2 / T
F / R / H

ωᵣ = ωₙ * sqrt(1 - ς * ς)

# δp = 1.10
# H, R, D, F, T = equivalentfrequencycoefficients[i, 1:5]
ωₙ = sqrt((D * R + 1) / (2 * H * R * T))
ς = (D * R * T + 2 * H * R + F * T) / (2 * (D * R + 1)) * ωₙ
ωᵣ = ωₙ * sqrt(1 - ς * ς)
α = sqrt((1 - 2 * T * ς * ωₙ + (T * ωᵣ) * (T * ωᵣ)) / (1 - ς * ς))
tnadir = abs((1 / ωᵣ) * (atan((ωᵣ * T) / (ς * ωᵣ * T - 1))))
fnadir = δp * (R) / (D * R + 1) * (1 + sqrt(1 - ς * ς) * α * exp(-ς * ωₙ * tnadir))

T = 20
1 - 2 * T * ς * ωₙ + (T * ωᵣ) * (T * ωᵣ)
(T * ωᵣ) * (T * ωᵣ)
2 * T * ς * ωₙ

α = sqrt((1 - 2 * T * ς * ωₙ + (T * ωᵣ) * (T * ωᵣ)) / (1 - ς * ς))
tnadir = abs((1 / ωᵣ) * (atan((ωᵣ * T) / (ς * ωᵣ * T - 1))))
fnadir = δp * (R) / (D * R + 1) * (1 + sqrt(1 - ς * ς) * α * exp(-ς * ωₙ * tnadir))

D, M, H, Tᵣ, Fₕ, R, K = Dg, Mg, Hg, Tg, Fg, Rg, Kg
R = Rg / Kg
δp = δp * 1.0
wₙ = sqrt((D + K / R) / (2 * H * Tᵣ))
ζ = (2 * H + (Fₕ * (K / R) + D) * Tᵣ) / (2 * sqrt(2 * H * Tᵣ) * sqrt(D + K / R))
# ζ = (D * R * Tᵣ + 2 * H * R + Fₕ * Tᵣ) / 2 / (D * R + 1) * wₙ
wᵣ = wₙ * sqrt(1 - ζ^2)
ψ = asin(sqrt(1 - ζ^2))

D / H
2 / Tᵣ
Fₕ / R / H
D * R

xdata = collect(0:δt:endtime)
ydata = zeros(size(xdata, 1), 1)
f_base = 50

for i in 1:size(xdata, 1)
	t = xdata[i, 1]
	# δf = R * δp / (D * R + 1)
	# δf = δf * (1 + α * exp(-1.0 * ζ * wₙ * t) * sin(wᵣ * 1.0 * t + ψ))
	δf = δp / (2 * H * Tᵣ * (wₙ^2))
	δf = δf +
		 δp / (2 * Hg * wᵣ) * exp(-ζ * wₙ * t) *
		 (sin(wᵣ * t) - 1 / (wₙ * Tᵣ) * sin(wᵣ * t + ψ))
	ydata[i, 1] = f_base - δf
end

Plots.plot(ydata)