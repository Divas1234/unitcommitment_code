include("src/pkg_enviroument.jl")
using Debugger
UnitsFreqParam, WindsFreqParam, StrogeData, DataGen, GenCost, DataBranch, LoadCurve, DataLoad = readxlssheet()
config_param, units, lines, loads, stroges, NB, NG, NL, ND, NT, NC = forminputdata(DataGen, DataBranch, DataLoad, LoadCurve, GenCost, UnitsFreqParam, StrogeData)

winds, NW = genscenario(WindsFreqParam, 1)

#* calculating different value along with different residual settings.
residual_scenarios_num = 5

e_x₀, e_p₀, e_pᵨ, e_pᵩ, e_seq_sr⁺, e_seq_sr⁻, e_pss_charge_p⁺, e_pss_charge_p⁻, e_su_cost, e_sd_cost, e_prod_cost, e_cost_sr⁺, e_cost_sr⁻ =
	enhance_FCUC_scucmodel(NT, NB, NG, ND, NC, units, loads, winds, lines, config_param)

#! NOTE basesline UC and FCUC
x₀, p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, pss_charge_p⁺, pss_charge_p⁻, su_cost, sd_cost, prod_cost, cost_sr⁺, cost_sr⁻ = FCUC_scucmodel(NT, NB, NG, ND, NC, units, loads, winds, lines, config_param)

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
fittingparameter_vector, whitenoise_parameter, whitenoise_parameter_probability =
	generatefreq_fittingparameters(units, winds, NG, NW, NN, flag_method_type, μ, σ)

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
current_Dw = sum(winds.Dw .* adjustablewindsVSCpower) / sum(winds.p_max) # Dw
current_Mw = sum(winds.Mw .* adjustablewindsVSCpower) / sum(winds.p_max) # Mw
current_Hw = current_Mw / 2
current_Rw =
	1 / (sum(winds.Kw .* inverse_winds_Rw .* (ones(NW, 1) - winds.Fcmode) .* winds.p_max) /
		 sum((ones(NW, 1) .* winds.p_max)))

#  powers for intia frequency response
localapparentpower = (sum(units.p_max[:, 1]) + sum(winds.p_max .* winds.Fcmode))
sumapparentpower = (localapparentpower - sum(winds.p_max .* winds.Fcmode) + sum(winds.p_max))
fittingparameter = fittingparameter_vector[1, :] * (-1)

fittingparameter[1] / sumapparentpower *
(sum(x .* units.Hg .* units.p_max) + sum(current_Mw .* adjustablewindsVSCpower)) +
fittingparameter[2] / sum(units.p_max) *
(sum(x .* units.Kg .* units.Fg ./ units.Rg .* units.p_max)) +
fittingparameter[3] / sum(units.p_max) *
(sum(x .* units.Kg ./ units.Rg .* units.p_max)) +
fittingparameter[4]
(f_base - f_nadir) * 2.0

#FIXME - sfr -----
whitenoise_parameter = rand(Normal(μ, σ), NN)
fittingparameter_vector[1, :] = creatfrequencyfittingfunction(units, winds, NG, NW, flag_method_type, whitenoise_parameter[1])
montecalrosimulation(units, winds, NG, NW, 1, whitenoise_parameter)


montecalrosimulation(units, winds, NG, NW, flag_method_type, whitenoise_parameter)

unitsamplestatues = rand(0:1, NG, 10)
whitenoise_parameter
Sampling_Statue = unitsamplestatues[:, 1]
f_nadir, t_nadir, H, δp, Kg, Fg, Rg, Dg = creatingfrequencyresponsesamplingdata(
	units, winds, NW, NG, Sampling_Statue, 2, 1, whitenoise_parameter,
)

sample_Num = 200
horizon = 30
flag, symflag = 1, 0
seed = rand()
δf_positor, δf_actual, δf_samplieddata = simulate(generate_data, particle_filter, sample_Num, horizon, flag, symflag, units, winds, Sampling_Statue, whitenoise_parameter, rand())
f_nadir = maximum(abs.(δf_positor)) * (-1)
t_nadir = findmax(abs.(δf_positor))[1]

# -----
# FIXME redebug
Sampling_Statue[1, 1] = 1
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
current_Dw = sum(winds.Dw .* adjustablewindsVSCpower) / sum(winds.p_max) # Dw
current_Mw = sum(winds.Mw .* adjustablewindsVSCpower) / sum(winds.p_max) # Mw
current_Hw = current_Mw / 2
current_Rw =
	1 / (sum(winds.Kw .* inverse_winds_Rw .* (ones(NW, 1) - winds.Fcmode) .* winds.p_max) /
		 sum((ones(NW, 1) .* winds.p_max)))

# units parameters
adjustabletheramlpower = units.p_max .* Sampling_Statue
current_Kg = 1.0
current_Tg = mean(units.Tg)
current_Fg_div_Rg =
	sum(units.Kg .* units.Fg ./ units.Rg .* adjustabletheramlpower) / sum(units.p_max)
current_Rg = 1 / (sum(units.Kg ./ units.Rg .* adjustabletheramlpower) / sum(units.p_max)) # Kg
current_Fg = current_Fg_div_Rg * current_Rg
current_Dg = sum(units.Dg .* adjustabletheramlpower) / sum(units.p_max)
current_Hg = sum(units.Hg .* adjustabletheramlpower) / sum(units.p_max)
current_Mg = current_Hg * 2

#  powers for intia frequency response
localapparentpower = (sum(units.p_max[:, 1]) + sum(winds.p_max .* winds.Fcmode))
sumapparentpower   = (localapparentpower - sum(winds.p_max .* winds.Fcmode) + sum(winds.p_max))
p_step             = maximum(units.p_max) * 1.0

# sumD and sumH
if flag == 1
	current_sumD = (sum(units.Dg .* adjustabletheramlpower)) / sumapparentpower
	current_sumH = (sum(current_Mg .* adjustabletheramlpower)) / sumapparentpower / 2
else
	current_sumD =
		(
			sum(units.Dg .* adjustabletheramlpower) + sum(
				winds.Dw .* winds.p_max .* winds.Fcmode +
				winds.Kw .* winds.p_max .* (ones(NW, 1) - winds.Fcmode),
			) +
			1 / 0.4
		) / sumapparentpower
	current_sumH =
		(
			sum(current_Mg .* adjustabletheramlpower) +
			sum(current_Mw .* adjustablewindsVSCpower)
		) / sumapparentpower / 2
end

D, H, F, R, T, K, δp = current_sumD, current_sumH, current_Fg, current_Rg, current_Tg, current_Kg, p_step

# FIXME - sfr -----
# @show [D, H, F, T, R, δp]
f_nadir, t_nadir = calculation_frequencynadirdata(H, R, D, F, T, δp)
f_nadir = f_nadir * (-1)


#TODO - PF SFR
sample_Num = 200
horizon = 30
flag, symflag = 1, 0
seed = rand()
δf_positor, δf_actual, δf_samplieddata = simulate(generate_data, particle_filter, sample_Num, horizon, flag, symflag, units, winds, Sampling_Statue, whitenoise_parameter, seed)
f_nadir = maximum(abs.(δf_positor)) * (-1)
t_nadir = findmax(abs.(δf_positor))[1]

data = generate_data(initialize, trans_fun_1, obs_fun_1, 50, flag, symflag, units, winds, Sampling_Statue, whitenoise_parameter)

Mg, Hg, Dg, Tg, Rg, Fg, Kg, δp, endtime = formparameter(units, winds, Sampling_Statue, flag)
tt= 50
δf, δp_add = fill(0.0, tt), fill(0.0, tt)
# FIXME - the kernel debug
for i in 2:tt
    δfₜ_cur = trans_fun_1(trans_fun_0, δf[i-1, 1], δp_add[i-1, 1], δt, i, flag, units, winds, Sampling_Statue, whitenoise_parameter, 1)
    δf[i] = δfₜ_cur[1]
    δpₜ_add = obs_fun_1(obs_fun_0, δf[1:i, 1], i, flag, symflag, units, winds, Sampling_Statue, whitenoise_parameter, 1)
    δp_add[i] = δpₜ_add[1]
end

i = 2
δfₜ_cur = trans_fun_1(trans_fun_0, δf[i-1, 1], δp_add[i-1, 1], δt, i, flag, units, winds, Sampling_Statue, whitenoise_parameter, 1)


tem = trans_fun_0(δf[i-1, 1], δp_add[i-1,1], δt, i, flag, units, winds, Sampling_Statue, whitenoise_parameter)
tem = tem[1] .+ rand(Normal(0, Q), length(𝚫f_pre)) * whitenoise_parameter * (5^(2 - 1))


function trans_fun_0(δfₜ, δp_add, δt, i, flag, units, winds, Sampling_Statue, whitenoise_parameter)
	Mg, Hg, Dg, Tg, Rg, Fg, Kg, δp₀, endtime = formparameter(units, winds, Sampling_Statue, flag)
	# current mismatch power
	temp_δpₜ = δp₀ - Dg * δfₜ[1] - δp_add
	# current ROCOF
	grad_δfₜ = temp_δpₜ / (2 * Hg)
	# current frequency change
	temp_δfₜ = δfₜ[1] + grad_δfₜ * δt
end