
function coi_generation_frequencyparameters(units, winds, NW, NG, Sampling_Statue, flag)

	#NOTE - version - 1
	# normalized winds parameters through COI
	# vsmFC_number = sum(winds.Fcmode[:, 1])
	# doopFC_number = length(winds.Fcmode[:, 1]) - vsmFC_number
	# adjustablewindsVSCpower = winds.Fcmode .* winds.p_max
	# inverse_winds_Rw = zeros(NW, 1)
	# for i in 1:NW
	# 	if winds.Fcmode[i, 1] == 0
	# 		inverse_winds_Rw[i, 1] = 1 / winds.Rw[i, 1]
	# 		end
	# 	end
	# 	current_Kw = 1.0
	# 	current_Dw = sum(winds.Dw .* adjustablewindsVSCpower) / sum(adjustablewindsVSCpower) # Dw
	# 	current_Mw = sum(winds.Mw .* adjustablewindsVSCpower) / sum(adjustablewindsVSCpower) # Mw
	# 	current_Hw = current_Mw / 2
	# 	current_Rw = 1 / (sum(winds.Kw .* inverse_winds_Rw .* (ones(NW, 1) - winds.Fcmode) .*
	# 					  winds.p_max) /
	# 				  sum((ones(NW, 1) .* winds.p_max)))

	# 	# units parameters
	# 	adjustabletheramlpower = units.p_max .* Sampling_Statue
	# 	current_Kg = 1.0
	# 	current_Tg = mean(units.Tg)
	# 	current_Fg_div_Rg = sum(units.Kg .* units.Fg ./ units.Rg .* adjustabletheramlpower) /
	# 						sum(units.p_max)
	# 	current_Rg = 1 /
	# 				 (sum(units.Kg ./ units.Rg .* adjustabletheramlpower) / sum(units.p_max)) # Kg
	# 	current_Fg = current_Fg_div_Rg * current_Rg
	# 	current_Dg = sum(units.Dg .* adjustabletheramlpower) / sum(units.p_max)
	# 	current_Hg = sum(units.Hg .* adjustabletheramlpower) / sum(units.p_max)
	# 	current_Mg = current_Hg * 2

	# 	#  powers for intia frequency response
	# 	localapparentpower = (sum(units.p_max[:, 1]) + sum(winds.p_max .* winds.Fcmode))
	# 	sumapparentpower   = (localapparentpower - sum(winds.p_max .* winds.Fcmode) + sum(winds.p_max))
	# 	p_step             = maximum(units.p_max) * 1.0

	# 	# sumD and sumH
	# 	if flag == 1
	# 		# current_sumD = (sum(units.Dg .* adjustabletheramlpower)) / sumapparentpower
	# 		# current_sumH = (sum(current_Mg .* adjustabletheramlpower)) / sumapparentpower / 2
	# 		current_sumD = (sum(units.Dg .* adjustabletheramlpower) +
	# 						sum(winds.Dw .* winds.p_max .* winds.Fcmode +
	# 							winds.Kw .* winds.p_max .* (ones(NW, 1) - winds.Fcmode))) /
	# 					   sumapparentpower
	# 		current_sumH = (sum(current_Mg .* adjustabletheramlpower) +
	# 						sum(current_Mw .* adjustablewindsVSCpower)) /
	# 					   sumapparentpower / 2
	# 	else
	# 		current_sumD = (sum(units.Dg .* adjustabletheramlpower) +
	# 						sum(winds.Dw .* winds.p_max .* winds.Fcmode +
	# 							winds.Kw .* winds.p_max .* (ones(NW, 1) - winds.Fcmode)) +
	# 						1 / 0.4) / sumapparentpower
	# 		current_sumH = (sum(current_Mg .* adjustabletheramlpower) +
	# 						sum(current_Mw .* adjustablewindsVSCpower)) / sumapparentpower / 2
	# 	end

	#NOTE - version = 2

	#? normalized winds parameters through COI
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

	#? units parameters
	# Sampling_Statue = ones(NG, 1)
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

	if flag == 1
		current_sumD = (sum(units.Dg .* adjustabletheramlpower) +
						sum(winds.Dw .* winds.p_max .* winds.Fcmode +
							winds.Kw .* winds.p_max .* (ones(NW, 1) - winds.Fcmode)) * 20) /
					   (sum(winds.p_max) + sum(adjustabletheramlpower))
		current_sumH = (sum(current_Mg .* adjustabletheramlpower) +
						sum(current_Mw .* adjustablewindsVSCpower)) /
					   (sum(adjustablewindsVSCpower) + sum(adjustabletheramlpower)) / 2
	else
		current_sumD = (sum(units.Dg .* adjustabletheramlpower) +
						sum(winds.Dw .* winds.p_max .* winds.Fcmode +
							winds.Kw .* winds.p_max .* (ones(NW, 1) - winds.Fcmode)) +
						1 / 0.4) / sumapparentpower
		current_sumH = (sum(current_Mg .* adjustabletheramlpower) +
						sum(current_Mw .* adjustablewindsVSCpower)) / sumapparentpower / 2
	end

	D, H, F, R, T, K, δp = current_sumD,
	current_sumH, current_Fg, current_Rg, current_Tg, current_Kg, p_step

	return D, H, F, R, T, K, δp
end