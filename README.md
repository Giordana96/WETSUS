# BufferModel
To check: excel sheet for the right sulphate concentration, CO2_void the K+ concentration is rounded to the right digit, sulphate concentration
Condition: [K+] = 1M, [S] = 0.1M, pCO2_an = 0.3bar, I_dens = 150A/m2
The standard unit for length is dm instead of m. 
Constants: q_area = Q/A = 9*10**-5dm/s (not fixed)
		   l_an = 0.5mm = 5*10**-3dm
		   l_m = 75um = 7.5*10**-4dm
		   l_cat = 0.5mm = 5*10**-3dm
		   l_tot = l_an + l_m + l_cat = 1.075*10**-2dm
		   mem_c = 2.67mol/dm3 (derived from IEC 0.8meq/g and water uptake of 30%)
		   Diffusion coefficients in solution: the values listed in BasicParSpeciesBook need to multiply by 10**-7 to get the unit dm2/s. 
		   Diffusion coefficients in membrane: 10 times lower than the values in solution. 