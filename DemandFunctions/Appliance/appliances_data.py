import pandas as pd
import numpy as np
import pym
import pym2
import os.path
from scipy.interpolate import griddata


def interpolation_function_1D(xp_data_points, interpolated_variable,  x_coordinates):
    # this function builds on the numpy 1D interpolation function
    # numpy.interp(x, xp, fp, left=None, right=None, period=None)
    fp_data_points = np.array(interpolated_variable[xp_data_points])
    y_coordinates = np.interp(x_coordinates, xp_data_points, fp_data_points)
    return y_coordinates


x = np.arange(1971, 2050)
data_path = os.path.join("demand_functions", "appliance", "imports", "data")

unit_energy_consumption_unit = pym.load_mym("UEC_unit.dat", data_path, x)

satur_price = pym.read_mym("Saturation.dat", data_path)[0]
# Fraction of cooling appliances which are air coolers
air_cooler_fraction = pym.read_mym("AirCoolFrac.dat", data_path)[0]
# Beta for cooling devices. Only 1-3 of 3rd dim used.
beta = pym.read_mym("Beta_global.dat", data_path)
# Gamma for cooling devices. Only 1-3 of 3rd dim used.
gamma = pym.read_mym("Gamma_global.dat", data_path)
# coefficient for quintile allocation of TV, washing machine and refrigerator
appliance_a = pym.read_mym("Appliance_a.dat", data_path)[0]
# Alpha coefficient for UEC development (kWh/normalising unit)
unit_energy_consumption_alpha = pym.read_mym("UEC_alpha.dat", data_path)[0]
# Beta coefficient (reduction rate) for UEC development (TV,washing,refrig)
unit_energy_consumption_beta = pym.read_mym("UEC_beta.dat", data_path)[0]
# Assumption on minimum UEC (kWh/normalising unit)
unit_energy_consumption_min_in = pym.read_mym("UEC_min.dat", data_path)[0]
unit_energy_consumption_min_scen = pym.read_mym("UEC_min.dat", data_path)[0]
# Normalising unit in which calculations are done (fridge L, washing machine L)
unit_energy_consumption_unit = pym.read_mym("UEC_unit.dat", data_path)[0]
# Normalising unit in which calculations are done (fridge L, washing machine L)
unit_energy_consumption_initial = pym.read_mym("UEC_ini.dat", data_path)[0]
# coefficients for elasticity equation, domain = annuity factor
# Phi1, Phi2, and Phi3 for gompertz for appliances 4-11
# x = np.arange(1971, 2101)
# y = pym.load_mym("Phi1.dat", "demand_functions/appliance/imports/data", x)
phi1 = pym.read_mym("Phi1.dat", data_path)[0]
phi2 = pym.read_mym("Phi2.dat", data_path)[0]
phi3 = pym.read_mym("Phi3.dat", data_path)[0]
elasticity_coefficient, af = pym.read_mym("elast_coeff.dat", data_path)
# For Microwaves, DVD/VCR and PC, time based income at which diffusion commences (@t=1970 = 10,000, @t=2000 = 700)
income_delay = pym.read_mym("IncomeDelay.dat", data_path)[0]
# coefficients for log curve describing miscellaneous appliances
coefficient_other = pym.read_mym("Other_coeff.dat", data_path)[0]
appliances_efficiency_per_ctax, ctax = pym.read_mym(
    "Appliance_eff.dat", data_path)

# GEA/FUTURE data

# Future scenario Beta coefficient
unit_energy_consumption_beta_future = pym.read_mym(
    "UEC_beta_fut.scn", data_path)[0]
# UEC standards in case of UEC policy
# unit_energy_consumption_policy = pym.read_mym("UECpolicy.dat", data_path)[0]
efficiency_cooling_future = pym.read_mym(
    "EffCooling_Global_fut.dat", data_path)[0]
# Energy efficiency of cooling equipment(numbers from MiniCAM input and LBL paper), Constant across regions
efficiency_cooling1 = pym.read_mym("EffCooling_Global.dat", data_path)[0]
# Ratio of ton clinker/ton cement
cooling_efficiency_per_ctax = pym.read_mym("EffCooling_Ctax.dat", data_path)[0]
# efficiency_per_ctax, ppsf = pym.read_mym("CTax_eff.dat", data_path)

# LIMITS LIFESTYLE DATA

standby_factor_in = pym.read_mym("StandByFac.scn", data_path)[0]
appliances_efficiency_factor_in = pym.read_mym("EffApplFac.scn", data_path)[0]

saturation_price_remi = np.ones(130)
saturation_price_remi[0:36] = satur_price[1:37, 0]

saturation_price_remi[36:49] = interpolation_function_1D(
    np.array([35, 49]), saturation_price_remi, np.arange(36, 49))

# # ***********************************************************************************************************
# annuity_factor = pym.read_mym("AnnuityFactor.out", data_path)[0]
# elasticity_coefficient_annuity = np.zeros((130, 27, 13, 2, 11))

# # for R in range(27):
# #     for i in range(13):
# #         elasticity_coefficient_annuity[:, R, i, ...] = griddata(
# #             af, elasticity_coefficient, annuity_factor[:, R, i], method='linear')


# elasticity_coefficient_annuity = griddata(
#     af, elasticity_coefficient, annuity_factor, method='linear')

# print(elasticity_coefficient_annuity[5, 25, 1, 0, 3])
# print(annuity_factor[5, 25, 1])
# print(elasticity_coefficient[:, 0, 3])
# print(elasticity_coefficient)
