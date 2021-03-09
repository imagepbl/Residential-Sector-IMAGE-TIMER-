# !********************************************************************
# ! Module Energy Use for Residential Appliances
# !
# ! Original Indian Version
# ! Author:   Bas van Ruijven
# ! Date:     May 2008
# !
# ! Global Version By: Vassilis Daioglou
# ! Date:     May 2010
# !********************************************************************

# Standard library imports
import numpy as np
import pandas as pd
import os.path
from scipy.interpolate import griddata

# Local application imports
# import pym

NREGIONS = 27
TURQ = 13
TUR = 3
YEARS = 130
NECN = 8
TSCEN = 49  # (2020)
# 1=fan, 2= Air Cooler, 3=Air Conditioner, 4=Refrigerator, 5=Microwave, 6=Washing Machine, 7=Clothes Dryer, 8=Dish Washer, 9=TV, 10=VCR/DVD, 11=PC/Other
# Four major end-uses: 1) Cooling (1-3) ; 2) Food storage and processing (4-5) ; 3) Washing & Cleaning (6-8) ; 4) Entertainment (9-11)
APPLIANCES = 11
APPLIANCE_LIFETIME = 15


def turq1_to_3(variable):
    # This function reculculates the values of TURQ 1, 2 and 3, by adding the appropriate sub-TURQs.
    turq_axis = 2

    if variable.ndim == 3:
        variable[..., 1] = variable[..., 3:8].sum(axis=turq_axis)
        variable[..., 2] = variable[..., 8:13].sum(axis=turq_axis)
        variable[..., 0] = variable[..., 1] + variable[..., 2]

    elif variable.ndim > 3:
        variable[:, :, 1, ...] = variable[:, :, 3:8, ...].sum(axis=turq_axis)
        variable[:, :, 2, ...] = variable[:, :, 8:13, ...].sum(axis=turq_axis)
        variable[:, :, 0, ...] = variable[:,
                                          :, 1, ...] + variable[:, :, 2, ...]

    return variable


def dimension_modify(variable1, variable2, prototype_variable):
    # This function modifies the variables adding dimensions where necessary. This is done to achieve a proper broadcasting later in the code.
    # prototype_variable represents the shape that is decired.

    min_dim = min(variable1.ndim, variable2.ndim)
    i = 0
    while i < min_dim:
        if prototype_variable.shape[i] != variable1.shape[i]:
            variable1 = np.expand_dims(variable1, axis=i)
        if prototype_variable.shape[i] != variable2.shape[i]:
            variable2 = np.expand_dims(variable2, axis=i)

        min_dim = min(variable1.ndim, variable2.ndim)
        i += 1

    while variable1.ndim != variable2.ndim:
        if variable1.ndim > variable2.ndim:
            variable2 = np.expand_dims(variable2, axis=-1)
        else:
            variable1 = np.expand_dims(variable1, axis=-1)

    return variable1, variable2


def division_with_exception_at_zero(variable1, variable2):
    # This function is a division of two arrays at which if an element of the denominator is zero it returns zero.

    division_result = np.divide(variable1, variable2, out=np.zeros_like(
        variable1.astype(float)), where=variable2 != 0)

    return division_result


def logarith_with_exception_at_zero(variable):
    # this function returns 0 as a result when the log requirments are not fulfilled (i.e., variable <=0)
    logarith = np.log(variable, out=np.zeros_like(
        variable), where=(variable != 0))
    return logarith


def capital_initial(variable1, variable2):
    # This function is a division of two arrays at which if an element of the denominator is zero it returns zero.

    division_result = np.divide(variable1, variable2, out=np.zeros_like(
        variable1.astype(float)), where=variable2 != 0)

    return division_result


def double_exponential_function(variable1, variable2, variable3):
    # This function is repeated.
    function_result = np.exp(- variable1 *
                             np.exp(- (variable2/1000) * variable3))

    return function_result


class ResidentialAppliances():

    def __init__(self, floorspace,  build_lifetime, CDD, household_expenditures_per_capita, households):

        self.floorspace = floorspace
        self.build_lifetime = build_lifetime
        self.population = population
        self.CDD = CDD
        self.household_expenditures_per_capita = household_expenditures_per_capita
        self.households = households
        self.carbon_tax = carbon_tax

        self.annuity_factor = np.zeros((YEARS, NREGIONS, TURQ))
        self.fuel_price_final = np.zeros((YEARS, NREGIONS, TURQ))
        self.levelised_cost_electricity_kWh1 = np.zeros((YEARS, NREGIONS))
        self.levelised_cost_electricity_kWh = np.zeros((YEARS, NREGIONS))
        self.efficiency_cooling = np.zeros((YEARS, NREGIONS))
        # self.climate_max_saturation = np.zeros((YEARS, NREGIONS))
        self.saturation_price = np.zeros((YEARS, NREGIONS))
        self.saturation = np.zeros((YEARS, NREGIONS, 3))
        self.availability = np.zeros((YEARS, NREGIONS, TURQ, 3))
        # Diffusion in TURQ for Fans, air coolers and air conditioners
        self.diffusion_cooling = np.zeros((YEARS, NREGIONS, TURQ, 3))
        # Diffusion in TUR for rest
        self.diffusion_rest = np.zeros((YEARS, NREGIONS, TUR, APPLIANCES))
        self.rough_diffusion = np.zeros((YEARS, NREGIONS, TUR, APPLIANCES))
        self.rough_diffusion_income_delay = np.zeros(
            (YEARS, NREGIONS, TUR, APPLIANCES))
        # Diffusion in Q for rest
        self.diffusion_quintile = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))
        self.diffusion_total = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))
        self.diffusion_adjusted = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))
        self.diffution_average = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))
        self.capital_requirement = np.zeros(
            (YEARS, NREGIONS, TURQ, APPLIANCES))
        self.capital_initial = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))
        self.capital_total = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))
        self.capital_new = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))
        self.capital_depreciation = np.zeros(
            (YEARS, NREGIONS, TURQ, APPLIANCES))
        self.new_capital_smoothing = np.zeros((YEARS, NREGIONS, APPLIANCES))
        self.capital_total_per_household = np.zeros(
            (YEARS, NREGIONS, TURQ, APPLIANCES))
        self.cooling_appliances_capital_per_household = np.zeros(
            (YEARS, NREGIONS, TURQ))
        self.unit_energy_consumption_min = np.zeros(
            (YEARS, NREGIONS, APPLIANCES))
        self.unit_energy_consumption = np.zeros(
            (YEARS, NREGIONS, TURQ, APPLIANCES))
        self.unit_energy_consumption_base1 = np.zeros(
            (YEARS, NREGIONS, APPLIANCES))
        self.unit_energy_consumption_base = np.zeros(
            (YEARS, NREGIONS, APPLIANCES))
        self.unit_energy_consumption_unit_new = np.zeros(
            (YEARS, NREGIONS, TURQ, APPLIANCES))
        # kWh/L after coe effects accounted for, baseline coe (levelised cost)
        self.unit_energy_consumption_unit_new1 = np.zeros(
            (YEARS, NREGIONS, TURQ, APPLIANCES))
        # kWh/L after coe effects accounted for, changed coe (levelised cost)
        self.unit_energy_consumption_unit_new2 = np.zeros(
            (YEARS, NREGIONS, TURQ, APPLIANCES))
        self.unit_energy_consumption_difference = np.zeros(
            (YEARS, NREGIONS, TURQ, APPLIANCES))
        # kWh/unit after coe effects accounted for
        self.unit_energy_consumption_new = np.zeros(
            (YEARS, NREGIONS, TURQ, APPLIANCES))
        # Multiplier for energy savings assuming appliances have no stand-by
        self.standby_factor = np.zeros((YEARS, NREGIONS, APPLIANCES))
        # Multiplier for energy savings assuming changes in appliance use(more efficient use)
        self.efficiency_factor = np.zeros((YEARS, NREGIONS, APPLIANCES))
        # Energy Requirement of all appliences, tmin
        self.energy_requirment = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))
        # Energy requirement of all appliances, all other t
        self.end_use_flow = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))
        # Energy use of all new(marginal) appliances at given t
        self.end_use_new = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))
        self.end_use_initial = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))
        # Energy use of all depreciated appliances at given t
        self.end_use_depreciation = np.zeros(
            (YEARS, NREGIONS, TURQ, APPLIANCES))
        # Energy requirement of all appliances, corrected for improvements in insulation(affects air conditioners/coorlers)
        self.end_use_total = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))
        # Energy requirement of all appliances
        self.end_use_final = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))
        # Aggregate UEC when energy standards are introduced
        self.unit_energy_consumption_aggregated = np.zeros(
            (YEARS, NREGIONS, TURQ, APPLIANCES))
        self.other = np.zeros((YEARS, NREGIONS, TURQ))
        self.other_share = np.zeros((YEARS, NREGIONS, TURQ))
        self.end_use_per_capita = np.zeros((YEARS, NREGIONS))
        self.end_use1 = np.zeros((YEARS, NREGIONS, TURQ))
        self.end_use3 = np.zeros((YEARS, NREGIONS, TURQ))
        self.end_use4 = np.zeros((YEARS, NREGIONS, TURQ))
        self.energy_aggregated = np.zeros((YEARS, NREGIONS, TUR))
        # GJ For ALL Appliances
        self.energy_aggregated2 = np.zeros((YEARS, NREGIONS, TURQ, NECN))
        # GJ for cooling appliances
        self.energy_aggregated3 = np.zeros((YEARS, NREGIONS, TURQ, NECN))
        # GJ for rest of appliances
        self.energy_aggregated4 = np.zeros((YEARS, NREGIONS, TURQ, NECN))
        # Aggregate improvement in UEC compared to tscen
        self.unit_energy_consumption_improvement = np.zeros((YEARS, NREGIONS))
        self.energy_aggregated3_per_capita = np.zeros((YEARS, NREGIONS, NECN))
        self.energy_aggregated4_per_capita = np.zeros((YEARS, NREGIONS, NECN))
        self.end_use_total = np.zeros((YEARS, NREGIONS, NECN))
        self.end_use = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))
        self.end_use_tuss = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES+1))
        self.fuel_price_final = np.zeros((YEARS, NREGIONS, TURQ))

    def economic_parameters(self, t, consumer_discount_rate, consumer_discount_rate_extra, price_secondary_fuel, fuel_subsidy):

        consumer_discount_rate_total = consumer_discount_rate + consumer_discount_rate_extra
        for R in range(NREGIONS):
            for i in range(TURQ):
                if consumer_discount_rate_total[t, R, i] != 0:
                    self.annuity_factor[t, R, i] = consumer_discount_rate_total[t, R, i] / (
                        1 - (1 + consumer_discount_rate_total[t, R, i])) ** (-APPLIANCE_LIFETIME)
                else:
                    self.annuity_factor[t, R, i] = 0

        # 2005$/ kWh, baseline
        self.levelised_cost_electricity_kWh1 = price_secondary_fuel * 1.252083 * 0.0036
        # 2005$/ kWh, future
        self.levelised_cost_electricity_kWh = price_secondary_fuel * 1.252083 * 0.0036

        # self.fuel_price_final = price_secondary_fuel - fuel_subsidy[..., 7]

        return self.annuity_factor, self.levelised_cost_electricity_kWh1, self.levelised_cost_electricity_kWh, self.fuel_price_final

    def difussion_cooling_appliances(self, t, saturation_price_remi, air_cooler_fraction, beta, gamma):

        # calculation for air conditioning from Morna/Detlef
        climate_max_saturation = 1.00 - 0.949 * \
            np.exp(-0.00187 * self.CDD[..., -1])

        # Diffusion of Fans as function of Floorspace per capita (Variation on McNeil/Letschert 2007)
        # Historic saturation level is a function of price development(Based on 1993/2002 difference for India)

        self.saturation_price[t, 2:10] = saturation_price_remi[t]
        self.saturation_price[t, 12:22] = saturation_price_remi[t]
        self.saturation_price[t, 24:NREGIONS] = saturation_price_remi[t]

        # saturation for 1-fans, 2-air coolers, 3-air conditioners
        self.saturation[..., 0] = self.saturation_price * \
            climate_max_saturation
        self.saturation[..., 1] = air_cooler_fraction * \
            self.saturation_price * climate_max_saturation
        self.saturation[..., 2] = (1-air_cooler_fraction) * \
            self.saturation_price * climate_max_saturation

        gamma_annual = np.zeros((YEARS, NREGIONS, TURQ, 7))
        betta_annual = np.zeros((YEARS, NREGIONS, TURQ, 7))

        household_expenditures_per_capita_mod, gamma_mod = dimension_modify(
            self.household_expenditures_per_capita, gamma_annual, gamma_annual)
        beta_mod = dimension_modify(
            betta_annual, gamma_annual, gamma_annual)[0]

        # Fitted diffusion levels of Fans, Air Coolers and Air Conditioners (TURQ)
        self.availability = double_exponential_function(
            beta_mod[..., :3], gamma_mod[..., :3], household_expenditures_per_capita_mod)
        for turq in range(3, 8):
            self.availability[:, :, turq, :] = double_exponential_function(
                beta_mod[..., 1, :3], gamma_mod[..., 1, :3], household_expenditures_per_capita_mod[:, :, turq, :])
        for turq in range(8, TURQ):
            self.availability[:, :, turq, :] = double_exponential_function(
                beta_mod[:, :, 2, :3], gamma_mod[:, :, 2, :3], household_expenditures_per_capita_mod[:, :, turq, :])

        # Fitted diffusion levels of Fans, Air Coolers and Air Conditioners (TURQ)

        saturation_mod = dimension_modify(
            self.saturation, self.availability, self.availability)[0]

        self.diffusion_cooling = saturation_mod * self.availability

        return self.diffusion_cooling

    def difussion_rest_appliances(self, t, phi1, phi2, phi3, income_delay, urbanization, appliance_a, household_expenditures_fraction_quintiles):

        # DIFFUSION FOR REST OF APPLIANCES FOLLOW 3 BASIC PATTERNS. In all cases Total averaged out from Urban and Rural.
        # 1. Constant Growth: Clothes Dryer, Dish Washer
        # 2. An 'income delay', income at which diffusion starts depends on time: Microwave, VCR, PC
        # 3. 'Surface', rate of diffusion and saturation point vary: Refrigerator, Washing Machine, TV
        # Appliances with income delay: !@1970 diffusion starts at 10000$(2005)ppp, !@2000 diffusion starts at 700$(2005)ppp.

        print(phi3.shape)
        print(phi2.shape)
        print(phi1.shape)

        self.rough_diffusion = phi1 * \
            double_exponential_function(
                phi2, phi3, self.household_expenditures_per_capita)
        self.rough_diffusion_income_delay = phi1 * \
            double_exponential_function(
                phi2, phi3, (self.household_expenditures_per_capita-income_delay))

        # FOOD PREPARATION AND STORAGE APPLIANCES
        # 3-Refrigerators, 4-Microwaves
        self.diffusion_rest[..., 1:TUR,
                            3] = self.rough_diffusion[..., 1:TUR, 3]
        self.diffusion_rest[..., 1:TUR,
                            4] = self.rough_diffusion_income_delay[..., 1:TUR, 4]
        # CLEANING APPLIANCES
        # 5-Washing Machine, 6-Clothes Dryer, 7-Dish Washer
        self.diffusion_rest[..., 1:TUR,
                            5:8] = self.rough_diffusion[..., 1:TUR, 5:8]
        # ENTERTAINMENT APPLIANCES
        # 8-Television, 9-DVD/VCR, 10-PC
        self.diffusion_rest[..., 1:TUR,
                            8] = self.rough_diffusion[..., 1:TUR, 8]
        self.diffusion_rest[..., 1:TUR,
                            9:APPLIANCES] = self.rough_diffusion_income_delay[..., 1:TUR, 9:APPLIANCES]
        # TOTAL allocation
        self.diffusion_rest[..., 0, 3:APPLIANCES] = (
            self.diffusion_rest[..., 2, 3:APPLIANCES] * urbanization) + (self.diffusion_rest[..., 3, 3:APPLIANCES] * (1-urbanization))
        # Quintile Allocation

        self.diffusion_quintile[..., 3:8, 3:APPLIANCES] = (appliance_a * logarith_with_exception_at_zero(
            abs(household_expenditures_fraction_quintiles))+1)*self.diffusion_rest[..., 1, 3:APPLIANCES]
        self.diffusion_quintile[..., 8:TURQ, 3:APPLIANCES] = (appliance_a * logarith_with_exception_at_zero(
            abs(household_expenditures_fraction_quintiles))+1)*self.diffusion_rest[..., 2, 3:APPLIANCES]

        return self.diffusion_quintile

    def difussion_indicators(self, t, electrification):

        self.diffusion_total[..., 0:3] = self.diffusion_cooling
        self.diffusion_total[..., 0:3,
                             3:APPLIANCES] = self.diffusion_rest[..., 0:3, 3:APPLIANCES]
        self.diffusion_total[..., 3:TURQ,
                             3:APPLIANCES] = self.diffusion_quintile[..., 3:TURQ, 3:APPLIANCES]

        # Adjust for Electrification

        self.diffusion_adjusted = self.diffusion_total * electrification

        if self.households[:, 0] > 0:
            self.diffution_average[..., 0, :] = (
                self.diffusion_adjusted[:, :, 3:TURQ, :]*self.households[:, :, 3:TURQ]).sum(2) / self.households[:, 0]
        else:
            self.diffution_average[..., 0, :] = 0

        if self.households[:, 1] > 0:
            self.diffution_average[..., 1, :] = (
                self.diffusion_adjusted[:, :, 3:8, :]*self.households[:, :, 3:8]).sum(axis=2) / self.households[:, 1]
        else:
            self.diffution_average[..., 1, :] = 0

        if self.households[:, 2] > 0:
            self.diffution_average[..., 2, :] = (
                self.diffusion_adjusted[:, :, 8:TURQ, :]*self.households[:, :, 8:TURQ]).sum(axis=2) / self.households[:, 2]
        else:
            self.diffution_average[..., 2, :] = 0

        return self.diffusion_total, self.diffusion_adjusted, self.diffution_average

    def vintage_capital_stock(self, t):
        # Vintage capital stock model
        households_mod = dimension_modify(
            self.households, self.diffusion_adjusted, self.diffusion_adjusted)[0]

        capital_factor = np.ones((YEARS, NREGIONS, TURQ, APPLIANCES))
        # Fans only
        capital_factor[..., 0] = 0.04 * self.floorspace

        self.capital_requirement = households_mod * \
            self.diffusion_adjusted * capital_factor

        self.capital_initial[t, ...] = self.capital_requirement[0, ...]
        # ApplCap_ini[R, i, j] = LAST(ApplCap_ini[R, i, j], ApplCapReq[R, i, j]), R = 1 to NRTim, i = 1 to TURQ, j = 1 to 11

        if t == 0:
            self.capital_new[t, ...] = max(
                0, self.capital_requirement[t, ...] + self.capital_depreciation[t, ...] - self.capital_initial[t, ...])
            self.capital_total[t, ...] = self.capital_initial[t, ...] - \
                self.capital_depreciation[t, ...] + \
                self.capital_new[t, ...]
        else:
            self.capital_new[t, ...] = max(
                0, self.capital_requirement[t, ...] + self.capital_depreciation[t, ...] - self.capital_total[t-1, ...])
            self.capital_total[t, ...] = self.capital_total[t-1, ...] - \
                self.capital_depreciation[t, ...] + \
                self.capital_new[t, ...]

        if t < APPLIANCE_LIFETIME - 2:
            self.capital_depreciation[t,
                                      ...] = self.capital_initial[t, ...] / APPLIANCE_LIFETIME
        elif t < APPLIANCE_LIFETIME + 2:
            self.capital_depreciation[t, ...] = 1/5 * self.capital_new[(t-APPLIANCE_LIFETIME-2):(t-APPLIANCE_LIFETIME+2), ...].sum(
                axis=0) + max(0, (APPLIANCE_LIFETIME+2-t)/5) * self.capital_initial[t, ...] / APPLIANCE_LIFETIME
        else:
            self.capital_depreciation[t, ...] = 1/5 * self.capital_new[(
                t-APPLIANCE_LIFETIME-2):(t-APPLIANCE_LIFETIME+2), ...].sum(axis=0)

        self.new_capital_smoothing = (self.capital_new[(
            t-5):(t-1), :, 3:TURQ, :].average(axis=0)).sum(axis=2)

        self.capital_total = turq1_to_3(
            self.capital_total)

        self.capital_total_per_household = division_with_exception_at_zero(
            self.capital_total, households_mod)

        # of air coolers and conditioners per household, for calibration purposes.
        self.cooling_appliances_capital_per_household = self.capital_total_per_household[..., 1:3].sum(
            axis=3)

    def energy_consumption(self, t, efficiency_cooling1, efficiency_cooling_future, cooling_efficiency_per_ctax, appliances_efficiency_per_ctax, ctax, unit_energy_consumption_min_in, unit_energy_consumption_min_scen, convergence2, unit_energy_consumption_alpha, unit_energy_consumption_beta, unit_energy_consumption_beta_future, elasticity_coefficient, af, unit_energy_consumption_unit, standby_factor_in, appliances_efficiency_factor_in, coefficient_other):
        # Energy Use ( in GJ/yr), UEC in kWh
        # Unit Energy Consumption for Fan, Air Cooler and Air Conditioner

        cooling_efficiency_ctax = griddata(
            ctax, cooling_efficiency_per_ctax, self.carbon_tax, method='linear')

        if t > TSCEN:
            self.efficiency_cooling = max(
                efficiency_cooling_future, cooling_efficiency_ctax)
        else:
            self.efficiency_cooling = efficiency_cooling1

        # Unit energy consumption for rest of appliances, autonomous decline over time(determine dby UEC_beta), described in normalised units:
        # UEC_beta is scenario dependent.
        # 4. Refrigerator: kWh/L(total volume)
        # 5. Microwave: kWh/unit
        # 6. Washing Machine: kWh/L(washing volume)
        # 7. Clothes Dryer: kWh/kg(load)
        # 8. Dish Washer: kWh/washing cycles(not cycles per year, but cycles per wash!!!)
        # 9/10/11. TV/DVD-VCR/PC: kWh/unit

        if t > TSCEN:
            self.unit_energy_consumption_min = convergence2 * unit_energy_consumption_min_scen + \
                ((1-convergence2) * unit_energy_consumption_min_in)
        else:
            self.unit_energy_consumption_min = unit_energy_consumption_min_in

        self.unit_energy_consumption_base1 = unit_energy_consumption_alpha * \
            unit_energy_consumption_beta**(t+1) + \
            self.unit_energy_consumption_min

        # These have to be set up in order to accomodate changing of UEC_beta, otherwise curve is discontinuous.
        beta_change1 = np.zeros((YEARS, NREGIONS, APPLIANCES))
        beta_change2 = np.zeros((YEARS, NREGIONS, APPLIANCES))

        appliances_efficiency_ctax = griddata(
            ctax, appliances_efficiency_per_ctax, self.carbon_tax, method='linear')

        if t <= TSCEN:
            self.unit_energy_consumption_base = self.unit_energy_consumption_base1
        else:
            beta_change1[t, ...] = self.unit_energy_consumption_base1[TSCEN, ...] - \
                self.unit_energy_consumption_min[t, ...]
            beta_change2 = (unit_energy_consumption_beta_future *
                            appliances_efficiency_ctax) ** (t-TSCEN)

            self.unit_energy_consumption_base = beta_change1 * \
                beta_change2 + self.unit_energy_consumption_min

        # Deviation from baseline due to coe.
        # Below variables irrelevant for cooling appliances

        prototype = np.zeros((YEARS, NREGIONS, TURQ, APPLIANCES))

        levelised_cost_electricity_kWh1_mod, levelised_cost_electricity_kWh_mod = dimension_modify(
            self.levelised_cost_electricity_kWh1, self.levelised_cost_electricity_kWh, prototype)

        cost_of_energy = np.array((
            levelised_cost_electricity_kWh1_mod, levelised_cost_electricity_kWh_mod))

        elasticity_coefficient_annuity = griddata(
            af, elasticity_coefficient, self.annuity_factor, method='linear')

        # look-up table from TIMER
        # UEC_unit_new1[R,i,j] 	= elast_coeff[1,j](AnnuityFactor[R,i]) * LOG(coe_kWh1[R]) + elast_coeff[2,j](AnnuityFactor[R,i]), R = 1 to NRTim, i = 1 to TURQ, j = 4 to 11; 	! kWh/unit, at baseline coe
        # UEC_unit_new2[R,i,j]	= elast_coeff[1,j](AnnuityFactor[R,i]) * LOG(coe_kWh[R]) + elast_coeff[2,j](AnnuityFactor[R,i]), R = 1 to NRTim, i = 1 to TURQ, j = 4 to 11;	! kWh/unit at new coe

        # kWh/unit, at baseline coe
        self.unit_energy_consumption_unit_new1[..., 3:APPLIANCES] = elasticity_coefficient_annuity[:, :, :, 0, 3:APPLIANCES] * np.log(
            cost_of_energy[0, ..., 3:APPLIANCES] + elasticity_coefficient_annuity[:, :, :, 1, 3:APPLIANCES])

        # kWh/unit at new coe
        self.unit_energy_consumption_unit_new2[..., 3:APPLIANCES] = elasticity_coefficient_annuity[:, :, :, 0, 3:APPLIANCES] * np.log(
            cost_of_energy[1, ..., 3:APPLIANCES] + elasticity_coefficient_annuity[:, :, :, 1, 3:APPLIANCES])

        self.unit_energy_consumption_difference[..., 3:APPLIANCES] = self.unit_energy_consumption_unit_new1[...,
                                                                                                            3:APPLIANCES] - self.unit_energy_consumption_unit_new2[..., 3:APPLIANCES]
        # We want the curve to go through UECbase in t = 2010
        self.unit_energy_consumption_unit_new[..., 3:APPLIANCES] = max(
            (self.unit_energy_consumption_base[..., 3:APPLIANCES] - self.unit_energy_consumption_difference[..., 3:APPLIANCES]), self.unit_energy_consumption_min[..., 3:APPLIANCES])

        unit_energy_consumption_unit_mod = dimension_modify(
            unit_energy_consumption_unit, prototype, prototype)

        # kWh/yr, of new refrigerators purchased
        self.unit_energy_consumption_new[..., 3:APPLIANCES] = self.unit_energy_consumption_unit_new[...,
                                                                                                    3:APPLIANCES] * unit_energy_consumption_unit_mod[..., 3:APPLIANCES]

        # Standyby and efficient appliance factors = 1 unless overridden from scenario

        FlagStandByMode = 1
        if FlagStandByMode == 1 and t > TSCEN:
            self.standby_factor = standby_factor_in
        else:
            self.standby_factor = 1

        # Final Unit Energy Consumption values:

        # fans
        self.unit_energy_consumption[...,
                                     0] = 0.0401 * self.CDD[..., -1] + 22.282

        # air cooler
        self.unit_energy_consumption[..., 1] = max(2.5/self.efficiency_cooling * CDD[..., -1] * (
            0.6053 * logarith_with_exception_at_zero(self.household_expenditures_per_capita) - 3.1897), 400) * 300/2160

        # air conditioner
        self.unit_energy_consumption[...,
                                     2] = self.unit_energy_consumption[..., 1] * 2160/300

        standby_factor_mod, appliances_efficiency_factor_mod = dimension_modify(
            self.standby_factor, appliances_efficiency_factor_in, prototype)

        self.unit_energy_consumption[..., 3:APPLIANCES] = self.unit_energy_consumption_new[..., 3:APPLIANCES] * \
            standby_factor_mod[..., 3:APPLIANCES] * \
            appliances_efficiency_factor_mod[..., 3:APPLIANCES]

        # Determine improvements in efficiency, to be used by service sector
        if t >= TSCEN:
            self.unit_energy_consumption_improvement = ((self.unit_energy_consumption[:, :, 1, 3] + self.unit_energy_consumption[:, :, 1, 5]) / 2) / (
                (self.unit_energy_consumption[TSCEN, :, 1, 3] + self.unit_energy_consumption[TSCEN, :, 1, 5] / 2))
        else:
            self.unit_energy_consumption_improvement = 1

        # Other unnacounted for appliance energy use, result from calibration to IEA data (30 years of Energy Use in IEA countries).
        # Related to income.
        # 'Other' Calculated on a per capita basis, since IEA data is on a per capita basis.
        # The coefficients of the "other" are based on regions which seem to have very high, when compared to default REMG, miscellaneous applies (Canada, USA, Japan)
        # and regions with low miscelaneous appliances (Western Europe, australia, and ROW by default).

        # Extra kwh/cap per region, FOR TOTAL
        self.other = max(0, coefficient_other[..., 0] * logarith_with_exception_at_zero(
            abs(household_expenditures_per_capita)) - coefficient_other[..., 1])

    def appliances_energy_use(self, t, insulation_cooling_reduction, outage_factor, unit_energy_consumption_initial):
        # Since UEC has been calculated on a marginal basis, and energy stock model is necessary to get aggregate UECs over appliance lifetime

        # NOT for Fans, Air Coolers and Air Conditioners which are climate based
        self.energy_requirment = self.capital_requirement * self.unit_energy_consumption
        self.end_use_initial[t, ...] = self.energy_requirment[0, ...]
        self.end_use_new = self.capital_new * self.unit_energy_consumption
        if t < APPLIANCE_LIFETIME - 2:
            self.end_use_depreciation = self.capital_depreciation * \
                unit_energy_consumption_initial
        elif t < APPLIANCE_LIFETIME - 2:
            self.end_use_depreciation = 1/5 * self.end_use_new[(t-APPLIANCE_LIFETIME-2):(t-APPLIANCE_LIFETIME+2), ...].sum(
                axis=0) + max(0, (APPLIANCE_LIFETIME+2-t)/5) * self.capital_depreciation[t, ...] * unit_energy_consumption_initial
        else:
            self.end_use_depreciation = 1/5 * \
                self.end_use_new[(t-APPLIANCE_LIFETIME-2)                                 :(t-APPLIANCE_LIFETIME+2), ...].sum(axis=0)
        self.end_use_flow[0, ...] = unit_energy_consumption_initial[0, ...] - \
            self.end_use_depreciation[0, ...] + self.end_use_new[0, ...]
        self.end_use_flow[t, ...] = self.end_use_flow[t-1, ...] - \
            self.end_use_depreciation[t, ...] + self.end_use_new[t, ...]

        # Make corrections of cooling appliances when there is insulation(via renovation). This has to be done separately in order to not interfere with the energy stock model
        insulation_reduction_mod = dimension_modify(
            insulation_reduction, self.end_use_total, self.end_use_total)[0]

        self.end_use_total = self.end_use_flow
        self.end_use_total[..., 1:2] = self.end_use_flow[...,
                                                         1:2] * insulation_reduction_mod

        self.end_use_total = turq1_to_3(self.end_use_total)

        outage_factor_mod, population_mod = dimension_modify(
            outage_factor, self.population, self.end_use_total)

        # kWh/cap
        self.end_use_per_capita = self.end_use_total * outage_factor / self.population
        self.end_use_per_capita[:, NREGIONS, ...] = 0

        self.unit_energy_consumption_aggregated = division_with_exception_at_zero(
            self.end_use_total, self.capital_total)
        self.unit_energy_consumption_aggregated[:, NREGIONS, ...] = 0

        # Electricity use for ALL appliances NOTE IT IS PER CAPITA!!!!!!!!!!!

        # kWh/cap, all appliances
        self.end_use1 = self.end_use.sum(axis=3) + self.other
        self.end_use3 = self.end_use[..., 0:3].sum(axis=3) + self.other
        self.end_use4 = self.end_use[..., 3:APPLIANCES].sum(
            axis=3) + self.other

        # output for TUSS

        # GJ
        self.end_use_tuss = self.end_use_total * outage_factor_mod * wh_to_MJ,

        outage_factor_mod = dimension_modify(
            outage_factor, self.population, self.population)[0]
        self.end_use_tuss[..., APPLIANCES+1] = self.other * \
            self.population * outage_factor_mod * wh_to_MJ
        self.end_use_tuss[..., APPLIANCES +
                          1] = turq1_to_3(self.end_use_tuss[..., APPLIANCES+1])

        self.other_share = max(0, division_with_exception_at_zero(
            self.other, self.end_use1)) * 100

        # GJ for appliances TUR
        self.energy_aggregated = turq1_to_3(
            self.end_use1 * self.population * wh_to_MJ)

        # ALL appliances - GJ
        # energy_aggregated2 equal to 0 for fuels/NECN =1 to 7.
        self.energy_aggregated2[..., -1] = self.end_use1 * \
            self.population * wh_to_MJ
        self.energy_aggregated2[:, NREGIONS, :, -1] = 0
        self.energy_aggregated2[..., -1] = turq1_to_3(
            self.energy_aggregated2[..., -1])

        # Cooling appliances - GJ
        # energy_aggregated3 equal to 0 for fuels/NECN =1 to 7.
        self.energy_aggregated3[..., -1] = self.end_use3 * \
            self.population * wh_to_MJ
        self.energy_aggregated3[..., -1] = turq1_to_3(
            self.energy_aggregated3[..., -1])

        self.energy_aggregated3[:, NREGIONS, :,
                                :] = self.energy_aggregated3[:, 0: NREGIONS, :, :].sum(axis=1)

        self.energy_aggregated3_per_capita = division_with_exception_at_zero(
            self.energy_aggregated3, self.population)

        # Rest of appliances
        # energy_aggregated4 equal to 0 for fuels/NECN =1 to 7.
        self.energy_aggregated4[..., -1] = self.end_use4 * \
            self.population * wh_to_MJ
        self.energy_aggregated4[..., -1] = turq1_to_3(
            self.energy_aggregated4[..., -1])

        self.energy_aggregated4[:, -1, :, -1] = 0

        self.energy_aggregated4_per_capita = division_with_exception_at_zero(
            self.energy_aggregated4, self.population)

        # In NECN for export - GJ
        self.end_use_final[:, :, -1] = self.energy_aggregated[:, :, 0]

        # Per capita for comparison - GJ/cap
        self.end_use_final_per_capita = division_with_exception_at_zero(
            self.energy_aggregated[:, :, 0], self.population[:, :, 0])

        # Expenditures for Appliance Energy
        self.energy_expenditures = self.fuel_price_final * \
            self.energy_aggregated2[..., -1]
        self.energy_expenditures[:, NREGIONS, :] = 0


if __name__ == "__main__":

    # REAL var1[NRC, 11](t)
    # !These have to be set up in order to accomodate changing of UEC_beta, otherwise curve is discontinuous.
    # REAL var2[NRC, 11](t)
    # ---------------------------------------------------------------------------- #
    # Import variables

    from appliances_data import *

    import_path = os.path.join("demand_functions", "appliance", "imports")

    # HHFloor[NR27,TURQ](t) - total floorspace
    floorspace, time = pym.read_mym("floorspace.out", "buildings/imports")
    population = pym.read_mym("pop_q.out", "buildings/imports")[0]
    build_lifetime = pym.read_mym("Dwel_LT.dat", "buildings/imports")[0]
    # CarbonTax[NR27](t) - From cookingwater.m
    carbon_tax = pym.read_mym("CarbonTax.out", import_path)[0]
    price_secondary_fuel = pym.read_mym("PriceSecFuel.out", import_path)[0]
    fuel_subsidy = pym.read_mym("Fuel_Subsidy.out", import_path)[0]
    convergence = pym.read_mym("convergence.out", import_path)[0]
    convergence2 = pym.read_mym("convergence2.out", import_path)[0]
    # PCOpcT_ppp[NR27,TURQ](t) - HHExp per capita for TURQ in pppUSD2005
    household_expenditures_per_capita = pym.read_mym(
        "PCOpcT_ppp.out", import_path)[0]
    consumer_discount_rate = pym.read_mym("CDR1.out", import_path)[0]
    households = pym.read_mym("Households.out", import_path)[0]
    # Electrification[NR27,TURQ](t) - From drivers/electrification.m
    electrification = pym.read_mym("Electrification.out", import_path)[0]
    # CDD[NR27,NMT](t) - from CoolHeat.m
    CDD = pym.read_mym("CDD.out", import_path)[0]
    # InsulationCoolingCor[NR27,TURQ](t) - Actual reduction in cooling demand (%) applied due to insulation of dwellings
    insulation_reduction = pym.read_mym(
        "InsulationCoolingCor.out", import_path)[0]
    urbanization = pym.read_mym("Urbanization.out", import_path)[0]
    # HHExpFac_ppp[NR27,TURQ](t) - from drivers.m - factor describingdescrbe differences between income quintile and average of urban or rural
    # (e.g. HHexp[R,4]/HHexp[R,2]).
    household_expenditures_fraction_quintiles = pym.read_mym(
        "HHExpFac_ppp.out", import_path)[0]
    outage_factor = pym.read_mym("OutageFactor.out", import_path)[0]
    payback_time_tc = pym.read_mym("PayBackTimeTC.out", import_path)[0]
    consumer_discount_rate_extra = 0
    #  MJperkWh/1000
    wh_to_MJ = 0.0036

    house_appliances = ResidentialAppliances(
        floorspace, build_lifetime, CDD, household_expenditures_per_capita, households)

    for t in range(YEARS):

        house_appliances.economic_parameters(
            t, consumer_discount_rate, consumer_discount_rate_extra, price_secondary_fuel, fuel_subsidy)
        house_appliances.difussion_cooling_appliances(
            t, saturation_price_remi, air_cooler_fraction, beta, gamma)
        # house_appliances.difussion_rest_appliances(
        #    t, phi1, phi2, phi3, income_delay, urbanization, appliance_a, household_expenditures_fraction_quintiles)
        # house_appliances.difussion_indicators(t, electrification)
        # print(t)

print("done")

print(population.shape)
