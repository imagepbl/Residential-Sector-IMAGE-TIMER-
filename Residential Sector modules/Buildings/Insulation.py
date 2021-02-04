# imports
import pandas as pd
import numpy as np
import pym
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from building_stock import BuildingStock

NREGIONS = 27
TURQ = 13
YEARS = 130
TSCEN = 49  # (2020)
YEAR_IN = 1971
CONSTRUCTION = 100
INSULATION_LT = 30
# (The calibration year refers to U-values and the costs of U-Values. Literature values used are from 2002 (YEARS = 31).
CALIBRATION_YEAR = 31
SURFACES = 4
LEVELS = 6
FUEL_TYPES = 8
APPS = 11


def turq1_to_3(variable):
    # This function reculculates the values of TURQ 1, 2 and 3, by adding the appropriate sub-TURQ.
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


def market_share_multinomial_logit(cost_per_level):

    level_axis = 3
    logit_par = 2

    cost_floor = np.amin(cost_per_level, axis=level_axis)
    cost_floor[cost_floor <= 0] = 0.0000001

    cost_floor = dimension_modify(
        cost_per_level, cost_floor, cost_per_level)[1]

    relative_cost = cost_per_level / cost_floor

    share_nominator = np.exp(-logit_par * relative_cost)
    share_denominator = (np.exp(-logit_par * relative_cost)).sum(axis=3)

    share_denominator = dimension_modify(
        share_nominator, share_denominator, share_nominator)[1]

    share = share_nominator / share_denominator

    return share


class InsulationInBuildings():

    def __init__(self, floorspace,  build_lifetime, household_floor, UValue, cost_UValue, fuel_price_heating, fuel_price_appliances, space_heating, appiances_end_use, cdr_I):

        # PI32=numpy.float32(PI)
        # , dtype=np.float64
        self.floorspace = floorspace
        self.floorspace = np.float32(self.floorspace)
        self.build_lifetime = build_lifetime
        self.build_lifetime = np.float32(self.build_lifetime)
        self.household_floor = household_floor
        self.household_floor = np.float32(self.household_floor)
        self.UValue = UValue
        self.UValue = np.float32(self.UValue)
        self.cost_UValue = cost_UValue
        self.cost_UValue = np.float32(self.cost_UValue)
        self.fuel_price_heating = fuel_price_heating
        self.fuel_price_heating = np.float32(self.fuel_price_heating)
        self.fuel_price_appliances = fuel_price_appliances
        self.fuel_price_appliances = np.float32(self.fuel_price_appliances)
        self.space_heating = space_heating
        self.space_heating = np.float32(self.space_heating)
        self.appiances_end_use = appiances_end_use
        self.appiances_end_use = np.float32(self.appiances_end_use)
        # self.new_buildings = new_buildings
        # self.building_stock = building_stock
        self.cdr_I = cdr_I
        # Subsidy for insulation appliance (%). Initially put as zero.
        self.insulation_subsidy = 0
        self.learning_factor = 1
        self.insulation_new_buildings_stock = np.zeros((YEARS, NREGIONS, TURQ))
        self.insulation_new_buildings_stock = np.float32(
            self.insulation_new_buildings_stock)
        self.renovated_total_stock = np.zeros((YEARS, NREGIONS, TURQ))
        self.renovated_total_stock = np.float32(self.renovated_total_stock)
        self.surface_area = np.zeros((YEARS, NREGIONS, TURQ, SURFACES))
        self.surface_area = np.float32(self.surface_area)
        self.heating_fuel_cost_per_level = np.zeros(
            (YEARS, NREGIONS, TURQ, LEVELS))
        self.heating_fuel_cost_per_level = np.float32(
            self.heating_fuel_cost_per_level)
        self.cooling_fuel_cost_per_level = np.zeros(
            (YEARS, NREGIONS, TURQ, LEVELS))
        self.cooling_fuel_cost_per_level = np.float32(
            self.cooling_fuel_cost_per_level)
        self.useful_energy_intensity = np.ones((YEARS, NREGIONS, TURQ))
        self.useful_energy_intensity = np.float32(self.useful_energy_intensity)
        self.useful_energy_intensity_per_level = np.ones(
            (YEARS, NREGIONS, TURQ, LEVELS))
        self.useful_energy_intensity_per_level = np.float32(
            self.useful_energy_intensity_per_level)
        self.suitable_buildings = np.zeros(
            (YEARS, NREGIONS, TURQ, CONSTRUCTION))
        self.suitable_buildings = np.float32(self.suitable_buildings)
        self.suitable_buildings_add = np.zeros((YEARS, NREGIONS, TURQ))
        self.suitable_buildings_add = np.float32(self.suitable_buildings_add)
        self.renovated_stock = np.zeros(
            (YEARS, NREGIONS, TURQ, LEVELS, CONSTRUCTION))
        self.renovated_stock = np.float32(self.renovated_stock)
        # self.renovated_stock = np.float32(self.renovated_stock)
        self.renovated_add_stock = np.zeros((YEARS, NREGIONS, TURQ, LEVELS))
        self.renovated_add_stock = np.float32(self.renovated_add_stock)
        self.netcost_per_level_and_constr_year = np.zeros(
            (YEARS, NREGIONS, TURQ, LEVELS, CONSTRUCTION))
        self.netcost_per_level_and_constr_year = np.float32(
            self.netcost_per_level_and_constr_year)
        self.renovated = np.zeros(
            (YEARS, NREGIONS, TURQ, LEVELS, CONSTRUCTION))
        self.renovated = np.float32(self.renovated)
        self.insulation_LT_renovated = np.zeros(
            (YEARS, NREGIONS, CONSTRUCTION))
        self.insulation_LT_renovated = np.float32(self.insulation_LT_renovated)
        self.capital_recovery_factor_renovated = np.zeros(
            (YEARS, NREGIONS, TURQ, CONSTRUCTION))
        self.capital_recovery_factor_renovated = np.float32(
            self.capital_recovery_factor_renovated)
        self.capital_cost_renovated = np.zeros(
            (YEARS, NREGIONS, TURQ, LEVELS, CONSTRUCTION))
        self.capital_cost_renovated = np.float32(self.capital_cost_renovated)

    def learning(self, t):
        # Learning rate of improvement of the building's shell (insulation). Learning rate is applied from 2002 onward, since values of this year were found in literature.
        # For the learning process, the whole word is taken in account. So, capacities refer to added global capacities.

        EXPERIENCE_INDEX = -0.2

        # This is equal to 31, which corresponds to year 2002
        if t >= CALIBRATION_YEAR:
            capacity_in = self.renovated_total_stock[30, :, 3:TURQ].sum(
            ) + self.insulation_new_buildings_stock[30, :, 3:TURQ].sum()
        else:
            capacity_in = 0

        capacity_fin = self.renovated_total_stock[t, :, 3:13].sum(
        ) + self.insulation_new_buildings_stock[t, :, 3:13].sum()

        if t > CALIBRATION_YEAR:
            self.learning_factor = (
                capacity_fin/capacity_in)**(EXPERIENCE_INDEX)
        else:
            self.learning_factor = 1.0

        return self.learning_factor

    def building_dimensions(self):

        # Description of a prototype house. These dimensions, analogies and properties are used to describe the average global building shape.
        # Inner surfaces are not accounted for insulation appliance (internal walls, floors between 2 appartments, etc.).
        # Roofs that are in contact with the outside environment (%).
        EXTERNAL_ROOF_PERC = 0.5
        # Floors that are in contact with soil (%).
        EXTERNAL_FLOOR_PERC = 0.5
        # Walls that are in contact with the outside environment (%).
        EXTERNAL_WALL_PERC = 0.75
        # Buildings average height. Assumed to be 3 meters.
        BUILDINGS_HEIGHT = 3
        # Buildings are assumed to be rectangular. The biger side is 1.6 times the smaller side (relative size).
        BUILDINGS_LARGE_SIDE = 1.6
        # 20% of the outter (external) walls are assumed to be covered by windows.
        WINDOWS_PERC = 0.2

        floors_to_floorspace = np.zeros((YEARS, NREGIONS, TURQ))
        roofs_to_floorspace = np.zeros((YEARS, NREGIONS, TURQ))
        walls_to_floorspace = np.zeros((YEARS, NREGIONS, TURQ))
        windows_to_floorspace = np.zeros((YEARS, NREGIONS, TURQ))

        floors_to_floorspace[...] = EXTERNAL_FLOOR_PERC

        roofs_to_floorspace[...] = EXTERNAL_ROOF_PERC

        walls_per_house = EXTERNAL_WALL_PERC * (1-WINDOWS_PERC) * ((2 * BUILDINGS_HEIGHT * np.sqrt(self.household_floor / BUILDINGS_LARGE_SIDE)) + (
            2 * BUILDINGS_HEIGHT * BUILDINGS_LARGE_SIDE * np.sqrt(self.household_floor / BUILDINGS_LARGE_SIDE)))

        walls_to_floorspace = division_with_exception_at_zero(
            walls_per_house, self.household_floor)

        windows_per_house = (walls_to_floorspace /
                             (1 - WINDOWS_PERC) * WINDOWS_PERC * self.household_floor)

        windows_to_floorspace = division_with_exception_at_zero(
            windows_per_house, self.household_floor)

        # Total surface areas (m^2)
        self.surface_area[..., 0] = walls_to_floorspace * self.floorspace
        self.surface_area[..., 1] = roofs_to_floorspace * self.floorspace
        self.surface_area[..., 2] = floors_to_floorspace * self.floorspace
        self.surface_area[..., 3] = windows_to_floorspace * self.floorspace
        self.surface_area = turq1_to_3(self.surface_area)

        return self.surface_area

    def useful_energy_intensity_and_investment_cost_per_level(self, t):
        # For SSP2 the autonomus energy efficiency indicator is put as 0.0075.
        AEEI = 0.0075

        self.UValue[t, ...] = self.UValue[t, ...] * \
            pow((1-AEEI), (t-CALIBRATION_YEAR))

        self.cost_UValue[t, ...] = self.cost_UValue[t, ...] * \
            self.learning_factor

        prototype_variable = np.zeros(
            (YEARS, NREGIONS, TURQ, SURFACES, LEVELS))
        UValue_mod, surface_area_mod = dimension_modify(
            self.UValue, self.surface_area, prototype_variable)
        # thermal conductance's unit is W/K
        thermal_conductance_sum = (UValue_mod * surface_area_mod).sum(axis=3)

        cost_UValue_mod = dimension_modify(
            self.cost_UValue, self.surface_area, prototype_variable)[0]

        # surfaces costs unit is $
        cost_all_surfaces = (cost_UValue_mod * surface_area_mod).sum(axis=3)

        prototype_variable = np.zeros((YEARS, NREGIONS, TURQ, LEVELS))
        thermal_conductance_sum, floorspace_mod = dimension_modify(
            thermal_conductance_sum, self.floorspace, prototype_variable)
        cost_all_surfaces = dimension_modify(
            cost_all_surfaces, self.floorspace, prototype_variable)[0]

        q_coefficient = division_with_exception_at_zero(
            thermal_conductance_sum, floorspace_mod)

        # Conversion of Q_coef to Usuful Energy Intensity. W/m^2/K to KJ/m^2/HDD.
        Watt_per_Kelvin_to_KJ_per_HHD = 86.4
        self.useful_energy_intensity_per_level = q_coefficient * \
            Watt_per_Kelvin_to_KJ_per_HHD

        self.investment_cost_per_level = division_with_exception_at_zero(
            cost_all_surfaces, floorspace_mod)

        return self.useful_energy_intensity_per_level, self.investment_cost_per_level

    def fuel_cost_buildings(self, t, building_stock):
        heating_fuel_cost = np.zeros((YEARS, NREGIONS, TURQ))

        fuel_price_heating, space_heating = dimension_modify(
            self.fuel_price_heating, self.space_heating, self.space_heating)

        heating_fuel_cost[t, ...] = (
            fuel_price_heating[t, ...] * space_heating[t-1, ...]).sum(axis=2)/building_stock[t, ...]

        useful_energy_intensity_mod = dimension_modify(
            self.useful_energy_intensity_per_level, self.useful_energy_intensity, self.useful_energy_intensity_per_level)[1]

        heating_fuel_cost = dimension_modify(
            self.useful_energy_intensity_per_level, heating_fuel_cost, self.useful_energy_intensity_per_level)[1]

        self.heating_fuel_cost_per_level[t, ...] = heating_fuel_cost[t, ...] * self.useful_energy_intensity_per_level[t, ...] / \
            useful_energy_intensity_mod[t-1, ...]

        # Cooling fuel Cost

        cooling_appliances_end_use = self.appiances_end_use[:, :, :, 1:3].sum(
            axis=3)

        cooling_fuel_cost = np.zeros((YEARS, NREGIONS, TURQ))
        wh_to_MJ = 0.0036
        cooling_fuel_cost[t, ...] = (self.fuel_price_appliances[t, ...] *
                                     cooling_appliances_end_use[t-1, ...]) * wh_to_MJ * population[t, ...] / building_stock[t, ...]

        cooling_fuel_cost = dimension_modify(
            self.useful_energy_intensity_per_level, cooling_fuel_cost, self.useful_energy_intensity_per_level)[1]

        self.cooling_fuel_cost_per_level[t, ...] = cooling_fuel_cost[t, ...] * \
            self.useful_energy_intensity_per_level[t, ...] / \
            useful_energy_intensity_mod[t-1, ...]

        return self.heating_fuel_cost_per_level, self.cooling_fuel_cost_per_level

    def insulation_net_cost_new_buildings(self, t):
        # cdr_I=Discount rate. The dwelling lifetime is biger than the insulation/renovation lifetime(INSULATION_LT) for all regions. Not the case for the renovations later on the model.
        capital_recovery_factor_denominator = (
            1-pow((1+cdr_I), -INSULATION_LT))

        capital_recovery_factor = division_with_exception_at_zero(
            cdr_I, capital_recovery_factor_denominator)

        capital_recovery_factor[capital_recovery_factor == 0] = 0.1

        capital_recovery_factor = dimension_modify(
            self.investment_cost_per_level, capital_recovery_factor, self.investment_cost_per_level)[1]

        # Capital recovery reduction factor represents an easier decision of people to invest in insulation in new buildings.
        REDUCTION_FACTOR = 0.9

        capital_cost_per_level = self.investment_cost_per_level * \
            (1 - self.insulation_subsidy) * \
            capital_recovery_factor * REDUCTION_FACTOR

        # Premium Factor ! Premium factors -> 0 as CTax -> 1000 AND LEVELS -> 6
        premium_factor = np.zeros((YEARS, NREGIONS, TURQ, LEVELS))
        # premium_factor = factor_insulation * \
        #     (1 - (insulation_discount) * insulation_tax_discount)

        # premium_factor =PremFacInsulation[R, j] * (1 - (InsulPremDiscount(j) * InsulCTaxDiscount(CarbonTax[R]))),
        # R=1 TO NR27, j=1 TO REN;

        self.netcost_per_level = capital_cost_per_level + \
            self.heating_fuel_cost_per_level + \
            self.cooling_fuel_cost_per_level + premium_factor

        return self.netcost_per_level

    def insulation_new_buildings_per_level(self, t, new_buildings):

        self.insulation_level_ms = market_share_multinomial_logit(
            self.netcost_per_level)

        new_buildings_mod = dimension_modify(
            self.insulation_level_ms, new_buildings, self.insulation_level_ms)[1]

        self.new_buildings_per_level = self.insulation_level_ms * new_buildings_mod
        self.new_buildings_per_level = turq1_to_3(self.new_buildings_per_level)

        self.insulation_investments_new_buildings = (
            self.new_buildings_per_level * self.investment_cost_per_level).sum(axis=3)/pow(10, 9)

        self.insulation_new_buildings_stock = self.insulation_new_buildings_stock[t -
                                                                                  1, ...] + self.new_buildings_per_level[..., 1: LEVELS].sum(axis=3)

        return self.insulation_level_ms, self.new_buildings_per_level, self.insulation_investments_new_buildings, self.insulation_new_buildings_stock

    def buidings_suitable_for_renovation(self, t, new_buildings):
        # Suitable for renovation are the buildings that were not built within the last INSULATION_LT (30) years (not 'new buildings') and the buildings that have not already been renovated in the last INSULATION_LT years.
        # The year limits  t > (Y+1970)+INSULATION_LT AND t < (Y+1970)+build_lifetime are taken mainly to reduce calculation time. Buildings outside this timelin can not renovate (too old or too new).

        # First get fraction of building deemed appropriate for renavation based on thier income
        # This is an ad-hoc behavioural parameter assuming that poorer households are less likely to chose to renovate (note: CDR used later only affects the renovation level, NOT the willingness to renovate)

        # multiplier was a concept added by Vassilis that after long discussions it will be removed
        multiplier = 1
        # if t <= TSCEN:
        #     multiplier = min(1, Suit_Renov_Mult_in(
        #         PCOpcT_ppp[R, i]) * RenovCTaxPremium(CarbonTax[R]))
        # else:
        #     multiplier = (Convergence * MIN(1, Suit_Renov_Mult_scen(PCOpcT_ppp[R, i]) * RenovCTaxPremium(
        #         CarbonTax[R]))) + ((1-Convergence)*multiplier[R, i](tscen))

        # The second dimension for time (years) represent the costruction year of the buildings. The first of course represent the year of the run.

        for construct_year in range(min(t, 100)):
            for R in range(NREGIONS):
                if (construct_year + INSULATION_LT) < t < (construct_year + self.build_lifetime[R]):
                    # The sum makes sure that the buildings constructed in a specific year can be renovated once.
                    self.suitable_buildings[t, :, :, construct_year] = multiplier * new_buildings[construct_year, ...] - \
                        self.renovated_stock[t-1, :, :, 1: 6,
                                             construct_year].sum(axis=2)
                else:
                    self.suitable_buildings[t, :, :, construct_year] = 0
        self.suitable_buildings[self.suitable_buildings < 0] = 0

        # Suitable buildings for renovation, built before 1971
        # Additional suitable floorspace/buildings for renovation. For the years before INSULATION_LT+1971 (usually 2001) the model cannot see the suitable buildings, as it cannot look so far in the past(values before 1971).
        for R in range(NREGIONS):
            if t <= INSULATION_LT:
                self.suitable_buildings_add[t, R, ...] = multiplier * \
                    ((self.build_lifetime[R] - (INSULATION_LT+1))*decommission_initial[R, :] -
                        self.renovated_add_stock[t-1, R, :, 1:6].sum(axis=1))
            else:
                self.suitable_buildings_add[t, R, ...] = multiplier * (self.build_lifetime[R]-(INSULATION_LT+1)) - (t - (
                    INSULATION_LT + 1971)) * decommission_initial[R, :] - self.renovated_add_stock[t-1, R, :, 1:6].sum(axis=1)
        # This part adds the residences build before 1971 only when it is necessary
        self.suitable_buildings_add[self.suitable_buildings_add < 0] = 0

        return self.suitable_buildings, self.suitable_buildings_add

    def renovation_costs_per_level_and_construction_year(self, t):
        for R in range(NREGIONS):
            for construct_year in range(min(t, 100)):
                if t < (construct_year + self.build_lifetime[R]):
                    self.insulation_LT_renovated[t, R, construct_year] = min(
                        INSULATION_LT, construct_year + self.build_lifetime[R] - t)
                else:
                    self.insulation_LT_renovated[t, R, construct_year] = 0

        prototype = np.zeros(
            (YEARS, NREGIONS, TURQ, CONSTRUCTION), dtype=np.float32)
        cdr_I_mod, insulation_LT_renovated_mod = dimension_modify(
            cdr_I, self.insulation_LT_renovated, prototype)

        capital_recovery_factor_denominator = 1 - \
            pow(1+cdr_I_mod, -insulation_LT_renovated_mod)

        cdr_I_mod2 = np.zeros_like(capital_recovery_factor_denominator)
        cdr_I_mod2[..., :] = cdr_I_mod[...]
        cdr_I_mod2[..., :1] = cdr_I_mod

        for R in range(NREGIONS):
            for construct_year in range(min(t, 100)):
                if t < (construct_year + self.build_lifetime[R]):
                    self.capital_recovery_factor_renovated[t, R, ...] = division_with_exception_at_zero(
                        cdr_I_mod2[t, R, ...], capital_recovery_factor_denominator[t, R, ...])
                else:
                    self.capital_recovery_factor_renovated[t, R, ...] = 0.1
        self.capital_recovery_factor_renovated[self.capital_recovery_factor_renovated == 0] = 0.1

        prototype = np.zeros(
            (YEARS, NREGIONS, TURQ, LEVELS, CONSTRUCTION), dtype=np.float32)
        investment_cost_mod, capital_recovery_factor_renovated_mod = dimension_modify(
            self.investment_cost_per_level, self.capital_recovery_factor_renovated, prototype)

        for R in range(NREGIONS):
            for construct_year in range(min(t, 100)):
                if t < (construct_year + self.build_lifetime[R]):
                    self.capital_cost_renovated[t, ...] = investment_cost_mod[t, ...] * (
                        1 - self.insulation_subsidy) * capital_recovery_factor_renovated_mod[t, ...]
                else:
                    self.capital_cost_renovated[t, ...] = 0

        premium_factor = np.zeros(
            (YEARS, NREGIONS, TURQ, LEVELS, CONSTRUCTION), dtype=np.float32)
        heating_fuel_cost_per_level_mod = dimension_modify(
            self.heating_fuel_cost_per_level, self.capital_cost_renovated, self.capital_cost_renovated)[0]
        cooling_fuel_cost_per_level_mod = dimension_modify(
            self.cooling_fuel_cost_per_level, self.capital_cost_renovated, self.capital_cost_renovated)[0]

        for R in range(NREGIONS):
            for construct_year in range(min(t, 100)):
                if (construct_year + INSULATION_LT) < t < (construct_year + self.build_lifetime[R]) or construct_year == 0:
                    self.netcost_per_level_and_constr_year[t, ...] = self.capital_cost_renovated[t, ...] + \
                        heating_fuel_cost_per_level_mod[t, ...] + \
                        cooling_fuel_cost_per_level_mod[t, ...] + \
                        premium_factor[t, ...]
                else:
                    self.netcost_per_level_and_constr_year[t, ...] = 0

        return self.netcost_per_level_and_constr_year

    def renovated_buildings(self, t):

        market_share_renovated = market_share_multinomial_logit(
            self.netcost_per_level_and_constr_year)

        suitable_buildings_mod = dimension_modify(
            self.suitable_buildings, market_share_renovated, market_share_renovated)[0]

        for R in range(NREGIONS):
            for construct_year in range(min(t, 100)):
                if (construct_year + INSULATION_LT) < t < (construct_year + self.build_lifetime[R]):
                    self.renovated[t, R, ...] = market_share_renovated[t,
                                                                       R, ...] * suitable_buildings_mod[t, R, ...]
                else:
                    self.renovated[t, R, ...] = 0

        self.renovated = turq1_to_3(self.renovated)

        prototype = np.zeros((YEARS, NREGIONS, TURQ, LEVELS), dtype=np.float32)

        suitable_buildings_add_mod = dimension_modify(
            self.suitable_buildings_add, prototype, prototype)[0]

        self.renovated_add = market_share_renovated[...,
                                                    0] * suitable_buildings_add_mod

        # if t < 2: #(2=1973?)
        #     renovated_add = market_share_renovated[...,0] * suitable_buildings_add
        # elif t = 2:
        #     renovated_add = market_share_renovated[1,:,:,:,0] * suitable_buildings_add
        # else:
        #     renovated_add = market_share_renovated[1,:,:,:,0] * suitable_buildings_add

        self.renovated_add = turq1_to_3(self.renovated_add)

        # For level = 1 (j=1) we actually have non-renovated buildings.
        self.renovated_per_level = self.renovated.sum(
            axis=4) + self.renovated_add

        return self.renovated, self.renovated_add, self.renovated_per_level

    def renovation_stocks(self, t, initial_stock):

        self.renovated_stock[t, ...] = self.renovated_stock[t -
                                                            1, ...] + self.renovated[t, ...]

        self.renovated_add_stock[t, ...] = self.renovated_add_stock[t -
                                                                    1, ...] + self.renovated_add[t, ...]

        self.renovated_total = self.renovated_per_level[:, :, :, 1:LEVELS].sum(
            axis=3) + self.renovated_add[:, :, :, 1:LEVELS].sum(axis=3)

        if t == 0:
            self.renovated_total_stock[0, ...] = initial_stock * 0.5/100
        else:
            self.renovated_total_stock[t, ...] = self.renovated_total_stock[t -
                                                                            1, ...] + self.renovated_total[t, ...]

        return self.renovated_stock, self.renovated_add_stock, self.renovated_total, self.renovated_total_stock

    def renovation_rate(self, t, building_stock):

        self.rate = division_with_exception_at_zero(
            self.renovated_total, building_stock)

        if t <= TSCEN:
            self.rate_average = self.rate
        else:
            self.rate_average = np.average(self.rate[0:t-TSCEN, ...], axis=0)

        if t == TSCEN+1:
            self.rate_average_historic = np.average(
                self.rate[0:(TSCEN+1-YEAR_IN), ...], axis=0)
        else:
            self.rate_average_historic = 0

        return self.rate, self.rate_average, self.rate_average_historic


# ---------------------------------------------------------------------------- #
if __name__ == "__main__":

    # Import variables
    pym.read_mym("floorspace.out")
    pym.read_mym("Dwel_LT.dat")
    pym.read_mym("pop_q.out")
    pym.read_mym("U_Values.dat")
    pym.read_mym("Costs_UValues.dat")
    floorspace = pym.read_mym("floorspace.out")[0]
    population = pym.read_mym("pop_q.out")[0]
    HDD = pym.read_mym("HDD.out")[0]
    HDD_annual = HDD[:, :, 12]
    household_floor = pym.read_mym("HHFloor.out")[0]
    u_value_help = pym.read_mym("U_Values.dat")
    u_value = np.reshape(u_value_help, ((5, 4, 6)))
    cost_u_value_help = pym.read_mym("Costs_UValues.dat")
    cost_u_value = np.reshape(cost_u_value_help, ((5, 4, 6)))
    lifetime_help = pym.read_mym("Dwel_LT.dat")
    build_lifetime = lifetime_help[0, :]
    cdr_I = pym.read_mym("CDR_I.out")[0]
    fuel_price_heating = pym.read_mym("FuelPrice_heating.out")[0]
    fuel_price_appliances = pym.read_mym("FuelPrice_Final.out")[0]
    space_heating = pym.read_mym("SpaceHeatingTOTAL_TURQ.out")[0]
    appiances_end_use = pym.read_mym("ApplEnUse.out")[0]

    fuel_price_heating = fuel_price_heating[..., 0:8]

    # Modifications
    initial_stock = floorspace[0, :, :]
    # U-values of Europe for 2002 are found and used in this model. The values are provided for 3 climate zones, Warm, Cold, and Moderate depended on the areas HDD.
    # Climate zones are linked to a specific HDD (heating degree day) value. For more precesion an interpolation through the climate zones takes place improving the U-Values for each region. Same idea is also applied to the costs.
    points = np.array((7000, 4500, 3000, 1800, 0))
    UValue = griddata(points, u_value, HDD_annual, method='linear')
    cost_UValue = griddata(points, cost_u_value, HDD_annual, method='linear')

    # Definitions

    buildings = BuildingStock(floorspace, build_lifetime)
    insulations = InsulationInBuildings(floorspace,  build_lifetime, household_floor, UValue, cost_UValue, fuel_price_heating,
                                        fuel_price_appliances, space_heating, appiances_end_use, cdr_I)

    for t in range(YEARS):
        print(t)
        # Import variables from the Class BuildingStock
        buildings.floorspace_change()
        decommission_initial, decommissioned, abandoned_residual = buildings.decommission(
            t, initial_stock)
        new_buildings, abandoned_buildings = buildings.stock_change(t)
        building_stock = buildings.calculate_stock(t, initial_stock)

        # Insulation module
        insulations.learning(t)
        insulations.building_dimensions()
        insulations.useful_energy_intensity_and_investment_cost_per_level(t)
        insulations.fuel_cost_buildings(t, building_stock)
        insulations.insulation_net_cost_new_buildings(t)
        insulations.insulation_new_buildings_per_level(t, new_buildings)
        insulations.buidings_suitable_for_renovation(t, new_buildings)
        insulations.renovation_costs_per_level_and_construction_year(t)
        insulations.renovated_buildings(t)
        insulations.renovation_stocks(t, initial_stock)
        rate = insulations.renovation_rate(t, building_stock)[0]

    print("END")

    print(rate[0:3, 10])
