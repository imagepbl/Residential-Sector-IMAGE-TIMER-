# imports
import pandas as pd
import numpy as np
import pym

NREGIONS = 27
TURQ = 13
YEARS = 130


class ResidentialAppliances():

    def __init__(self, floorspace,  build_lifetime):

        self.floorspace = floorspace
        self.build_lifetime = build_lifetime
        self.x = np.zeros((YEARS, NREGIONS, TURQ))


# ---------------------------------------------------------------------------- #
if __name__ == "__main__":

    # !**************************************************************************************************
    # ! Module Energy Use for Residential Appliances
    # !
    # ! Original Indian Version
    # ! Author:   Bas van Ruijven
    # ! Date:     May 2008
    # !
    # ! Global Version By: Vassilis Daioglou
    # ! Date:     May 2010
    # !**************************************************************************************************

    # ! 1=fan, 2= Air Cooler, 3=Air Conditioner, 4=Refrigerator, 5=Microwave, 6=Washing Machine, 7=Clothes Dryer, 8=Dish Washer, 9=TV, 10=VCR/DVD, 11=PC/Other
    # ! Four major end-uses:
    # ! 1) Cooling (1-3)
    # ! 2) Food storage and processing (4-5)
    # ! 3) Washing & Cleaning (6-8)
    # ! 4) Entertainment (9-11)

    # INCLUDE ../Data/DatAppliance_TIMER.m	! Read datafiles

    # Import variables

    # HHFloor[NR27,TURQ](t) - total floorspace
    floorspace, time = pym.read_mym(
        "floorspace.out", "Buildings/imports")
    population = pym.read_mym("pop_q.out", "Buildings/imports")[0]
    build_lifetime = pym.read_mym("Dwel_LT.dat", "Buildings/imports")[0, :]
    # CarbonTax[NR27](t) - From cookingwater.m
    carbon_tax = pym.read_mym(
        "CarbonTax.out", "DemandFunctions/Appliance/imports")[0]
    price_secondary_fuel = pym.read_mym(
        "PriceSecFuel.out", "DemandFunctions/Appliance/imports")[0]
    fuel_subsidy = pym.read_mym(
        "Fuel_Subsidy.out", "DemandFunctions/Appliance/imports")[0]
    convergence = pym.read_mym(
        "convergence.out", "DemandFunctions/Appliance/imports")[0]
    convergence2 = pym.read_mym(
        "convergence2.out", "DemandFunctions/Appliance/imports")[0]
    # PCOpcT_ppp[NR27,TURQ](t) - HHExp per capita for TURQ in pppUSD2005
    household_expenditures_per_capita = pym.read_mym(
        "PCOpcT_ppp.out", "DemandFunctions/Appliance/imports")[0]
    consumer_discount_rate1 = pym.read_mym(
        "CDR1.out", "DemandFunctions/Appliance/imports")[0]
    households = pym.read_mym(
        "Households.out", "DemandFunctions/Appliance/imports")[0]
    # Electrification[NR27,TURQ](t) - From drivers/electrification.m
    electrification = pym.read_mym(
        "Electrification.out", "DemandFunctions/Appliance/imports")[0]
    # CDD[NR27,NMT](t) - from CoolHeat.m
    cdd = pym.read_mym("CDD.out", "DemandFunctions/Appliance/imports")[0]
    # InsulationCoolingCor[NR27,TURQ](t) - Actual reduction in cooling demand (%) applied due to insulation of dwellings
    insulation_cooling_reduction = pym.read_mym(
        "InsulationCoolingCor.out", "DemandFunctions/Appliance/imports")[0]
    urbanization = pym.read_mym(
        "Urbanization.out", "DemandFunctions/Appliance/imports")[0]
    # HHExpFac_ppp[NR27,TURQ](t) - from drivers.m - factor describingdescrbe differences between income quintile and average of urban or rural
    # (e.g. HHexp[R,4]/HHexp[R,2]).
    household_expenditures_fraction_for_quintiles = pym.read_mym(
        "HHExpFac_ppp.out", "DemandFunctions/Appliance/imports")[0]
    outage_factor = pym.read_mym(
        "OutageFactor.out", "DemandFunctions/Appliance/imports")[0]
    payback_time_tc = pym.read_mym(
        "PayBackTimeTC.out", "DemandFunctions/Appliance/imports")[0]

    appliances = ResidentialAppliances(floorspace, build_lifetime)

    for t in range(YEARS):
        #     buildings.floorspace_change()
        #     decommission_initial, decommissioned, abandoned_residual = buildings.decommission(
        #         t, initial_stock)
        #     new_buildings, abandoned_buildings = buildings.stock_change(t)
        #     building_stock = buildings.calculate_stock(t, initial_stock)
        print(t)

    print(convergence)
    print("done")

print("done")
