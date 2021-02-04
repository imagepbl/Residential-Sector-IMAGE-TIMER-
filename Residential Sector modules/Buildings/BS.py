# imports
import pandas as pd
import numpy as np
import pym

NREGIONS = 27
TURQ = 13
YEARS = 130


def turq1_to_3(variable):
    # This function reculculates the values of TURQ 1, 2 and 3, by adding the appropriate sub-TURQ.
    if variable.shape == (130, 27, 13):
        variable[:, :, 1] = variable[:, :, 3:8].sum(axis=2)
        variable[:, :, 2] = variable[:, :, 8:13].sum(axis=2)
        variable[:, :, 0] = variable[:, :, 1]+variable[:, :, 2]

    elif variable.shape == (27, 13):
        variable[:, 1] = variable[:, 3:8].sum(axis=1)
        variable[:, 2] = variable[:, 8:13].sum(axis=1)
        variable[:, 0] = variable[:, 1]+variable[:, 1]
    return variable


class BuildingStock():

    def __init__(self, t, floorspace, build_lifetime):
        self.t = t
        self.floorspace = floorspace
        self.build_lifetime = build_lifetime
        self.new_buildings = new_buildings

    def floorspace_change(self):
        floor_change = np.diff(floorspace, axis=0, prepend=0)
        floor_change[0, :] = 0
        increased_floor = floor_change.copy()
        increased_floor[increased_floor < 0] = 0
        decreased_floor = floor_change.copy()
        decreased_floor[decreased_floor > 0] = 0
        decreased_floor = -decreased_floor
        self.increased_floor = increased_floor
        return(increased_floor, decreased_floor)

    def decommission(self):

        decom_distribution = np.zeros((YEARS, NREGIONS, TURQ))
        decom_help = np.zeros((YEARS, NREGIONS, TURQ))

        decommission_initial = (initial_building_stock.T / build_lifetime).T

        for R in range(NREGIONS):
            if t > build_lifetime[R]:
                decom_distribution[t, R, ...] = new_buildings[int(
                    t - (build_lifetime[R])), R, ...]
            else:
                decom_distribution[t, R, ...] = decommission_initial[R, ...]

            decom_help[t, ...] = decom_distribution[t, ...] + \
                abandoned_buildings[t-1, ...] + abandoned_residual[t-1, ...]

        decommissioned[t, ...] = decom_help[t, ...].copy()
        decommissioned[decommissioned < 0] = 0

        abandoned_residual[t, ...] = decom_help[t, ...].copy()
        abandoned_residual[abandoned_residual > 0] = 0

        return (decommission_initial, decommissioned, abandoned_residual)

    def stock_change(self):
        change = np.zeros((YEARS, NREGIONS, TURQ))

        change[t, ...] = increased_floor[t, ...] - \
            decreased_floor[t, ...] + decommissioned[t, ...]

        self.new_buildings[t, ...] = change[t, ...].copy()
        self.new_buildings[self.new_buildings < 0] = 0

        abandoned_buildings[t, ...] = change[t, ...].copy()
        abandoned_buildings[abandoned_buildings > 0] = 0

        return (self.new_buildings, abandoned_buildings)

    def calculate_stock(self):
        if t == 0:
            building_stock[t, ...] = initial_building_stock
        else:
            building_stock[t, ...] = building_stock[t-1, ...] + new_buildings[t, ...] + \
                abandoned_buildings[t, ...] - decommissioned[t, ...]
        return building_stock


# ---------------------------------------------------------------------------- #
if __name__ == "__main__":

    # Import variables
    pym.read_mym("floorspace.out")
    pym.read_mym("Dwel_LT.dat")
    pym.read_mym("pop_q.out")
    floorspace_per_capita, time = pym.read_mym("floorspace.out")
    population, time = pym.read_mym("pop_q.out")
    lifetime_help = pym.read_mym("Dwel_LT.dat")
    build_lifetime = lifetime_help[0, :]

    floorspace = floorspace_per_capita*population
    initial_building_stock = floorspace[0, :, :]

    # Definitions
    new_buildings = np.zeros((YEARS, NREGIONS, TURQ))
    decommissioned = np.zeros((YEARS, NREGIONS, TURQ))
    abandoned_buildings = np.zeros((YEARS, NREGIONS, TURQ))
    abandoned_residual = np.zeros((YEARS, NREGIONS, TURQ))
    building_stock = np.zeros((YEARS, NREGIONS, TURQ))

    for t in range(YEARS):
        BS = BuildingStock(t, floorspace, build_lifetime)
        BS.floorspace_change()
        increased_floor, decreased_floor = BS.floorspace_change()
        increased_floor = turq1_to_3(increased_floor)
        decreased_floor = turq1_to_3(decreased_floor)
        BS.decommission()
        decommission_initial, decommissioned, abandoned_residual = BS.decommission()
        decommission_initial = turq1_to_3(decommission_initial)
        decommissioned = turq1_to_3(decommissioned)
        BS.stock_change()
        new_buildings, abandoned_buildings = BS.stock_change()
        new_buildings = turq1_to_3(new_buildings)
        abandoned_buildings = turq1_to_3(abandoned_buildings)
        BS.calculate_stock()
        building_stock = BS.calculate_stock()
    print("done")

# ---------------------------------------------------------------------------- #
