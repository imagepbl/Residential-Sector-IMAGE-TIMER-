# **************************************************************************************************
#  Module Energy Use for Residential Appliances

#  Original Indian Version
#  Author:   Bas van Ruijven
#  Date:     May 2008

#  Global Version By: Vassilis Daioglou
#  Date:     May 2010
#  Modification for python: Efstratios Mikropoulos
# **************************************************************************************************

# imports
import pandas as pd
import numpy as np
import pym

NREGIONS = 27
TURQ = 13
YEARS = 130


def turq1_to_3(variable):
    # This function reculculates the values of TURQ 1, 2 and 3, by adding the appropriate sub-TURQ.
    last_axis = variable.ndim - 1

    variable[..., 1] = variable[..., 3:8].sum(axis=last_axis)
    variable[..., 2] = variable[..., 8:13].sum(axis=last_axis)
    variable[..., 0] = variable[..., 1] + variable[..., 2]
    return variable


class BuildingStock():

    def __init__(self, floorspace,  build_lifetime):

        self.floorspace = floorspace
        self.build_lifetime = build_lifetime
        self.increased_floor = np.zeros((YEARS, NREGIONS, TURQ))
        self.decreased_floor = np.zeros((YEARS, NREGIONS, TURQ))
        self.new_buildings = np.zeros((YEARS, NREGIONS, TURQ))
        self.decommissioned = np.zeros((YEARS, NREGIONS, TURQ))
        self.decommission_initial = np.zeros((NREGIONS, TURQ))
        self.abandoned_buildings = np.zeros((YEARS, NREGIONS, TURQ))
        self.abandoned_residual = np.zeros((YEARS, NREGIONS, TURQ))
        self.stock = np.zeros((YEARS, NREGIONS, TURQ))

    def floorspace_change(self):
        floor_change = np.diff(self.floorspace, axis=0, prepend=0)
        floor_change[0, :] = 0
        self.increased_floor = floor_change.copy()
        self.increased_floor[self.increased_floor < 0] = 0
        self.decreased_floor = floor_change.copy()
        self.decreased_floor[self.decreased_floor > 0] = 0
        self.decreased_floor = -self.decreased_floor
        self.increased_floor = turq1_to_3(self.increased_floor)
        self.decreased_floor = turq1_to_3(self.decreased_floor)

        return(self.increased_floor, self.decreased_floor)

    def decommission(self, t, initial_stock):

        decom_distribution = np.zeros((YEARS, NREGIONS, TURQ))
        decom_total = np.zeros((YEARS, NREGIONS, TURQ))

        self.decommission_initial = (
            initial_stock.T / self.build_lifetime).T
        self.decommission_initial = turq1_to_3(self.decommission_initial)

        for R in range(NREGIONS):
            if t > self.build_lifetime[R]:
                decom_distribution[t, R, ...] = self.new_buildings[int(
                    t - (self.build_lifetime[R])), R, ...]
            else:
                decom_distribution[t,
                                   R, ...] = self.decommission_initial[R, ...]

            decom_total[t, ...] = decom_distribution[t, ...] + \
                self.abandoned_buildings[t-1, ...] + \
                self.abandoned_residual[t-1, ...]

        self.decommissioned[t, ...] = decom_total[t, ...].copy()
        self.decommissioned[self.decommissioned < 0] = 0
        self.decommissioned = turq1_to_3(self.decommissioned)

        self.abandoned_residual[t, ...] = decom_total[t, ...].copy()
        self.abandoned_residual[self.abandoned_residual > 0] = 0
        self.abandoned_residual = turq1_to_3(self.abandoned_residual)

        return (self.decommission_initial, self.decommissioned, self.abandoned_residual)

    def stock_change(self, t):
        change = np.zeros((YEARS, NREGIONS, TURQ))

        change[t, ...] = self.increased_floor[t, ...] - \
            self.decreased_floor[t, ...] + self.decommissioned[t, ...]

        self.new_buildings[t, ...] = change[t, ...].copy()
        self.new_buildings[self.new_buildings < 0] = 0
        self.new_buildings = turq1_to_3(self.new_buildings)

        self.abandoned_buildings[t, ...] = change[t, ...].copy()
        self.abandoned_buildings[self.abandoned_buildings > 0] = 0
        self.abandoned_buildings = turq1_to_3(self.abandoned_buildings)

        return (self.new_buildings, self.abandoned_buildings)

    def calculate_stock(self, t, initial_stock):
        if t == 0:
            self.stock[t, ...] = initial_stock
        else:
            self.stock[t, ...] = self.stock[t-1, ...] + self.new_buildings[t, ...] + \
                self.abandoned_buildings[t, ...] - self.decommissioned[t, ...]
        return self.stock


# ---------------------------------------------------------------------------- #
if __name__ == "__main__":

    # Import variables
    floorspace, time = pym.read_mym("floorspace.out", "Buildings/imports")
    population = pym.read_mym("pop_q.out", "Buildings/imports")[0]
    build_lifetime = pym.read_mym("Dwel_LT.dat", "Buildings/imports")[0]

    initial_stock = floorspace[0, :, :]
    print(build_lifetime.shape)
    print(floorspace)
    # Definitions
    # increased_floor = np.zeros((YEARS, NREGIONS, TURQ))
    # decreased_floor = np.zeros((YEARS, NREGIONS, TURQ))

    buildings = BuildingStock(floorspace, build_lifetime)

    for t in range(YEARS):
        buildings.floorspace_change()
        decommission_initial, decommissioned, abandoned_residual = buildings.decommission(
            t, initial_stock)
        new_buildings, abandoned_buildings = buildings.stock_change(t)
        building_stock = buildings.calculate_stock(t, initial_stock)

    print("done")
