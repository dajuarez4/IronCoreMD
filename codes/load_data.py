import numpy as np
import json

data = np.load("/scratch/dajuarez4/testing/non-magnetic/2.41_5000K.npz", allow_pickle=False)

print(data.files)

initial_positions_alat = data["initial_positions_alat"]
initial_cell_alat = data["initial_cell_alat"]

positions = data["positions"]                    # QE units, usually crystal
positions_unit = data["positions_unit"]

cell_parameters = data["cell_parameters"]        # QE units if printed
cell_parameters_unit = data["cell_parameters_unit"]
species = data["species"]
forces_ry_au = data["forces_ry_au"]
energy_ry = data["energy_ry"]
temperature_K = data["temperature_K"]
pressure_kbar = data["pressure_kbar"]

print(initial_positions_alat.shape)   # (natoms, 3)
print(initial_cell_alat.shape)        # (3, 3)
print(positions.shape)                # (nsteps, natoms, 3)
print(positions_unit[:5])
# print(initial_positions_alat)
# print(cell_parameters_unit[:5])
# print(cell_parameters_unit)
# print(data["input_file"])
# print(data["input_cell_parameters"])
# print(energy_ry)
# print(temperature_K)
# print(np.mean(energy_ry[0:399]))
# print(np.mean(temperature_K))
# print(pressure_kbar)
# print(np.mean(pressure_kbar[0:399])*0.1)  # convert to GPa
# print(len(species))