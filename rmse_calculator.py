import os
from Bio.PDB import PDBParser
import pandas as pd
import numpy as np

parser = PDBParser(QUIET=True)

experiment_directory = "relabeled/"
generated_directory = "3dStructs/"

results = []

for file in os.listdir(generated_directory):
    if not file.endswith(".pdb"):
        continue

    file_id = file.split("_")[0]
    experiment_file = os.path.join(experiment_directory, file)
    generated_file = os.path.join(generated_directory, file)

    if not os.path.exists(experiment_file):
        continue

    experiment_structure = parser.get_structure("experiment", experiment_file)
    generated_structure = parser.get_structure("generated", generated_file)

    experiment_atoms = [atom.get_coord() for atom in experiment_structure.get_atoms()]
    generated_atoms = [atom.get_coord() for atom in generated_structure.get_atoms()]

    if len(experiment_atoms) != len(generated_atoms):
        rmse = "Atom mismatch"
    else:
        experiment_values = np.array(experiment_atoms)
        generated_values = np.array(generated_atoms)

        rmse = np.sqrt(np.mean(((generated_values - experiment_values) ** 2).sum(axis=1)))


    results.append({
        "id": file_id,
        "rmse": rmse
    })

df = pd.DataFrame(results)
df.to_csv("final_rmse.csv", index=False)