import os
import requests
from Bio.PDB import PDBParser
import pandas as pd
import numpy as np

parser = PDBParser(QUIET=True)

generated_directory = "relaxed/"
experiment_directory = "ground_truths/"
os.makedirs(experiment_directory, exist_ok=True)

results = []

def rcsb_search(hla):
    if hla == "HLA-B0801":
        return "7NUI"
    if hla == "HLA-A0101":
        return "5BRZ"
    if hla == "HLA-A0201":
        return "7KGP"
    if hla == "HLA-A1101":
        return "8I5C"
    if hla == "HLA-B0702":
        return "7LGD"
    if hla == "HLA-B3501":
        return "8V50"

    return None
    

def download_pdb(pdb_id):
    ground_truth_path = os.path.join(experiment_directory, f"{pdb_id}.pdb")
    if os.path.exists(ground_truth_path):
        return ground_truth_path

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    request = requests.get(url)

    if request.status_code == 200:
        with open(ground_truth_path, "wb") as f:
            f.write(request.content)
        return ground_truth_path

    return None

for file in os.listdir(generated_directory):
    if not file.endswith(".pdb"):
        continue

    file_id = file.split("_")[0]
    file_HLA = file.split("_")[3]
    file_peptide = file.split("_")[4]
    generated_file = os.path.join(generated_directory, file)

    hla = file_HLA[:1] + "*" + file_HLA[1:3] + ":" + file_HLA[3:]
    hla = file_HLA.replace(" ", "")
    if hla.startswith("HLA-"):
        hla = hla
    else:
        hla = "HLA-" + hla

    experiment_pdb = rcsb_search(hla)

    if experiment_pdb is None:
        continue

    experiment_pdb_file = download_pdb(experiment_pdb)

    if experiment_pdb_file is None:
        print("no pdb file")
        continue

    if not os.path.exists(experiment_pdb_file):
        print(experiment_file)
        continue

    experiment_structure = parser.get_structure("experiment", experiment_pdb_file)
    generated_structure = parser.get_structure("generated", generated_file)

    experiment_atoms = [atom.get_coord() for atom in experiment_structure.get_atoms()]
    generated_atoms = [atom.get_coord() for atom in generated_structure.get_atoms()]

    n = min(len(experiment_atoms), len(generated_atoms))

    experiment_values = np.array(experiment_atoms[:n])
    generated_values = np.array(generated_atoms[:n])

    rmse = np.sqrt(np.mean(((generated_values - experiment_values) ** 2).sum(axis=1)))

    print(rmse)
    print("pdb done")

    results.append({
        "run number": file_id,
        "peptide": file_peptide,
        "hla": hla,
        "rmse": rmse
    })

df = pd.DataFrame(results)
df.to_csv("final_rmse.csv", index=False)
