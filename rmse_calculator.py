import os
import requests
from Bio.PDB import PDBParser
import pandas as pd
import numpy as np

parser = PDBParser(QUIET=True)

experiment_directory = "ground_truths/"
generated_directory = "relaxed/"
os.makedirs(experiment_directory, exist_ok=True)

results = []

def rcsb_search(peptide, hla):
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "sequence",
                    "parameters": {
                        "evalue_cutoff": 0.001,
                        "identity_cutoff": 1,
                        "sequence_type": "protein",
                        "value": peptide
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "struct_keywords.pdbx_keywords",
                        "operator": "contains_phrase",
                        "value": hla
                    }
                }
            ]
        },
        "return_type": "entry"
    }

    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    response = requests.post(url, json=query)

    if response.status_code != 200:
        return []

    result = response.json().get("result_set", [])
    if not result:
        return None


    results_sorted = sorted(
        result,
        key=lambda x: x.get("score", 0),
        reverse=False
    )

    return results_sorted[0]["identifier"]

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

    hla = file_HLA.replace(" ", "")
    if hla.startswith("HLA-"):
        hla = hla
    else:
        hla = "HLA-" + hla

    experiment_pdb = rcsb_search(file_peptide, hla)

    if experiment_pdb is None:
        continue

    experiment_pdb_file = download_pdb(experiment_pdb)

    if experiment_pdb_file is None:
        continue

    experiment_file = os.path.join(experiment_directory, experiment_pdb_file)

    if not os.path.exists(experiment_file):
        continue

    experiment_structure = parser.get_structure("experiment", experiment_file)
    generated_structure = parser.get_structure("generated", generated_file)

    experiment_atoms = [atom.get_coord() for atom in experiment_structure.get_atoms()]
    generated_atoms = [atom.get_coord() for atom in generated_structure.get_atoms()]

    if experiment_file != "":
        n = min(len(experiment_atoms), len(generated_atoms))

        experiment_values = np.array(experiment_atoms[:n])
        generated_values = np.array(generated_atoms[:n])

        rmse = np.sqrt(np.mean(((generated_values - experiment_values) ** 2).sum(axis=1)))


        results.append({
            "id": file_id,
            "rmse": rmse
        })

df = pd.DataFrame(results)
df.to_csv("final_rmse.csv", index=False)
