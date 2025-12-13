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

experiment_pdb = rcsb_search("RAKFKQLL", "B*08:01")

print(experiment_pdb)
