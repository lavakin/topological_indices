#!/usr/bin/env python3
import requests
from pysmiles import read_smiles
import calculate_indices as ind
import sys

def get_smiles(ids):
    try:
        response = requests.get('https://opsin.ch.cam.ac.uk/opsin/' + ids + '.smi')
        resp = response.text
        if len(resp) < 2000:
            return resp
        else:
            print("molecule not found")
            exit(1)
    except:
        print("molecule not found")
        exit(1)

G = read_smiles(get_smiles(sys.argv[1]), explicit_hydrogen=False)
print(ind.Connectivity.zagreb(G))
print(ind.Connectivity.general_sum_connectivity(G,1))
