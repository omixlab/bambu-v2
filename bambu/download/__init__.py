from argparse import ArgumentParser
from rdkit import Chem
from tqdm import tqdm
import pandas as pd 
import requests
import json 

API_ENDPOINT = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'

def main():
    
    argument_parser = ArgumentParser(description="Utilitary script to download sustance/compound data for a PubChem BioAssay")
    argument_parser.add_argument('--pubchem-assay-id', required=True)
    argument_parser.add_argument('--pubchem-molecule-type', choices=["compounds","substances"], default="compounds")
    argument_parser.add_argument('--pubchem-InchI-chunksize', default=100, type=int)
    argument_parser.add_argument('--output', required=True)

    arguments = argument_parser.parse_args()

    download_pubchem_assay_data(
        pubchem_assay_id=arguments.pubchem_assay_id, 
        pubchem_molecule_type=arguments.pubchem_molecule_type,
        pubchem_InchI_chunksize=arguments.pubchem_InchI_chunksize,
        output=arguments.output
    )

def download_pubchem_assay_data(pubchem_assay_id, pubchem_molecule_type, output, pubchem_InchI_chunksize):

    df = pd.DataFrame(columns=['pubchem_molecule_id', 'pubchem_molecule_type', 'InChI', 'activity'])

    for activity in ['active', 'inactive']:

        pubchem_ids = get_assay_molecules_ids(
            pubchem_assay_id=pubchem_assay_id, 
            pubchem_molecule_type=pubchem_molecule_type,
            activity=activity
        )

        print(f'Downloading {activity} molecules ...')

        InChIs = get_molecules_InChIs(pubchem_ids, pubchem_molecule_type=pubchem_molecule_type, pubchem_InchI_chunksize=pubchem_InchI_chunksize)

        for pubchem_molecule_id, InChI in InChIs.items():

            df = df.append({
                'pubchem_molecule_id':pubchem_molecule_id, 
                'pubchem_molecule_type': pubchem_molecule_type,
                'InChI': InChI,
                'activity': activity
            }, ignore_index=True)

    df.to_csv(output, index=False)        

def get_assay_molecules_ids(pubchem_assay_id, pubchem_molecule_type='substances', activity="all"):
    
    if pubchem_molecule_type == "substances":
        subset = "sid"
    elif pubchem_molecule_type == "compounds":
        subset = "cid"
    else:
        raise Exception('molecule type not available (use "substances" or "compounds")')

    response = requests.get(f'{API_ENDPOINT}/assay/aid/{pubchem_assay_id}/{subset}s/JSON?{subset}s_type={activity}').text
    response_data = json.loads(response)
    return response_data['InformationList']['Information'][0][subset.upper()]

def get_molecules_InChIs(pubchem_molecules_ids, pubchem_molecule_type='substances', pubchem_InchI_chunksize=50):
    
    InChI_dict = {}
    
    if pubchem_molecule_type == "substances":
        subset = "sid"
    elif pubchem_molecule_type == "compounds":
        subset = "cid"
    else:
        raise Exception('molecule type not available (use "substances" or "compounds")')    

    for i in tqdm(range(0, len(pubchem_molecules_ids), pubchem_InchI_chunksize)):
        
        pubchem_molecules_ids_chunk = pubchem_molecules_ids[i:i+pubchem_InchI_chunksize]
        response = requests.get(f'{API_ENDPOINT}/{pubchem_molecule_type.strip("s")}/{subset}/{",".join([str(pubchem_id) for pubchem_id in pubchem_molecules_ids_chunk])}/property/InChI/JSON').text
        response_data = json.loads(response)
        
        for pubchem_molecule in response_data['PropertyTable']['Properties']:
            InChI_dict[pubchem_molecule[subset.upper()]] = pubchem_molecule['InChI']

    return InChI_dict

if __name__ == '__main__':
    main()