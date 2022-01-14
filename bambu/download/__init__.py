from bambu.logo import logo
from argparse import ArgumentParser
from rdkit import Chem
from tqdm import tqdm
from retrying import retry
import pandas as pd 
import requests
import json 

API_ENDPOINT = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'

def main():

    print(logo)
    
    argument_parser = ArgumentParser(prog="bambu-download", description="Downloads data from a PubChem BioAssays dataset")
    argument_parser.add_argument('--pubchem-assay-id', required=True, help="Assay ID")
    argument_parser.add_argument('--pubchem-InchI-chunksize', default=100, type=int, help="Number of InchI formulas to be downloaded per request")
    argument_parser.add_argument('--output', required=True, help="output file with the InchI formula and activity label for each compound/substance tested in the assay")

    arguments = argument_parser.parse_args()

    download_pubchem_assay_data(
        pubchem_assay_id=arguments.pubchem_assay_id, 
        pubchem_InchI_chunksize=arguments.pubchem_InchI_chunksize,
        output=arguments.output
    )

def download_pubchem_assay_data(pubchem_assay_id, output, pubchem_InchI_chunksize, pubchem_molecule_type='compounds'):

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
        request_uri = f'{API_ENDPOINT}/{pubchem_molecule_type.strip("s")}/{subset}/{",".join([str(pubchem_id) for pubchem_id in pubchem_molecules_ids_chunk])}/property/InChI/JSON'
        response_data = request_json_properties_data(request_uri)
        for pubchem_molecule in response_data:
            InChI_dict[pubchem_molecule[subset.upper()]] = pubchem_molecule['InChI']

    return InChI_dict

@retry(stop_max_attempt_number=100, wait_random_min=10000, wait_random_max=30000)
def request_json_properties_data(request_uri):
    response = requests.get(request_uri).text
    return json.loads(response)['PropertyTable']['Properties']

if __name__ == '__main__':
    main()