from bambu.logo import logo
from .run import model_predict_wrapper
from rdkit import Chem
from argparse import ArgumentParser
import pandas as pd
import pickle
import os

def main():

    print(logo)

    argument_parser = ArgumentParser(description="Utilitary script to download sustance/compound data for a PubChem BioAssay")
    argument_parser.add_argument('--input', required=True, help="path of the input file containing molecules to be analyzed in .sdf, .mol2 or .smi/.smiles format")
    argument_parser.add_argument('--output', required=True, help="path of the output CSV file")
    argument_parser.add_argument('--preprocessor', required=True, help="path of the preprocessor produced by bambu-preprocess")
    argument_parser.add_argument('--model', required=True, help="path of the trained model produced by bambu-train")
    arguments = argument_parser.parse_args()

    predict(
        arguments.input, 
        arguments.output, 
        arguments.preprocessor, 
        arguments.model,
    )

def predict(input_file, output_file, preprocessor_file, model_file):

    with open(model_file, 'rb') as model_reader:
        model = pickle.loads(model_reader.read())
    with open(preprocessor_file, 'rb') as preprocessor_reader:
        preprocessor = pickle.loads(preprocessor_reader.read())

    _, input_fileext = os.path.splitext(input_file)

    if input_fileext == ".sdf":
        mols = Chem.SDMolSupplier(input_file)
    elif input_fileext == ".mol2":
        mols = [Chem.MolFromMolFile(input_file)]
    elif input_fileext in [".smiles", "smi"]:
        mols = Chem.rdmolfiles.SmilesMolSupplier(input_file)
    else:
        raise Exception(f"Unsuported extension '{input_fileext}' (use .sdf, .mol2 or .smi/.smiles files)")

    for m, mol in enumerate(mols):

        if m % 100 == 0 and m > 0:
            print(f"{m} molecules processed ...")

        predicted_activity, predicted_activity_proba = model_predict_wrapper(mol, model, preprocessor)

        df_output = pd.DataFrame(
            [{
                "query":m+1,
                "activity": predicted_activity,
                "activity_probability": predicted_activity_proba
            }], columns=["query", "activity", "activity_probability"]
        ).to_csv(
            output_file, 
            index=False,
            mode="w" if m == 0 else "a",
            header=True if m == 0 else False
        )

if __name__ == "__main__":
    main()