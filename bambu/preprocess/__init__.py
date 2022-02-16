from bambu.logo import logo
from bambu.preprocess.preprocessors.descriptors import DescriptorsPreprocessor
from bambu.preprocess.preprocessors.morgan import MorganPreprocessor
from bambu.preprocess.preprocessors.mol2vec import Mol2VecPreprocessor
from rdkit import Chem
from pathlib import Path
from rdkit.Chem import AllChem
from argparse import Action, ArgumentParser, _HelpAction
from tqdm import tqdm
from sklearn.model_selection import train_test_split
from imblearn.under_sampling import RandomUnderSampler 
import argparse
import pandas as pd
import requests
import numpy as np
import pickle
import os
import sys

def main():

    if '--list-descriptors' in sys.argv:
        for feature in DescriptorsPreprocessor().features:
            print(feature)
        exit(0)

    if '--download-mol2vec-model' in sys.argv:
        user_home_directory = os.path.expanduser("~")
        download_file_directory = os.path.join(user_home_directory, "bambu", "mol2vec")
        os.makedirs(download_file_directory, exist_ok=True)
        download_file_path = os.path.join(download_file_directory, "mol2vec.model")

        response = requests.get("https://github.com/samoturk/mol2vec/blob/master/examples/models/model_300dim.pkl?raw=true")
        with open(download_file_path, 'wb') as handle:
            handle.write(response.content)
        exit(0)

    print(logo)

    argument_parser = ArgumentParser(prog="bambu-preprocess", description="Computes descriptors or fingerprints for molecules from a BioAssays")
    argument_parser.add_argument('--input', required=True, help="path of the input CSV file containing molecules InchI string and biological activity")
    argument_parser.add_argument('--output', required=True, help="path of the output CSV file containing the computed features")
    argument_parser.add_argument('--output-preprocessor', required=True, help="path of the output preprocessor object to be used by bambu-predict")
    argument_parser.add_argument('--feature-type', choices=['descriptors', 'mol2vec', 'morgan-32', 'morgan-64', 'morgan-128', 'morgan-256', 'morgan-512', 'morgan-1024', 'morgan-2048'], default='morgan-1028', help="type of feature to be computed")
    argument_parser.add_argument('--mol2vec-model-path', type=str, help="pre-trained mol2vec model", default=None)
    argument_parser.add_argument('--train-test-split-percent', type=float, default=0.75, help="percent of the dataset to be used for training")
    argument_parser.add_argument('--list-descriptors', default=False, action='store_true', help="list all descriptors available when using --feature-type descriptors")
    argument_parser.add_argument('--descriptors', nargs="+", default=None, help="list of descriptors to be used when using --feature-type descriptors. If not specified, all RDKit descriptors will be used")
    argument_parser.add_argument('--undersample', action='store_true', default=False, help="balance dataset using random undersampling before computing features")
    arguments = argument_parser.parse_args()

    if arguments.list_descriptors:
        for descriptor in DescriptorsPreprocessor().features:
            print(descriptor)
        exit(0)

    preprocess(
        arguments.input, 
        arguments.output, 
        arguments.output_preprocessor, 
        arguments.feature_type, 
        train_test_split_percent=arguments.train_test_split_percent, 
        undersample=arguments.undersample,
        descriptors=arguments.descriptors,
        mol2vec_model_path=arguments.mol2vec_model_path
    )

def preprocess(input_file, output_file, output_preprocessor_file, feature_type, train_test_split_percent=None, undersample=False, descriptors=None, mol2vec_model_path=''):

    df_input = pd.read_csv(input_file)
    
    X = df_input.drop(['activity'], axis=1)
    y = df_input['activity']

    if undersample:
        X, y = RandomUnderSampler().fit_resample(X, y)

    df_input = pd.DataFrame(X)
    df_input['activity'] = y

    if feature_type.startswith('morgan'):
        bits = int(feature_type.split('-')[1])
        preprocessor = MorganPreprocessor(bits=bits, radius=2)

    elif feature_type == "descriptors":
        preprocessor = DescriptorsPreprocessor()

    elif feature_type == "mol2vec":
        if mol2vec_model_path is None:
            raise Exception("The path to a pretrained model must be passed when using mol2vec features")
        if not os.path.isfile(mol2vec_model_path):
            raise Exception(f"The path '{mol2vec_model_path}' is not a valid mol2vec model")
        preprocessor = Mol2VecPreprocessor(pretrained_model=mol2vec_model_path)
    
    for r, row in tqdm(df_input.iterrows(), total=df_input.shape[0]):
    
        mol = Chem.inchi.MolFromInchi(row.InChI)

        try:
            mol_features = preprocessor.compute_features(mol)
            mol_features['activity'] = 1 if row.activity == "active" else 0
        except:
            continue

        pd.DataFrame([mol_features], columns=[*preprocessor.features, 'activity']).to_csv(
            output_file, 
            index=False, 
            mode='w' if r == 0 else 'a',
            header=True if r == 0 else False
        )

    with open(output_preprocessor_file, 'wb') as preprocessor_writer:
        preprocessor_writer.write(pickle.dumps(preprocessor, protocol=pickle.HIGHEST_PROTOCOL))
    
    if train_test_split_percent is not None:

        df_output = pd.read_csv(output_file)
        filepath, fileext = os.path.splitext(output_file)
        train_filepath = f'{filepath}_train{fileext}'
        test_filepath  = f'{filepath}_test{fileext}'

        X = df_output.drop(['activity'], axis=1)
        y = df_output['activity']

        X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=train_test_split_percent)

        X_train = X_train.astype(np.float32)
        X_test  = X_test.astype(np.float32)

        df_output_train = pd.DataFrame(X_train)
        df_output_train['activity'] = y_train
        df_output_train = clean_dataset(df_output_train)
        df_output_train.to_csv(train_filepath, index=False)

        df_output_test = pd.DataFrame(X_test)
        df_output_test['activity'] = y_test
        df_output_test = clean_dataset(df_output_test)
        df_output_test.to_csv(test_filepath, index=False)

def clean_dataset(df):
    df = df.astype(np.float32)
    df.dropna(inplace=True)
    indices_to_keep = ~df.isin([np.nan, np.inf, -np.inf]).any(1)
    return df[indices_to_keep].astype(np.float32)

if __name__ == "__main__":
    main()