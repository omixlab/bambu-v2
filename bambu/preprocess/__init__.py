from bambu.preprocess.preprocessors.descriptors import DescriptorsPreprocessor
from bambu.preprocess.preprocessors.morgan import MorganPreprocessor
from rdkit import Chem
from rdkit.Chem import AllChem
from argparse import ArgumentParser
from tqdm import tqdm
from sklearn.model_selection import train_test_split
from imblearn.under_sampling import RandomUnderSampler 
import pandas as pd
import numpy as np
import pickle
import os

def main():
    
    argument_parser = ArgumentParser(description="Utilitary script to download sustance/compound data for a PubChem BioAssay")
    argument_parser.add_argument('--input', required=True)
    argument_parser.add_argument('--output', required=True)
    argument_parser.add_argument('--output-preprocessor', required=True)
    argument_parser.add_argument('--feature-type', choices=['descriptors', 'morgan-1024', 'morgan-2048'], default='morgan-1028')
    argument_parser.add_argument('--train-test-split-percent', type=float, default=0.75)
    argument_parser.add_argument('--undersample', action='store_true', default=False)
    arguments = argument_parser.parse_args()

    preprocess(
        arguments.input, 
        arguments.output, 
        arguments.output_preprocessor, 
        arguments.feature_type, 
        train_test_split_percent=arguments.train_test_split_percent, 
        undersample=arguments.undersample
    )

def preprocess(input_file, output_file, output_preprocessor_file, feature_type, train_test_split_percent=None, undersample=False):

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

        df_output_train = pd.DataFrame(X_train)
        df_output_train['activity'] = y_train
        df_output_train.to_csv(train_filepath, index=False)

        df_output_test = pd.DataFrame(X_test)
        df_output_test['activity'] = y_test
        df_output_test.to_csv(test_filepath, index=False)


if __name__ == "__main__":
    main()