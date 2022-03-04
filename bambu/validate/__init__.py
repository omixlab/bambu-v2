from bambu.logo import logo
from argparse import ArgumentParser
from y_scramble import Scrambler
from sklearn.metrics import SCORERS
import numpy as np
import scipy
import pandas as pd
import json
import pickle
import copy
import tqdm

def main():

    print(logo)

    argument_parser = ArgumentParser(prog="bambu-validate", description="validates a classification model generated by bambu-train using y-randomization")
    argument_parser.add_argument('--input-train', required=True, help="path to CSV file containing training set generated by bambu-preprocess")
    argument_parser.add_argument('--input-test', required=True, help="path to CSV file containing test set generated by bambu-preprocess")
    argument_parser.add_argument('--randomizations', type=int, help="number of randomizations to rum", default=100)
    argument_parser.add_argument('--model', required=True, help="path of the trained model produced by bambu-train")
    argument_parser.add_argument('--output', required=True, help="path to the output JSON file containing the validation results")
    arguments = argument_parser.parse_args()

    validate(input_train=arguments.input_train, input_test=arguments.input_test, model_path=arguments.model, output=arguments.output, randomizations=arguments.randomizations)

def validate(input_train, input_test, model_path, output, randomizations=100):

    df_train = pd.read_csv(input_train)
    df_test  = pd.read_csv(input_test)

    X_train = df_train.drop(['activity'], axis=1)
    y_train = df_train['activity']
    X_test  = df_test.drop(['activity'], axis=1)
    y_test  = df_test['activity'] 

    with open(model_path, 'rb') as model_reader:
        automl = pickle.load(model_reader)
        model  = automl._trained_estimator.estimator

    report = validate_model(model, X_train, y_train, X_test, y_test, randomizations=randomizations)

    with open(output, 'w') as report_writer:
        report_writer.write(json.dumps(report))

def validate_model(model, X_train, y_train, X_test, y_test, metrics=["accuracy", "recall", "precision", "f1", "roc_auc"], randomizations=100):
    
    model_validation_data = {"raw_scores": {}}
    model = model.fit(X_train, y_train)

    for metric in metrics:
        
        scorer = SCORERS.get(metric)
        model_validation_data["raw_scores"][metric]  = [scorer(model, X_test, y_test)]

    for _ in tqdm.tqdm(range(randomizations)):
       
        y_train_scrambled = copy.copy(y_train)
        np.random.shuffle(y_train_scrambled)
        model.fit(X_train, y_train_scrambled)

        for metric in metrics:
            scorer = SCORERS.get(metric)
            model_validation_data["raw_scores"][metric].append(scorer(model, X_test, y_test))

    model_validation_data["zscores"] = {}
    model_validation_data["pvalues"] = {}

    for metric in metrics:
        
        model_validation_data["zscores"][metric] = scipy.stats.zscore(model_validation_data["raw_scores"][metric])
        model_validation_data["pvalues"][metric] = scipy.stats.norm.sf(abs(model_validation_data["zscores"][metric]))*2
        model_validation_data["zscores"][metric] = list(model_validation_data["zscores"][metric])
        model_validation_data["pvalues"][metric] = list(model_validation_data["pvalues"][metric])

    return model_validation_data
