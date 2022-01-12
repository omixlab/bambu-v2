from argparse import ArgumentParser
from flaml import AutoML
from sklearn.metrics import classification_report
import pandas as pd
import pickle

def main():
    
    argument_parser = ArgumentParser(description="")
    argument_parser.add_argument('--input-train', required=True)
    argument_parser.add_argument('--input-test', required=True)
    argument_parser.add_argument('--estimator-list', required=True)
    argument_parser.add_argument('--output', required=True)
    arguments = argument_parser.parse_args()

    train(
        arguments.input_train, 
        arguments.input_test, 
        arguments.output
    )

def train(input_train, input_test, output):

    df_train = pd.read_csv(input_train)
    df_test  = pd.read_csv(input_test)

    X_train = df_train.drop(['activity'], axis=1)
    y_train = df_train['activity']
    X_test  = df_test.drop(['activity'], axis=1)
    y_test  = df_test['activity']

    automl = AutoML()
    automl.fit(X_train, y_train, task="classification", estimator_list=["rf"])
    y_pred = automl.predict(X_test)
    report = classification_report(y_test, y_pred)

    print(report)
    
    with open(output, 'wb') as model_writer:
        model_writer.write(pickle.dumps(automl, protocol=pickle.HIGHEST_PROTOCOL))

    with open(output+'_classification_report.txt', 'w') as report_writer:
        report_writer.write(report)

if __name__ == "__main__":
    main()