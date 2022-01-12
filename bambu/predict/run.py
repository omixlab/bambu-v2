import pandas as pd
import numpy as np

def model_predict_wrapper(mol, model, preprocessor):
    mol_features_computed = True
    
    try:
        mol_features = preprocessor.compute_features(mol)
    except:
        mol_features_computed = False
    
    if mol_features_computed:
        df_features  = pd.DataFrame([mol_features], columns=preprocessor.features)
        predicted_activity = model.predict(df_features)[0]
        predicted_activity_proba = model.predict_proba(df_features)

        if predicted_activity_proba.shape[1] == 2:
            predicted_activity_proba = predicted_activity_proba[0][1]
        else:
            predicted_activity_proba = predicted_activity_proba[0]
    else:
        predicted_activity = None
        predicted_activity_proba = None

    return predicted_activity, predicted_activity_proba
