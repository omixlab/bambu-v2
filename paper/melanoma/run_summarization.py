from tqdm import tqdm
import pandas as pd
import json
import numpy
import glob

df = pd.DataFrame(columns=["assay", "model", "iterations", "recall", "precision", "accuracy", "f1-score"])

all_assays = [assay_id.strip("\n") for assay_id in open("assay_ids.txt")]
processed_assays = []

for json_file_path in tqdm(glob.glob("assays/*/*/*.json")):

    try:
        _, assay, model, json_file_basename = json_file_path.split("/")
        iterations = model.split("_")[-1]
        model      = "_".join(model.split("_")[:-1:])
        json_data  = json.loads(open(json_file_path).read())
        recall     = json_data["1.0"]["recall"]
        precision  = json_data["1.0"]["precision"]
        accuracy   = json_data["accuracy"]
        f1         = json_data["1.0"]["f1-score"]

        df_assay = pd.read_csv(f"assays/{assay}/{assay}_preprocess_test.csv")
        positive_examples = df_assay.query("activity == 1.0").shape[0]
        negative_examples = df_assay.query("activity == 0.0").shape[0]

    except:
        continue

    df = df.append({
        "assay": assay,
        "model": model,
        "iterations": iterations,
        "positive_examples": positive_examples,
        "negative_examples": negative_examples,
        "total_examples": positive_examples + negative_examples,
        "recall": recall,
        "precision": precision,
        "accuracy": accuracy,
        "f1-score": f1
    }, ignore_index=True)

    if assay not in processed_assays:
        processed_assays.append(assay)

df = df.query("total_examples >= 500")

print(df.groupby(["model", "iterations"]).mean())
print(df.sort_values(["accuracy"], ascending=False).head(10))

print(
    "Done: ", 
    len(processed_assays), 
    "/", 
    len(all_assays), 
    "(%.2f)"%((float(len(processed_assays))/float(len(all_assays)))*100)
)