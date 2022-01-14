# Bambu

Bambu (*BioAssays Model Builder*), is a simple tool to generate QSAR models based on PubChem BioAssays datasets. It uses mainly on RDKit and the FLAML AutoML library and provides utilitaries for downloading and preprocessing datasets, as well training and running the predictive models.

## Try it!

[Try Bambu on this Google Colab Notebook ^^](https://colab.research.google.com/github/omixlab/bambu-v2/blob/main/notebooks/Bambu%20Google%20Colab%20Tutorial.ipynb)

## Installing

### Installing as a conda package using `conda`:

*coming soon*

### Installing as a PyPI package using `pip`:

```
$ pip install bambu-qsar
```

**Note:** RDKit must be installed separately.

### Intalling as an environment using `conda`:

```
$ git clone 
$ cd bambu-qsar
$ conda env create --file environment.yml
$ conda activate bambu-qsar
```

## Downloading a PubChem BioAssays data

Downloads a PubChem BioAssays data and save in a CSV file, containing the InchI representation and the label indicating molecules that were found to be active or inactive against a given target.

```
$ bambu-download \
    --pubchem-assay-id 29 \
    --output 29_raw.csv
```

The generated output contains the columns `pubchem_molecule_id` (Substance ID or Compound ID, depending on the option selected during download), `InChI` and `activity`. Only the fields `InchI` and `activity` are used
in futher steps.

## Computing descriptors or fingerprints

Computes molecule descriptors or Morgan fingerprints for a given datasets produced by `bambu-download` (or following the same format). The output also contains a `train` and `test` subsets, whose sizes are defined based on the `--train-test-split-percent` argument. The argument `--resample` might be used to perform a random undersampling in the dataset, as most HTS datasets are heavily umbalanced.

```
$ bambu-preprocess \
    --input 29_raw.csv \
    --train-test-split-percent 0.75 \
    --feature-type descriptors \
    --undersample \
    --output 29_preprocessed.csv \
    --output-preprocessor 29_descriptor_preprocessor.pickle
``` 

## Train

Trains a classification model using the FLAML AutoML framework based on the `bambu-preprocess` output datasets. The user may adjust most of the `flaml.automl.AutoML` parameters using the command line arguments.
CLI arguments.

```
$ bambu-train \
	--input-train 29_preprocess_train.csv \
	--input-test 29_preprocess_test.csv \
	--output 29_model.pickle \
	--model-history \
	--max-iter 10 \
	--time-budget 10 \
	--estimators rf extra_tree
``` 

A list of all available estimators can be accessed using the command `bambu-train --list-estimators`. Currently, only `rf` (*Random Forest*) and `extra_tree` are available.

## Predict

Receives an inputs, preprocess it using a preprocess object (generated using `bambu-preprocess`) and then runs a classification model (generated using `bambu-train`). Results are saved in a CSV file.

```
$ bambu-predict \
	--input pubchem_compounds.sdf \
	--preprocessor 29_preprocessor.pickle \
	--model 29_model.pickle \
	--output 29_predictions.csv
``` 

