# Bambu

Bambu (*BioAssays Model Builder*), is a simple tool to generate QSAR models based on PubChem BioAssays datasets. It uses mainly on RDKit and the FLAML AutoML library and provides utilitaries for downloading, preprocessing, training and provisioning of predictive models in distributed environments.


## Installing

### Installing as a PyPI package using `pip`:

```
$ pip install bambu-qsar
```

* Note: * RDKit must be installed separately.

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
    --pubchem-molecule-type compounds \
    --output 29_raw.csv
```

The generated output contains the columns `pubchem_molecule_id` (Substance ID or Compound ID, depending on the option selected during download), `pubchem_molecule_type` (`compounds` or `substances`), `InChI` and `activity`. Only the fields `InchI` and `activity` are used
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

Trains a classification model using the FLAML AutoML framework based on the `bambu-preprocess` output datasets.

```
$ bambu-train \
    --input-train \
    --input-test \
    --output 29_model.picle \
``` 



## Predict

```
$ bambu-predict
``` 

