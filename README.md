# Bambu

![](https://img.shields.io/docker/pulls/omixlab/bambu-qsar)

![](https://img.shields.io/pypi/dm/bambu-qsar)

![](media/logo.png)

Bambu (*BioAssays Model Builder*), is a simple tool to generate QSAR models based on [PubChem BioAssays](https://pubchem.ncbi.nlm.nih.gov/) datasets. It relies on [RDKit](https://rdkit.org/) and the [FLAML AutoML framework](https://github.com/microsoft/FLAML) and provides utilitaries for downloading and preprocessing datasets, as well training and running the predictive models.

## Try it!

[Try Bambu on this Google Colab Notebook ^^](https://colab.research.google.com/github/omixlab/bambu-v2/blob/main/notebooks/Bambu%20Google%20Colab%20Tutorial.ipynb)

## Installing

### Installing as a PyPI package using `pip`:

```
$ pip install bambu-qsar
```

**Note:** RDKit must be installed separately.

### Intalling as an environment using `conda` on Linux:

```
$ git clone https://github.com/omixlab/bambu-v2
$ cd bambu-qsar
$ make setup PLATFORM=linux
$ conda activate bambu-qsar
```

### Running it with Docker

```
$ docker run -ti omixlab/bambu-qsar:latest
```

## Downloading a PubChem BioAssays data

Downloads a PubChem BioAssays data and save in a CSV file, containing the InchI representation and the label indicating molecules that were found to be active or inactive against a given target.

```
$ bambu-download \
	--pubchem-assay-id 29 \
	--pubchem-InchI-chunksize 100 \
	--output 29_raw.csv
```

The generated output contains the columns `pubchem_molecule_id` (Substance ID or Compound ID, depending on the option selected during download), `InChI` and `activity`. Only the fields `InchI` and `activity` are used in futher steps. 

|pubchem_molecule_id|pubchem_molecule_type|InChI                                                                                                                                                                                                                                                                                                                                                                                                              |activity|
|-------------------|---------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------|
|596                |compounds            |InChI=1S/C9H13N3O5/c10-5-1-2-12(9(16)11-5)8-7(15)6(14)4(3-13)17-8/h1-2,4,6-8,13-15H,3H2,(H2,10,11,16)                                                                                                                                                                                                                                                                                                              |active  |
|1821               |compounds            |InChI=1S/C9H11FN2O6/c10-3-1-12(9(17)11-7(3)16)8-6(15)5(14)4(2-13)18-8/h1,4-6,8,13-15H,2H2,(H,11,16,17)                                                                                                                                                                                                                                                                                                             |active  |
|2019               |compounds            |InChI=1S/C62H86N12O16/c1-27(2)42-59(84)73-23-17-19-36(73)57(82)69(13)25-38(75)71(15)48(29(5)6)61(86)88-33(11)44(55(80)65-42)67-53(78)35-22-21-31(9)51-46(35)64-47-40(41(63)50(77)32(10)52(47)90-51)54(79)68-45-34(12)89-62(87)49(30(7)8)72(16)39(76)26-70(14)58(83)37-20-18-24-74(37)60(85)43(28(3)4)66-56(45)81/h21-22,27-30,33-34,36-37,42-45,48-49H,17-20,23-26,63H2,1-16H3,(H,65,80)(H,66,81)(H,67,78)(H,68,79)|active  |
|2082               |compounds            |InChI=1S/C12H15N3O2S/c1-3-6-18-8-4-5-9-10(7-8)14-11(13-9)15-12(16)17-2/h4-5,7H,3,6H2,1-2H3,(H2,13,14,15,16)                                                                                                                                                                                                                                                                                                        |active  |
|2569               |compounds            |InChI=1S/C15H19N3O5/c1-8-11(17-3-4-17)14(20)10(9(22-2)7-23-15(16)21)12(13(8)19)18-5-6-18/h9H,3-7H2,1-2H3,(H2,16,21)                                                                                                                                                                                                                                                                                                |active  |
|2674               |compounds            |InChI=1S/C29H26O10/c1-10(30)5-12-18-19-13(6-11(2)31)29(37-4)27(35)21-15(33)8-17-23(25(19)21)22-16(38-9-39-17)7-14(32)20(24(18)22)26(34)28(12)36-3/h7-8,10-11,30-31,34-35H,5-6,9H2,1-4H3                                                                                                                                                                                                                            |active  |
|2693               |compounds            |InChI=1S/C31H30N6O6S4/c1-33-25(42)30(15-38)34(2)23(40)28(33,44-46-30)12-17-13-36(21-11-7-4-8-18(17)21)27-14-29-24(41)35(3)31(16-39,47-45-29)26(43)37(29)22(27)32-20-10-6-5-9-19(20)27/h4-11,13,22,32,38-39H,12,14-16H2,1-3H3                                                                                                                                                                                       |active  |


## Computing descriptors or fingerprints

Computes molecule descriptors or Morgan fingerprints for a given datasets produced by `bambu-download` (or following the same format). The output also contains a `train` and `test` subsets, whose sizes are defined based on the `--train-test-split-percent` argument. The argument `--resample` might be used to perform a random undersampling in the dataset, as most HTS datasets are heavily umbalanced. The path passed to `--output` is used as template to generate the `train` and `test` file. In this case, `29_preprocess_train.csv` and `29_preprocess_test.csv` respectively.

```
$ bambu-preprocess \
	--input 29_raw.csv \
	--output 29_preprocess.csv \
	--output-preprocessor 29_preprocessor.pickle \
	--feature-type morgan-2048 \
	--train-test-split 0.75 \
	--undersample
``` 

## Train

Trains a classification model using the FLAML AutoML framework based on the `bambu-preprocess` output datasets. The user may adjust most of the `flaml.automl.AutoML` parameters using the command line arguments. In this case we are using an Extra Trees Classifier.

```
$ bambu-train \
	--input-train 29_preprocess_train.csv \
	--output 29_model.pickle \
	--time-budget 3600 \
	--estimators extra_tree
``` 

A list of all available estimators can be accessed using the command `bambu-train --list-estimators`. Currently, only `rf` (*Random Forest*) and `extra_tree` are available.

## Validation

An y-randomization validation can be performed using the command `bambu-validate`, which will compute accuracy, recall, precision, f1-score and ROC AUC training with the original training dataset and by validating 
with the test one, and furtherly randomizing the training labels and several times (`--randomizations`). For each randomization, classification metrics
are computed again and significancy value (p-value) is computed based
on the z-score-normalized metrics.

```
$ bambu-validate \
	--input-train 29_preprocess_train.csv \
	--input-test 29_preprocess_test.csv \
	--model 29_model.pickle \
	--output 29_model.validation.json \
	--randomizations 100
```

## Predict

Receives an inputs, preprocess it using a preprocess object (generated using `bambu-preprocess`) and then runs a classification model (generated using `bambu-train`). Results are saved in a CSV file.

```
$ bambu-predict \
	--input pubchem_compounds.sdf \
	--preprocessor 29_preprocessor.pickle \
	--model 29_model.pickle \
	--output 29_predictions.csv
``` 

# Contact

Feel free to open issues or pull requests! You may also contact us by email.

Dr. Frederico Schmitt Kremer, PhD. E-mail: [fred.s.kremer@gmail.com](fred.s.kremer@gmail.com).
