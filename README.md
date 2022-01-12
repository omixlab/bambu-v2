# Bambu

Bambu (*BioAssays Model Builder*), is a simple tool to generate QSAR models based on PubChem BioAssays datasets. It uses mainly on RDKit and the FLAML AutoML library and provides utilitaries for downloading, preprocessing, training and provisioning of predictive models in distributed environments.

## Setup 

```
$ conda env create --file environment.yml
```

## Downloading datasets

```
$ bambu-download
```

## Preprocessing

```
$ bambu-preprocess
```

## Training

```
$ bambu-train
```

## Using the model (CLI)

```
$ bambu-predict
```

## Using the model (Web, single worker)

```
$ bambu-server --host --port --basic-web-auth
```
