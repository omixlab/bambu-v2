#!/usr/bin/env bash

PUBCHEM_BIOASSAY_ID=$1
PROJECT_DIRECTORY=$2

mkdir -p $2
cd $2

# download assay data

bambu-download \
	--pubchem-assay-id $PUBCHEM_BIOASSAY_ID \
	--pubchem-InchI-chunksize 100 \
	--output $PUBCHEM_BIOASSAY_ID\_raw.csv

# compute fingerprints

bambu-preprocess \
	--input $PUBCHEM_BIOASSAY_ID\_raw.csv \
	--output $PUBCHEM_BIOASSAY_ID\_preprocess.csv \
	--output-preprocessor $PUBCHEM_BIOASSAY_ID\_preprocessor.pickle \
	--feature-type morgan-1024 \
	--train-test-split 0.75 \
	--undersample

# train a predictive model

bambu-train \
	--input-train $PUBCHEM_BIOASSAY_ID\_preprocess_train.csv \
	--input-test $PUBCHEM_BIOASSAY_ID\_preprocess_test.csv \
	--output $PUBCHEM_BIOASSAY_ID\_model.pickle \
	--model-history \
	--max-iter 10 \
	--time-budget 10 \
	--estimators rf extra_tree 
