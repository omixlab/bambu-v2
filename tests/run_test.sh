#!/usr/bin/env bash

PUBCHEM_BIOASSAY_ID=29
TEST_INPUT_FILE=pubchem_sample.sdf

if [ -f $TEST_INPUT_FILE.gz ]; then
	gzip -d -f $TEST_INPUT_FILE.gz
fi

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
	--max-iter 10000 \
	--time-budget 3600 \
	--estimators rf extra_tree

# runs the predictive model against a new set of molecules

bambu-predict \
	--input $TEST_INPUT_FILE \
	--preprocessor $PUBCHEM_BIOASSAY_ID\_preprocessor.pickle \
	--model $PUBCHEM_BIOASSAY_ID\_model.pickle \
	--output $PUBCHEM_BIOASSAY_ID\_predictions.csv
