#!/usr/bin/env bash

PUBCHEM_BIOASSAY_ID=$1
PROJECT_DIRECTORY=$2
BAMBU_MAX_ITERATION=$3
BAMBU_TIME_BUDGET=$4
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

for estimator in rf extra_tree decision_tree svm logistic_regression gradient_boosting; do
	mkdir $estimator
	cd $estimator
	bambu-train \
		--input-train ../$PUBCHEM_BIOASSAY_ID\_preprocess_train.csv \
		--input-test ../$PUBCHEM_BIOASSAY_ID\_preprocess_test.csv \
		--output $PUBCHEM_BIOASSAY_ID\_$estimator\_model.pickle \
		--model-history \
		--max-iter $BAMBU_MAX_ITERATION \
		--time-budget $BAMBU_TIME_BUDGET \
		--estimators $estimator
	cd ..
done