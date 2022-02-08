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

for bambu_max_iteration in 1 10 100 1000 10000 100000; do 
	for estimator in rf extra_tree decision_tree logistic_regression gradient_boosting; do
		mkdir $estimator\_$bambu_max_iteration
		cd $estimator\_$bambu_max_iteration
		bambu-train \
			--input-train ../$PUBCHEM_BIOASSAY_ID\_preprocess_train.csv \
			--input-test ../$PUBCHEM_BIOASSAY_ID\_preprocess_test.csv \
			--output $PUBCHEM_BIOASSAY_ID\_$estimator\_model.pickle \
			--model-history \
			--max-iter $bambu_max_iteration \
			--estimators $estimator \
			--threads 16 \
			--metric f1
		cd ..
	done
done
