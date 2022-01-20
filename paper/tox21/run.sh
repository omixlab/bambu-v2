#!/usr/bin/env bash

BAMBU_MAX_ITERATION=$1
BAMBU_TIME_BUDGET=$2

if [ -z "$BAMBU_MAX_ITERATION" ]; then
	BAMBU_MAX_ITERATION=100
fi

if [ -z "$BAMBU_TIME_BUDGET" ]; then
        BAMBU_TIME_BUDGET=3600
fi

while read assay_id; do
	if [ -d assays/$assay_id ]; then
		rm -rf assays/$assay_id
	fi
	bash run_bambu_for_assay.sh $assay_id assays/$assay_id $BAMBU_MAX_ITERATION $BAMBU_TIME_BUDGET
done < assay_ids.txt
