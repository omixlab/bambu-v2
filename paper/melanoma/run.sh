#!/usr/bin/env bash

mkdir -p assays/

while read assay_id; do
	bash run_bambu_for_assay.sh $assay_id assays/$assay_id
done < assay_ids.txt
