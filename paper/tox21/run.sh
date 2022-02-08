#!/usr/bin/env bash

while read assay_id; do
	echo $assay_id
	bash run_bambu_for_assay.sh $assay_id assays/$assay_id #&> /dev/null
done < assay_ids.txt