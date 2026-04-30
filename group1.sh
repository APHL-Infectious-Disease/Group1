#!/bin/bash

nextflow run . -profile singularity -resume --kraken2_db assets/kraken2db_v2/

rm -r -d ./work/*

rm -r -d ./results/fastqdl/*
rm -r -d ./results/pipeline_info/*