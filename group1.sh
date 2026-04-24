#!/bin/bash

nextflow run . -profile singularity -resume

rm -r -d ./work/*

rm -r -d ./results/fastqdl/*
rm -r -d ./results/pipeline_info/*