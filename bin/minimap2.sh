#!/bin/bash

FILE=${1}
REFERENCE=$2

minimap2 -ax map-ont ${REFERENCE} ${FILE} | samtools sort 
