#!/bin/sh

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
mkdir -p $SCRIPT_DIR/../../data/exp_raw
curl https://datadryad.org/bitstream/handle/doi:10.5061/dryad.4539d/ExperimentalData.zip 
#unzip $SCRIPT_DIR/../../data/exp_raw/Wakamoto
