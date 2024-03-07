#!/bin/sh

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
mkdir -p $SCRIPT_DIR/../../data/exp_raw
git clone https://gitlab.com/MEKlab/growth-sos-msb-2021.git $SCRIPT_DIR/../../data/exp_raw/growth-sos-msb-2021

