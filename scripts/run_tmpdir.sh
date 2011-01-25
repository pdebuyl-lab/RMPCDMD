#!/bin/bash

BASE=`pwd`

exefile=$1
configfile=$2

exename=`basename $exefile`

tmpdir=`mktemp -d`

cp $exefile $configfile $tmpdir/

cd $tmpdir

./$exename

$BASE/scripts/plot_energy.py

