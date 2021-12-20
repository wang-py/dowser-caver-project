#!/bin/bash

cp $1 files
cd files
./main.sh $1 $2
cp PredictedInternal.pdb ..


