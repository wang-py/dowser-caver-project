#!/bin/bash

DOWSER_O_INPUT=$1

NUM_OF_WATER=$(awk 'END {print NR}' $DOWSER_O_INPUT)

echo $NUM_OF_WATER