#!/bin/bash

R --vanilla --slave --args --infile "regions.bed" --release "mm9" -v < ../main.R;
