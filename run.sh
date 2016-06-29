#!/bin/bash

# Copyright 2016 Alberto Termanini


R --vanilla --slave --args --infile "regions.bed" --release "mm9" -v < main.R;
