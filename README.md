# doGREAT

This R script performs a GREAT analysis starting from a set of genomic regions (BED file).
You can use this script to automate the analysis one or more set of genomic regions without manually download the results from the GREAT web site.
The script saves one xls file for each available ontology ("GO Molecular Function", "InterPro" etc.) in a given output directoy.
Results can be filtered by FDR and binomial fold enrichment.

For details about GREAT click here: http://bejerano.stanford.edu/great/public/html/

## Installation

Just dowload and run! 

Required packages:
- getopt (CRAN)
- rGREAT (Bioconductor)


## Usage

Run `R --vanilla --slave --args --help < main.R` to have the complete list of parameters.

To run a GREAT analysis with genome release "mm9" and to filter results by FDR <= 0.05 and binomial fold enrichment >= 2,
run `R --vanilla --slave --args --infile "regions.bed" --release "mm9" -v < main.R`.
See run.sh in example folder.

## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## History

This is the first release!

## Credits

Alberto Termanini

## License

See LICENCE file