![alt tag](https://raw.github.com/jfortin1/TCGA_AB_Compartments/master/figures/try.png)

# HiC_AB_Compartments_Mouse

Analysis of Hi-C data has shown that the genome can be divided into two compartments
called A/B compartments. These compartments are cell-type specific and are
associated with open and closed chromatin. This GitHub repo contains the genome-wide A/B compartments computed for different mouse Hi-C datasets. All the data are at at resolution 40kb, for mm9. 

Citation for using these data:

Jean-Philippe Fortin, Kasper D. Hansen. *Reconstructing A/B compartments as revealed by Hi-C using long-range correlations in epigenetic data*. Genome Biology,16:180, 2015, doi:10.1186/s13059-015-0741-y

and 

Dixon JR, Selvaraj S, Yue F, Kim A, Li Y, Shen Y, et al. Topological domains in mammalian genomes identified by analysis of chromatin interactions. Nature. 2012;485:376–80. doi:10.1038/nature11082.


## Data format

The data are saved as tab-delimited text files. The first three columns represent the genomic coordinates of the bin, the fourth column reports the value of the eigenvector for the bin, and the fifth column categorizes the bin as "open" or "closed" (A and B compartment respectively). 
```{}
chr    start    end    eigen    domain
chr1	500000	599999	-0.448558697994734	open
chr1	600000	699999	-0.317063576958963	open
chr1	700000	799999	-0.0518659472420911	open
chr1	800000	899999	0.0711965336724111	closed
chr1	900000	999999	0.098443079780858	closed
chr1	1000000	1099999	-0.051011778069161	open
chr1	1100000	1199999	-0.318273718360445	open
chr1	1200000	1299999	-0.53951736779136	open
chr1	1300000	1399999	-0.585600098287715	open
```

## Description of the files

The files are in the directory `data`

| Files        | Description  | Reference Dataset
| ------------- |-------------| ------------- |
| hic_compartments_40kb_mESC.txt     | Calculated A/B compartments for Hi-C mESC | [1] | 
| hic_compartments_40kb_mCortex.txt | Calculated A/B compartments for Hi-C mCortex | [1] |

The methodology to compute the A/B compartments is described in 

Jean-Philippe Fortin, Kasper D. Hansen. *Reconstructing A/B compartments as revealed by Hi-C using long-range correlations in epigenetic data*. Genome Biology,16:180, 2015, doi:10.1186/s13059-015-0741-y

## Visualization 

The script [visualization.R](https://github.com/Jfortin1/TCGA_AB_Compartments/blob/master/scripts/visualization.R) is a script to load the data into R and create a barplot of the eigenvector values (see image above), as well as a little function called `convert2GRanges()` that will convert the data to a GRanges object for convenience. Positive and negative values represent closed and open domains respectively. 

## Note about preprocessing

To decide of the sign of the eigenvectors, we used GC content averaged within each bin. We chose the sign of the eigenvectors so that positive values of the eigenvectors (closed domain) correlate negatively with GC content. 

## References

[1] Dixon JR, Selvaraj S, Yue F, Kim A, Li Y, Shen Y, et al. Topological domains in mammalian genomes identified by analysis of chromatin interactions. Nature. 2012;485:376–80. doi:10.1038/nature11082.




