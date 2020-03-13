### Machine Learning Applications in Genetics
***
#### CS224 Final Project

Project entails haplotype phasing from genotypes containing ~40000 SNPs for 50 individuals.
Genotype data contains masked SNPs and therefore must be imputed.

Imputation was performed by using a max likelihood approach and assigning based on SNP transition frequencies.
Phasing was performed using a tiled variant of Clark's algorithm.

To impute masked data, run the following command.
Imputed data will be in text file *output*.
```bash
./phaser impute masked_data.txt
```

To see accuracy of imputation, divide the outputs of the following commands.
```bash
./phaser compare true_geno.txt output
./phaser total_mask masked_data.txt
```

Finally, to phase the genotype run the following. 
Phased data will be in text file *phased*.
```bash
./phaser phase output
```

Imputation and phasing can also be done all at once:
```bash
./phaser run masked_data.txt
```

To see accuracy of phasing run the Rscript
```bash
Rscript calculate_switch_accuracy.R true_phases.txt phased
```
