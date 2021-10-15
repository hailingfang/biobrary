# Biobrary

## Introduction
Biobrary is a python library, which contain data and methods for biological computation.

## Install  
* using pip  
`pip install biobrary --user`

* manally install  
`git clone https://github.com/benjaminfang/biobrary.git`  
`cd biobrary`  
`mv biobrary PYTHONPATH`  
where PYTHONPATH is python library searching path.

## Usage  
```
import biobrary  
dir(biobrary)  
```

## Data and Method  

* biopaser  
    * Fasta_parser  
        class for fasta file.
    * Gff_parser  
        class for gff file.


* biocondon  
    * CODON_AA  
        python dictionary of codon and amino acids.  
    * base_complement
    * start_codon
    * stop_codon



* amino_acids_mw  



* CircleNode  
    class for phylogenic tree traverse and operations. And divide tree to circle node according
    to phylogenic distance.
