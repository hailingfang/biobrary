# Biobrary

## Introduction

Biobrary is a python library, which contain data and methods for biological computation.

## Requirement  

* ete3

## Install  

* using pip  

`pip install biobrary --user`

* manally install  

`git clone https://github.com/benjaminfang/biobrary.git`  

`cd biobrary`  

`mv biobrary $PYTHONPATH`  

where PYTHONPATH is python library searching path.

## Usage  

```
import biobrary  
dir(biobrary)  
```

## Data and Method  

* bioparse  

    * FASTA  

        class for fasta file.

    * GTF  

        class for gtf file.


* tree

    * CircleNode  

        class for phylogenic tree traverse and operations. And divide tree to circle node according
        to phylogenic distance.
    
