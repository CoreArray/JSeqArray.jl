JSeqArray: data manipulation of whole-genome sequencing variants with SeqArray files in Julia
===

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
[![Build Status](https://travis-ci.org/CoreArray/JSeqArray.png)](https://travis-ci.org/CoreArray/JSeqArray)

pre-release version: v0.1.0


## Features

Data management of whole-genome sequence variant calls with thousands of individuals: genotypic data (e.g., SNVs, indels and structural variation calls) and annotations in SeqArray files are stored in an array-oriented and compressed manner, with efficient data access using the Julia programming language.

The SeqArray format is built on top of Genomic Data Structure (GDS) data format, and defines required data structure. GDS is a flexible and portable data container with hierarchical structure to store multiple scalable array-oriented data sets. It is suited for large-scale datasets, especially for data which are much larger than the available random-access memory. It also offers the efficient operations specifically designed for integers of less than 8 bits, since a diploid genotype usually occupies fewer bits than a byte. Data compression and decompression are available with relatively efficient random access.


## Installation

* Development version from Github, requiring `julia >= v0.5`
```julia
Pkg.status()
# install package dependencies
Pkg.clone("https://github.com/CoreArray/jugds.jl.git")
Pkg.build("jugds")

Pkg.clone("https://github.com/CoreArray/JSeqArray.jl.git")
Pkg.build("JSeqArray")
```


## Package Maintainer

Dr. Xiuwen Zheng ([zhengxwen@gmail.com](zhengxwen@gmail.com))


## Tutorials

* Learn X in Y minutes (where X=Julia): [http://learnxinyminutes.com/docs/julia/](http://learnxinyminutes.com/docs/julia/)


## Citation

#### Original paper (implemented in an [R/Bioconductor](http://bioconductor.org/packages/SeqArray) package):

[SeqArray](http://bioconductor.org/packages/SeqArray)

Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance data format for WGS variant calls. *Bioinformatics*. [DOI: 10.1093/bioinformatics/btx145](http://dx.doi.org/10.1093/bioinformatics/btx145).

#### Python package

[PySeqArray](https://github.com/CoreArray/PySeqArray)



## SeqArray File Download

* [1000 Genomes Project](http://bochet.gcc.biostat.washington.edu/seqarray/1000genomes)


## Examples

More examples will be posted here.
