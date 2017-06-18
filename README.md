JSeqArray: data manipulation of whole-genome sequencing variants with SeqArray files in Julia
===

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)

[![Build Status](https://travis-ci.org/CoreArray/JSeqArray.jl.png)](https://travis-ci.org/CoreArray/JSeqArray.jl)

pre-release version: v0.1.0


## Features

Data management of whole-genome sequence variant calls with thousands of individuals: genotypic data (e.g., SNVs, indels and structural variation calls) and annotations in SeqArray files are stored in an array-oriented and compressed manner, with efficient data access using the Julia programming language.

The SeqArray format is built on top of Genomic Data Structure (GDS) data format, and defines required data structure. GDS is a flexible and portable data container with hierarchical structure to store multiple scalable array-oriented data sets. It is suited for large-scale datasets, especially for data which are much larger than the available random-access memory. It also offers the efficient operations specifically designed for integers of less than 8 bits, since a diploid genotype usually occupies fewer bits than a byte. Data compression and decompression are available with relatively efficient random access.


## Installation

* Development version from Github, requiring `julia >= v0.5`
* require [jugds](https://github.com/CoreArray/jugds.jl)

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

```julia
using JSeqArray

fn = joinpath(Pkg.dir(), "JSeqArray", "demo", "data", "1KG_phase1_release_v3_chr22.gds")
f = seqOpen(fn)
```
```
SeqArray File: JSeqArray/demo/data/1KG_phase1_release_v3_chr22.gds (1.1M)
+    [  ] *
|--+ description   [  ] *
|--+ sample.id   { Str8 1092 LZMA_ra(10.5%), 914B } *
|--+ variant.id   { Int32 19773 LZMA_ra(8.39%), 6.6K } *
|--+ position   { Int32 19773 LZMA_ra(52.0%), 41.1K } *
|--+ chromosome   { Str8 19773 LZMA_ra(0.28%), 166B } *
|--+ allele   { Str8 19773 LZMA_ra(22.7%), 111.9K } *
|--+ genotype   [  ] *
|  |--+ data   { Bit2 2x1092x19773 LZMA_ra(8.17%), 882.5K } *
|  |--+ extra.index   { Int32 3x0 LZMA_ra, 19B } *
|  \--+ extra   { Int16 0 LZMA_ra, 19B }
|--+ phase   [  ]
|  |--+ data   { Bit1 1092x19773 LZMA_ra(0.02%), 550B } *
|  |--+ extra.index   { Int32 3x0 LZMA_ra, 19B } *
|  \--+ extra   { Bit1 0 LZMA_ra, 19B }
|--+ annotation   [  ]
|  |--+ id   { Str8 19773 LZMA_ra(35.2%), 77.0K } *
|  |--+ qual   { Float32 19773 LZMA_ra(3.62%), 2.9K } *
|  |--+ filter   { Int32,factor 19773 LZMA_ra(0.21%), 170B } *
|  |--+ info   [  ]
|  \--+ format   [  ]
\--+ sample.annotation   [  ]
   |--+ Family.ID   { Str8 1092 LZMA_ra(15.3%), 1.1K }
   |--+ Population   { Str8 1092 LZMA_ra(5.08%), 222B }
   |--+ Gender   { Str8 1092 LZMA_ra(5.85%), 386B }
   \--+ Ancestry   { Str8 1092 LZMA_ra(2.43%), 233B }
```

```julia
# get genotype data (ploidy × sample × variant), 0xFF is missing value
seqGetData(f, "genotype")

seqClose(f)
```
```
2×1092×19773 Array{UInt8,3}:
[:, :, 1] =
 0x00  0x01  0x01  0x00  0x01  …  0x00  0x00  0x00  0x00  0x00
 0x00  0x00  0x00  0x01  0x00     0x01  0x01  0x01  0x01  0x01
[:, :, 2] =
 0x00  0x00  0x00  0x00  0x00  …  0x00  0x00  0x00  0x00  0x00
 0x00  0x01  0x01  0x00  0x01     0x00  0x00  0x00  0x00  0x00
[:, :, 3] =
 0x00  0x00  0x00  0x00  0x00  …  0x00  0x00  0x00  0x00  0x00
 0x00  0x00  0x00  0x00  0x00     0x00  0x00  0x00  0x00  0x00

...

[:, :, 19771] =
 0x00  0x00  0x00  0x00  0x00  …  0x00  0x00  0x00  0x00  0x00
 0x00  0x00  0x00  0x00  0x00     0x01  0x00  0x00  0x00  0x00
[:, :, 19772] =
 0x00  0x01  0x01  0x00  0x01  …  0x01  0x01  0x00  0x01  0x00
 0x01  0x00  0x00  0x00  0x00     0x01  0x01  0x00  0x00  0x00
[:, :, 19773] =
 0x00  0x00  0x00  0x00  0x00  …  0x00  0x00  0x00  0x00  0x01
 0x00  0x00  0x00  0x00  0x00     0x00  0x00  0x00  0x01  0x00
```

### More examples

Julia tutorial with SeqArray files: [demo/tutorial.ipynb](demo/tutorial.ipynb)

Julia tutorial with built-in parallel programming: [demo/tutorial_parallel.ipynb](demo/tutorial_parallel.ipynb)
