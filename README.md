# R-looper sim

Rlooper with additional feature: ssDNA nicking and simulate peaks from a given DNA

## Synopsis

```
rlooper_sim --i <fastaFile.fa> --outdir <outdir>
```

## Installation

Prerequisites:
- Unix based system (Ubuntu, OSX, etc)
- Install git if required: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git. Use git to pull the repository. 
- Git clone the repository and `make`:

```
git clone git@github.com:chedinlab/rlooper_sim.git
cd rlooper_sim
make
```

## Examples

```
cd rlooper_sim

# Example, no nick
bin/rlooper_sim --i example/pFC9_small.fa --outdir example/pFC9_small_outdir

# Example, no nick, with junc energy -a 20 and sigma -0.07 (7%), getting peaks that are at least 50bp
bin/rlooper_sim --i example/pFC9_small.fa --outdir example/pFC9_small_outdir_a20_sigmamin7_minlen50 --a 20 --sigma -0.070 --minlength 50

# Example as above, but with nick at 300bp position and nicklen 1
bin/rlooper_sim --i example/pFC9_small.fa --outdir example/pFC9_small_outdir_a20_sigmamin7_minlen50_nick300_nicklen1 --a 20 --sigma -0.070 --minlength 50 --nick 300 --nicklen 1
```

## Usage

```
Usage: bin/rlooper_sim --i <fastaFile.fa> --outdir <outdir> --o <outprefix>

----------------
Options

--i/--fasta     [file.fa]  : fasta input.fa. IMPORTANT: HEADER HAS TO BE IN THIS FORMAT (Example from 1-3908):


>mydef range=mygene:1-3908 5'pad=0 3'pad=0 strand=- repeatMasking=none
ACGTCGT...
-o/--out        [string]   : output prefix name
--outdir        [string]   : output directory
--N             [int]      : specifies length of superhelical region in bp [1500]
--a [int]                  : junction energy in kcal/mol, for B DNA it's 10-11 [10]
--sigma         [float]    : superhelical density of the region [-0.07]
--minlength     [int]      : minimum rloop peak length in bp [50]
--nick          [int]      : specifies the ssDNA nick region [N/A]
--nicklen       [int]      : specifies the length of ssDNA nick region [N/A]
--beg           [int]      : get seq only from this bp from the input fasta [1]
--end           [int]      : get seq only until this bp from the input fasta [all]
--help/--h                 : print this message
--seed                     : Seed for RNG [42]

----------------
```

## Contributors

The original version of this software (rlooper) is developed by Robert Stolz (rstolzATucdavis.edu) at https://github.com/chedinlab/rlooper/
rlooper_sim was developed and maintained by Stella R Hartono (srhartonoATucdavis.edu)

Development of this software was funded in part by NIH grant GM120607 and NSF CAREER grant DMS1057284. 

## License
    rlooper_sim
    Copyright (C) 2025 Stella R Hartono and Robert Stolz

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    The license can be found at https://www.gnu.org/licenses/gpl-3.0.en.html. 
    By downloading or using this software, you agree to the terms and conditions of the license. 
