TEpeaks
=========

Version: 0.1

TEpeaks takes ChIP-seq (and similar data) alignment files (BAM or BED),
identiifes narrow peaks, and is also able to do differential analysis over
peaks of two sets of libraries. It is an extension of MACS by adding the
funcionality of taking into account multi-reads, another normalization
method, bin correlation, and differential analysis. The differential
analysis is performed using DESeq. 


`MHammell Lab <http://hammelllab.labsites.cshl.edu/software>`_

Created by Ying Jin, Yuan Hao, Molly Hammell, June 2017

Contact: Ying Jin (yjin@cshl.edu)

Requirements
------------

gcc:       4.8 or greater 
boost:     1.53 or greater (please define variable BOOST_ROOT as the folder where boost header files were installed.)  


Installation
------------

1. Download compressed tarball.
2. Unpack tarball.
3. Navigate into unpacked directory.
4. Run the following:

    for linux like system:

    $ make -f Makefile.LINUX

    for MAC OS:

    $ make -f Makefile.MACOS  

5. copy TEpeaks to bin directory, for example

    $ cp TEpeaks ~/bin

*NOTE* In the above example, you must add

    ~/bin

to the PATH variable.

TEpeaks
=========

Usage
---------

::

    Usage: TEpeaks <CMD> [arguments] ..

    Where <CMD> can be one of:

        narrow        Call puntate peaks 
        broad         Call diffused peaks (not available yet)
        diff          Call differential peak analysis (will release soon)
        version       Prints version information

    Running TEpeaks <CMD> without arguments prints usage information for <CMD>

    usage: TEpeaks narrow -t treatment sample [treatment sample ] 
                        --tinput treatment input
                        -s genome  
                        [optional arguments]

    Required arguments:
      -t | --treatment [treatment sample 1 ]
      --tinput    treatment input 
      -s genome  (hg: human, mm: mouse, dm: d. melanogaster)

    Optional arguments:
        --format                  format of the files, BAM, BED, or BEDPE
        -a, --auto                    auto detect shiftsize for single-end reads (default: off)
        -f, --fraglen=INT             Fragment size 
        --ratio=DOUBLE            ratio between IP sample and input sample, IP/Input. By default, it's calculated from the data, but can also be set by user.  
        --keepDup=STRING          How to deal with duplicate reads. The valid values are 'auto', 'all', or 1 (default: auto)
        --shift=INT               Shift reads towards 3' end, if positive, or 5' end if negative. (default: 0)
        --lmfold=INT              lower limit of fold ratio against background to build
                              model(default: 10)
        --hmfold=INT              higher limit of fold ratio against background to build
                              model (default: 30)
        -n, --prjname=STRING          name of the prject (default: NONAME)
        --norm=STRING             normalization methods. sd (library size) or bc (bin
                              correlation) (default: sd)
        -p, --pval=DOUBLE             p-value cutoff (default: 1e-5)
        --fdr=DOUBLE              false discovery rate cutoff (default: 0.05)
        --toLarge                 Scale library size to large sample (default: off)
        -g, --gsize=INT               effective genome size (default: human genome 2.7e9)
                              (default: value is estimated from the input data)
        --threads=INT             Number of threads to use (default: 1)
        --pileup=INT              the minuim pileup required for peaks with multi-reads (default: 20)
        --fe=DOUBLE              the minuim fold enrichment required for peaks with multi-reads (default: 3)
        -i, --numItr=INT              Number of iterations (default: 0)


Example Command Lines
----------------------

::

    TEpeaks narrow -t IP.bam --tinput input.bam -s mm -o test --format BAM -i 20 --threads 4


*NOTE* BAM files must be either unsorted or sorted by queryname. 


Copying & distribution
======================


TEpeaks is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but *WITHOUT ANY WARRANTY*; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE*.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with TEToolKit.  If not, see `this website <http://www.gnu.org/licenses/>`_.


