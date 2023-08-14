TEpeaks
=========

Version: 0.1

TEpeaks takes ChIP-seq (and similar data) alignment files and
identifies narrow ChIP peaks. It is an extension of MACS by adding the
functionality of handling reads with multiple alignments. It also allows
normalization by bin correlation. While it is currently optimized for 
narrow peak calling (e.g. TF binding sites), future updates will include
broad peak calling (e.g. repressive histone marks), and differential peak
analysis using DESeq2.

`MHammell Lab <http://hammelllab.labsites.cshl.edu/software>`_

Created by Ying Jin, Yuan Hao, Molly Gale Hammell, June 2017

Contact: mghcompbio@gmail.com

Requirements
------------

- gcc:       4.8 or greater
- boost:     1.53 or greater (please define variable BOOST_ROOT as the folder where boost header files were installed.)  


Installation
------------

1. Download compressed tarball.
2. Unpack tarball.
3. Navigate into unpacked directory.
4. Run the following:

      $ make -f Makefile.LINUX

5. copy TEpeaks to bin directory, for example:

      $ cp TEpeaks ~/bin

   *NOTE*: In the above example, you must add the binary directory
   to the PATH variable

      $ export PATH=~/bin:$PATH

TEpeaks
=========

Usage
---------

.. code::

    Usage: TEpeaks <CMD> [arguments] ..

    Where <CMD> can be one of:

        narrow        Call puntate peaks 
        version       Prints version information

    Running TEpeaks <CMD> without arguments prints usage information for <CMD>

    usage: TEpeaks narrow -t treatment sample [treatment sample ...]
                          -c control sample [control sample ...]
                          -o output directory
                          [optional arguments]

    Required arguments:
      -t | --treatment=STRING    IP sample(s) [BAM]
      -c | --control=STRING      Control (Input) sample(s) [BAM]
      -o | --outputdir=STRING    Directory to write output to

    Optional arguments:
      -f | --fraglen=INT         Fragment size (default: 200)
     --keepDup=STRING            How to deal with duplicate reads. The valid values are 'auto', 'all', or 1 (default: auto)
     --shift=INT                 Shift reads towards 3' end, if positive, or 5' end if negative. (default: 0)
     --lmfold=INT                Lower limit of fold ratio against background to build model (default: 10)
     --hmfold=INT                Higher limit of fold ratio against background to build model (default: 30)
     -n | --prjname=STRING       Project name used in output files (default: NONAME)
     -p | --pval=DOUBLE          P-value cutoff (default: 1e-5)
     --fdr=DOUBLE                False discovery rate cutoff (default: 0.05)
     --toLarge                   Scale library size to large sample (default: off)
     -s | --species=STRING       Species e.g., hs (Human hg19),  mm (Mouse mm9). (default: hs)
     -g | --gsize=INT            Effective genome size (default: human genome 2.7e9)
     --threads=INT               Number of threads to use (default: 1)
     --pileup=INT                The minuim pileup required for peaks with multi-reads (default: 20)
     --fe=DOUBLE                 The minuim fold enrichment required for peaks with multi-reads (default: 3)
     -i | --numItr=INT           Number of iterations (default: 50)


Example Command Lines
----------------------

.. code::

    TEpeaks narrow -t Pol2_IP_rep1.bam Pol2_IP_rep2.bam -c input_rep1.bam input_rep2.bam -o mouse_ChIP -s mm -i 20 -n mouse_Pol2


*NOTE*: BAM files must be either unsorted or sorted by queryname. 


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
along with TEpeaks.  If not, see `this website <http://www.gnu.org/licenses/>`_.


