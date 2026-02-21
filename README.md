News
----
Please check out our new Google newsgroup for TAPIS: https://groups.google.com/forum/?hl=en#!forum/tapis

***July 13, 2016***

**Version 1.2.1 released!**

GMAP version of 2015-07-23 is not compatible with TAPIS as a bug in their software does not provide strand information for all alignments.  We recommend using a more recent version.  Thank you to Lianfeng Gu from Haixia Institute of Science and Technology (HIST) for bringing this to our attention!

Summay
------
TAPIS (Transriptome Analysis Pipeline from Isoform Sequencing)
is a program for correcting and aligning long reads with/without
second generation reads, transcript clustering, novel and full-length
splice isoform detection, and identification and analysis of
polyadenylation (poly(A)) and alternative poly(A) (APA).

Install
-----------
 python setup.py build

 python setup.py install


**Note**:

To install to a local set of python libraries
use --home=/PATH/TO/LOCAL/LIBRARY

Required Packages
-------------------------
- SpliceGrapher
- bx-python
- pysam
- matplotlib
- NumPy

Running
---------
 Align reads to reference genome using *alignPacBio.py*

 Run analysis on aligned reads using *run_tapis.py*