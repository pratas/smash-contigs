# SMASH-CONTIGS #

![ScreenShot](/logo.png)

Smash-contigs is a fast tool to handle the detection of similarity between pairs of DNA sequences.
According to Smash, it is the lighter and faster approach to compute rearrangements at expense of some more memory.
This is a Linux server version. Current implementation is limited to a reference with maximum of 4 GB.

## INSTALLATION ##

<pre>
wget https://github.com/pratas/smash-contigs/archive/master.zip
unzip master.zip
cd smash-contigs-master
make
cp smash-contigs ..
cd ..
</pre>

## USAGE ##

To see the possible options type
<pre>
./smash-contigs
</pre>
or
<pre>
./smash-contigs -h
</pre>
These will print the following options:
<pre>
Usage: smash-contigs [OPTION]... [FILE] [FILE]                           
Similarity contigs mapping according to a reference sequence.            
                                                                         
Non-mandatory arguments:                                                 
                                                                         
  -h                         give this help,                             
  -V                         display version number,                     
  -v                         verbose mode (more information),            
  -k &#60k-mer&#62                 k-mer size [1;20],                          
  -m &#60minimum&#62               minimum similar block size,                 
  -r &#60repeats&#62               maximum repeats number,                     
  -e &#60editions&#62              number of allowed editions,                 
  -i                         do NOT map inversions,                      
  -n &#60nThreads&#62              number of threads,                          
  -o &#60FILE&#62                  output filename with positions,             
                                                                         
Mandatory arguments:                                                     
                                                                         
  &#60FILE&#62                     contigs filename (FASTA),                   
  &#60FILE&#62                     reference sequence filename (FASTA).        
                                                                         
Report bugs to &#60{pratas,ap,pjf}@ua.pt&#62.        
</pre>

By default, smash-contigs has many parameters assigned in order to avoid the estimation, enabling only to set both reference and target files. However, for specific purposes you might need to adjust several parameters. 

#### Options meaning

| Parameters     | Meaning                                                                              |
|----------------|:-------------------------------------------------------------------------------------|
| -h             | It will print the parameters menu (help menu).                                        |
| -V             | It will print the smash version number, license type and authors.                    |
| -v             | It will print progress information such as positions of the patterns, times, etc.    |
| -k &#60;k-mer&#62;   | Size of the matching word  (interval [1;20]). 
| -m &#60;minimum&#62;     | Minimum similarity size.
| -r &#60;repeats&#62;     | Maximum number of repeats (common value is between 50 and 500).
| -e &#60;editions&#62;     | Maxmim number of substitutions allowed in the repeat without being discarded.
| -i             | It will not detect inverted pattern regions. |
| -n             | Number of threads for computation. |
| -o &#60;outFile&#62;     | Output file with positions. The default name uses the match as prefix. |
| [refFile]     | The contigs filename. This should only be a FASTA format file. |
| [tarFile]     | The reference filename. This should only be a FASTA format file. |

#### Positions file

A file will be created, reporting each similarity position, with the following fields: 
CONTIG-ID, CONTIG-INIT, CONTIG-END, REFERENCE-ID, REFERENCE-INIT, REFERENCE-END. Reverse complements (inverted repeats) will be represented in the REFERENCE-INIT and REFERENCE-END, where REFERENCE-INIT is larget than REFERENCE-END.

## CITATION ##

On using this software/method please cite:

Diogo Pratas, Raquel M. Silva, Armando J. Pinho, Paulo J. S. G. Ferreira. An alignment-free method to find and visualise rearrangements between pairs of DNA sequences. Sci. Rep. 5, 10203; doi: 10.1038/srep10203 (2015).

## ISSUES ##

For any issue let us know at [issues link](https://github.com/pratas/smash-contigs/issues).

## LICENSE ##

GPL v2.

For more information:
<pre>http://www.gnu.org/licenses/gpl-2.0.html</pre>

