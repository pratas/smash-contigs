# SMASH-CONTIGS #

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

By default, Smash has many parameters assigned in order to avoid the estimation, enabling only to set both reference and target files. However, these defaults are only estimated to detect rearrangements in primates. Therefore, for other purposes you might need to adjust context and threshold parameters. Moreover, for custom image maps you might also need to set other parameters, such as width and/or ratio (scale).
Only [refFile] and [TarFile] are mandatory.

#### Options meaning

| Parameters     | Meaning                                                                              |
|----------------|:-------------------------------------------------------------------------------------|
| -h             | It will print the parameters menu (help menu)                                        |
| -V             | It will print the smash version number, license type and authors.                    |
| -v             | It will print progress information such as positions of the patterns, times, etc.    |
| -f             | It will force to write over files already created.                                   |
| -c &#60;context&#62;   | Size of the FCM (Markov) context order (interval [1;28]). Contexts above 14 will be handled with a hash-table, where the implementation is approximately linear in memory relatively to the size of the sequence. When the sequence is very fragmented, or the species are somehow distant, or the sequencing/assembly process has low quality this value show not be very high. |
| -t &#60;threshold&#62; | It will be used to segment the high from the low regions of information content (interval [0;2]). For distant species this value might be slightly below 2 (such as 1.9). |
| -m &#60;mSize&#62;     | Minimum size of the block considered as a valid patters after each segmentation process. Values below 1 Million for primate chromosomes might emerge excessive valid patterns. However for other purposes, such as gene scale analysis, this value should be set almost to 1. |
| -i             | It will not show the information map regarding to inverted pattern regions. |
| -n             | It will not show the information map regarding to regular pattern regions (normal regions). |
| -r &#60;ratio&#62;     | Sets the ratio size of the image. Currently is fixed to 1000000 which is an estimated value to the medium of the primate chromosomes sizes relatively to the medium of the screen resolution. This parameter is not automatically adaptad since a fixed value will bring different size chromosomes to the same scale. Nevertheless, to use it in small sequences, namely bacterial genomes, this parameter might be adjusted. |
| -a &#60;alpha&#62;       | Probabilities estimator. This value relates a linear interpolation between the maximum likelihood estimator and the uniform distribution. This also shows that when the total number of events is large, the estimator behaves as a maximum likelihood estimator. Default value is set to 1000. |
| -s &#60;seed&#62;        | This is a parameter to the pseudo-random generation function. Different seed values ensure different generated values. |
| -w &#60;wSize&#62;       | The window size among with the sub-sampling is calculated automatically, nevertheless this value might be adjusted for special needs. |
| -wt &#60;wType&#62;      | Window filtering type. Types: 0, 1, 2 or 3. Type 0 stands for Hamming, 1 for Hann, 2 for Blackman, while 3 represents a rectangular window. |
| -d &#60;dSize&#62;       | Sub-sampling value. This value among with the window size is calculated automatically. Nevertheless, for special purposes this value might be adjusted. | 
| -wi &#60;width&#62;      | Thickness of the image for each sequence. Default value is set to 25. |
| -p &#60;posFile&#62;     | The positions from all the rearrangements detected on the run. |
| -o &#60;outFile&#62;     | The output SVG image filename. The default uses the concatenation of reference with the target filenames (adding the "svg" extension). Beware: if the files are not in the working directory this might have problems due to several types of characters (such as '/'). |
| [refFile]     | The reference filename. |
| [tarFile]     | The target filename. |

#### Positions file

A file will be created, reporting each rearrangement position, with the following characteristics: type, id, start position, end position, direction.

An example ca be seen below, where columns are separated by a tab (ascii code:9)

|type            |id     |start        | end            | direction |
|:---------------|:------|:------------|:---------------|:----------|
|TARGET          |1      |12890        |  9068115       | 0-regular |
|REFERENCE       |1      |5542700      |  14243450      | 0-regular |
|TARGET          |2      |9100340      |  13134910      | 0-regular |
|REFERENCE       |2      |14243450     |  16350965      | 0-regular |
|REFERENCE       |2      |16383190     |  18129785      | 0-regular |
|TARGET          |3      |13154245     |  18883850      | 0-regular |
|REFERENCE       |3      |18136230     |  19083645      | 0-regular |
|REFERENCE       |3      |19160985     |  23666040      | 0-regular |
|TARGET          |4      |18961190     |  20920470      | 1-inverted|
|REFERENCE       |4      |23840055     |  23975400      | 0-regular |
|REFERENCE       |4      |24001180     |  24639235      | 0-regular |
|REFERENCE       |4      |24697240     |  25754220      | 0-regular |

The first column reports if the regions are in the target or reference sequences and the second column sets an id for each similar region. The third and fourth columns, repectively, indicate the beginning and the end of each similar region, while the last column reports the direction (if was inverted or regular).

## CITATION ##

On using this software/method please cite:

Diogo Pratas, Raquel M. Silva, Armando J. Pinho, Paulo J. S. G. Ferreira. An alignment-free method to find and visualise rearrangements between pairs of DNA sequences. Sci. Rep. 5, 10203; doi: 10.1038/srep10203 (2015).

## ISSUES ##

For any issue let us know at [issues link](https://github.com/pratas/smash/issues).

## LICENSE ##

GPL v2.

For more information:
<pre>http://www.gnu.org/licenses/gpl-2.0.html</pre>

