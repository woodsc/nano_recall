# Nano Recall  

 Oxford Nanopore HIV caller.  Still in early development, so please get in touch.  There is likely a later version available.  This software is considered RESEARCH USE ONLY and that we make no warranties of its fitness for any purpose.

## Requirements

1)  You need to install the program Ruby (version 2.6 or greater).  You can check your version with:
    ruby --version

If running in windows, you'll need a version of ruby with the devkit installed,
you can find those at https://rubyinstaller.org/downloads/ .

2)  You need to run this program in the "Command Prompt" - the black box with the c:\prompt  

3) The first thing you need to do, (only once when you install the program, or again each time you update it
is to run the command:
    ruby setup.rb
This will install the alignment algorithm.

## Usage of the program

    Usage: nano_recall.rb [options]
        -h, --help                       Prints help
        -i, --input FILE                 Input file
            --batch-input FOLDER         Input folder
        -o, --output FILE                Output file
            --batch-output FOLDER        Output folder
            --use-sampleids              Output format SAMPLEID+BARCODE.ext
        -r  --regions "pr:1-99 rt:1-440 int:1-288"  Define aa regions for analysis (optional)

To process a single fastq (or fastq.gz)type :

    ruby nano_recall.rb -i example.fastq.gz -o output

To process a whole folder of fastq files:

    ruby nano_recall.rb --batch-input INPUTFOLDER --batch-output SAVEFOLDER --use-sampleids


So, for example, in Windows, if you installed the program into c:\nanorecall and your data lives in d:\data, and you want your data in a d:\output folder, you type

c:\nanorecall\ruby nano_recall.rb --batch-input d:/data/ --batch-output d:/output/

**NOTICE THAT THE "/" go the opposite way that you may expect for the data and output folders**

The --use-sampleids option only works in combination with --batch-output, and
will output as SAMPLEID+BARCODE.ext, UNLESS it finds a duplicate SAMPLEID+BARCODE
in another file, in which case it will output as FILENAME.ext.

If your sequencing doesn't reach the extent of each gene, you can narrow the
amino acid region to analyze with --regions "pr:1-99 rt:1-440 int:1-288".  A
gene can be removed entirely by setting the region to 0-0, such as "int:0-0".
If modified regions are used, a note will be added to the resistance report.

Settings are in the config/settings.txt file.  The only settings currently
interesting are:

    debug-genes=pr,rt,int
    optimization-coverage-limit-target=400
    alignment-optimize=true
    alignment-optimize-pad-size=6

The debug-genes setting lets you choose which genes you want to process.   For example, you can choose only pr,rt

The optimization-coverage-limit-target controls how much coverage is needed at
each base before it decides it has enough.   400 is a good choice to mimic Sanger data.

alignment-optimize turns on indel alignment optimizations.  At every insertion
or deletion in an alignment, it re-aligns a small region
(alignment-optimize-pad-size * 2) based on the most common indel-free sequence
in the sample.  This helps prevent many common misalignments, at the cost of
longer processing time.





### Alignments

Because of the large number of deletions, aligning one way or another can sometimes
make a big difference.  For example, fastq_SA_NGS059.fastq (13).gz  Amino 333
was interpreted as having G68% E29%.


    The common sequence around that region was:
    CAGGGAGATGAT

    However often a gap appeared in the alignment like:
    CAGGAGA-TGAT

    Which should have been aligned like:
    CAGGGAGATGAT
    CAGG-AGATGAT

This misalignment caused the false E call at a high percentage.

Also we are currently simply removing any insertions, which may also be causing
errors that we should be addressing.  We will need to investigate ways to make
our alignments as accurate as we can, or maybe even use the deletion information
to estimate a confidence value for each base?



### Speed optimizations

We could try to do some parallel processing to speed things up.  Ruby isn't
very good at this sort of thing by default, but there are improvements we
could do anyway.
