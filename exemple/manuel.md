# Instructions for using « main.py » and « function.py »

Main.py (subsequently called “main”) is the pilot file from which functions defined in functions.py (subsequently called “functions”) are invoked.

DATA (containing the data source files) and fig (containing the output figures) folders needs to be created.

The current “main” file provides a list of examples drawing figures from the data files.

The scripts have been developed and used within Ubuntu 14.04, but have been run successfully on windows10 with python 2.7 and dependencies such as numpy and matplotlib as can for instance be provided by the free version of the Enthought canopy package.

Launch the script as follows:

- within the Enthought Canopy (64-bit) list of programs, launch the “Canopy 64-bit command prompt”
- in the command prompt, change directory to script location and run the command (adapt as necessary)

<pre>
python.exe main.py
</pre>



The data source files are exported from Geneious (tested up to Geneious R11) as indicated below.
- Reads (or portions of reads) are selected. Typically within the Mulikevirus project, these could correspond to reads containing host DNA after trimming of the phage DNA (left-ends or right-ends). Reads would have been mapped on the phage genome and parts of the protruding reads selected and exported. Exporting as contig will allow to keep the alignment and subsequently annotate conveniently all first flanking bases.
- One base (typically the first and last) will be annotated defining a new (specific) annotation type i.e. “start” or “stop” (it is convenient to orient the annotation to show which side was integrated)
- the contig is then saved as a list of sequences
- the annotated reads in the list are mapped on the target sequence (in the example, would be the host bacterial genome, Pseudomonas aeruginosa PAO1, PA14 etc. related to a given project)
- Option 1: from the resulting contig, annotations (limited to the relevant ones, here start and/or stop) are exported. Columns labelled “Name”, ”Type”, ”Minimum”, ”Direction”, ”Transferred From” are needed. The resulting csv file needs to keep only the transferred annotations. For trimming, edit file with notepad ++ or another text editor, not excel as excel may include additional characters etc when saving the csv.
- Option 2: from the resulting contig, transfer the annotations (the orange arrow “live annotate and predict” in “Geneious” (check the annotation type name as Geneious may remember the last name used which may not be in use). Save, select the option to not save to original reference sequence. A new sequence is created, containing the transferred annotations in addition to the source reference genome. The relevant annotations can then be exported from this file, with no subsequent need to trim the exported file. If you do not see a content in the “transferred from” column, check that the content does not appear lower when scrolling the list.
