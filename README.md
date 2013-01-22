

**mRNA_Markup**


1. Description
2. Installation requirements
3. Formatting a local nucleotide and peptide database using BLAST+
4. Integration with a local version of GALAXY



<ol>
<li>Bird</li>
<li>McHale</li>
<li>Parish</li>
</ol>

1. Description

This Galaxy tool is based on BioExtract mRNA Markup workflow.
This workflow represents comprehensive annotation 
and primary analysis of a set of transcripts.

It consists of several "steps" each using several analytic tools:

0. Submit input: initial mRNA (query file) 
		 bacterial contamination database in FASTA format (a file representing typical bacterial hosts) 
		 reference protein database in FASTA format (proteins most likely to have homologs in the mRNA translations of the input)
		 comprehensive protein database in FASTA format (to be searched when the reference protein set did not give hits)
1. Eliminate Vector Contamination
2. Eliminate bacterial contamination
3. Find matches in a reference protein database 
3.1. Identify potential full-length coding sequences 
3.2. Identify potential chimeric sequences
4. Find matches in a Comprehensive Protein database
5. Find matches in Protein Domain Database
6. Produce summary report


For more information see:
(1) http://www.bioextract.org/scripts/index.jsp
(2) http://bioservices.usd.edu/mrnamarkup.html




2. Installation requirements

Before you start makes sure that 
(1) ncbi-blast+
(2) biopython
(3) MuSeqBox
are already installed on your computer.

Note: Version 4.1 of MuSeqBox (available from 
http://www.plantgdb.org/MuSeqBox/MuSeqBox.php)
is based on BLAST+ applications and was tested with ncbi-blast-2.2.24, 
available from http://blast.ncbi.nlm.nih.gov/.  
As there are numerous changes in the output
format between BLAST+ and the legacy BLAST applications, MuSeqBox will fail
to parse output from older BLAST versions.




3. Formatting a local nucleotide and peptide database using BLAST+

The workflow uses a sample input file consisting of Arabidopsis mRNA and searches the following BLAST databases by default:

BacteriaDB and RefProtDB	 	http://www.bioextract.org/Download?action=zip&file=AllOutput.zip&src=/usr/local/BioStreamServer/tmpFiles/mRNA_Markup_Start/1291915878836
AllProtDB 				http://www.uniprot.org/uniref/?query=identity:0.9+taxonomy:33090&format=*&compress=yes
Smart_LE/Smart				ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/big_endian/
UniVec					ftp://ftp.ncbi.nih.gov/pub/UniVec/


Command for formatting nucleotide databases:
makeblastdb -dbtype nucl -in DATABASE -parse_seqids


Command for formatting for protein databases:
makeblastdb -dbtype prot -in DATABASE -parse_seqids


Note: Smart databases do not need formatting. This databaes is already formatted and packed and available from the website mentioned above.



4. Integration with a local version of GALAXY

Simply copy the galaxy_dist directory from the repository onto the GALAXY directory tree.
The mRNA_Markup was tested only with a local GALAXY version.


********************
Katarzyna Wręczycka 
email: kw292555@students.mimuw.edu.pl
18.01.2013

