The mRNA Markup workflow represents comprehensive annotation and primary analysis of a set of transcripts. Such sets are currently abundantly generated from assembly of EST sequences, and with increasing read lengths of novel sequencing technologies, there will be similar assemblies of RNA-Seq data from many species and sampling conditions. Using NCBI BLAST+, MuSeqBox (Multi-query sequence BLAST output examination with MuSeqBox developed by the Brendel Group), and Linux shell scripting, an input transcript set is partitioned in many ways, including into contaminants (sequencing artifacts), potential chimeras, likely full-length protein-coding mRNAs, miRNAs, and potentially novel transcripts for further analysis with other programs. (more informations about mRNA Markup workflow, you can find on the this site: http://bioservices.usd.edu/mrnamarkup.html)

Before you start makes sure that ncbi-blast+ package (sudo apt-get install ncbi-blast+) and MuSeqBox (http://www.plantgdb.org/MuSeqBox/download.php) is already installed in your computer.

First, in folder where you have projekt1.py, you have to make folders:

	1) bazy (there are all databases: 
		Vector Database: UniVec (ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec)
		Bacteria Database: E. Coli from NCBI Nucleotide (http://www.bioextract.org/ViewFile?filename=/usr/local/BioStreamServer/tmpFiles/mRNA_Markup_Start/1291915878836/BacteriaDB)
		Reference Database: ATpepTAIR10 (A set of Arabidopsis protein sequences available at http://www.plantgdb.org/XGDB/phplib/download.php?GDB=At)
		All Protein Database: UniRef90-Viridiplantae (http://www.uniprot.org/uniref/?query=identity:0.9+taxonomy:33090&format=*&compress=yes)
		Protein Domain Database: NCBI’s CDD (ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz)

	2) input (your input, you should name your input file 'Input')

	3) output (your all outputs (without summary, because it will be created in folder where is projekt1.py))

Also, folder with MuSeqBox should be here.

Now, you can run program ./projekt1.py. 
