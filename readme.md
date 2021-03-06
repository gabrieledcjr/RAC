Introduction
------------
The RAC project uses the *alignment-assembly algorithm* to detect single-nucleotide polymorphisms (SNPs) in antibiotic-resistant (AR) genes that were identified from a specific strain. 

In detail, this program filters out the unmapped sequences from a strain that was aligned through a sequencer (Illumina MiSeq). The program will have to identify what AR genes are present by finding similarities against a locally maintained Antibiotic Resistant Database (ARDB) using the Basic Local Alignment Search Tool (BLAST). After identifying the AR genes present in the strain, the high scoring sequences are then assembled under the particular AR gene sequence it belongs, using the CLC Genomics Workbench. Besides assembly, CLC provides a visual mapping of the assembled sequences as well as a summarized tabular format. Thus making it particularly easy to detect SNPs and interpret the results. 


Files
-----
**run.py**
- It runs through all the individual steps of the project

**bam_to_unmapped_fasta.py**
- Filters unmapped reads from BAM file (MiSeq output) and generates a sorted BAM file and a FASTA file

**blast_fasta.py**
- Identifies similarities between a FASTA file against the ARDB using Basic Local Alignment Search Tool (BLAST)

**parse_xml_specimen.py**
- Parses the results of the BLAST

**parse_xml_ARDB.py**
- Parses the results of the BLAST while finding similarities within the ARDB to create antibiotic resistant groups

**assemble_hit_files.py**
- Generate files to be assembled against specific AR genes identified by BLAST

**/Data/Database/create_database.py**
- Creates standalone BLAST Antibiotic Resistance Database (ARDB)


System Requirements
-------------------
All of the applications or software listed below are required, unless it is marked as optional, to be able to run the program. It is important to read the installation guide for each software as some of the softwares listed require other softwares to be installed in order for it to work properly. Lastly, select the correct version for your computer’s operating system.

**Python**
- Current release is Python 2.7.5. It is highly recommended to remain in Python 2.7 series as other python modules and softwares used in this program does not support the Python 3 series.
- This program was written entirely in Python.
- Download Python at http://www.python.org/download/releases/ to be able to launch the program.

**Biopython**
- Biopython is a free tool for computational biology written in Python that includes tools and Python libraries for projects in bioinformatics.
- For links to tutorial and more about Biopython, visit http://biopython.org/wiki/Main_Page.
- Download the installer at http://biopython.org/wiki/Download.

**SAMtools**
- Current release is SAMtools 0.1.19
- SAMtools is used to manipulate alignments in SAM/BAM file format.
- For links to tutorial and more about SAMtools, visit http://samtools.sourceforge.net.
- Download the installer at http://sourceforge.net/projects/samtools/files/ and save it under /RAC/Util folder.

**Pysam**
- Pysam is a Python module for reading and manipulating SAM/BAM files. It is basically a wrapper for the SAMtools program.
- For links to tutorial and more about Pysam, visit http://code.google.com/p/pysam/.
- Download the installer at http://code.google.com/p/pysam/downloads/list.

**BamUtil**
- Current release is BamUtil 1.0.9
- BamUtil is a repository of programs that manipulate SAM/BAM files
- For links to tutorial and more about BamUtil, visit http://genome.sph.umich.edu/wiki/BamUtil.
- Download the installer at http://genome.sph.umich.edu/wiki/BamUtil#Releases and save it under /RAC/Util folder.

**BLAST+** (Basic Local Alignment Search Tool)
- Current release is NCBI BLAST 2.2.27
- BLAST+ is a widely used application from the National Center for Biotechnology Information (NCBI) to find similarities in nucleotide or protein sequences. BLAST application is available online but for the purpose of this project, we are using it’s command line application.
- BLAST Command Line Application User Manual at http://www.ncbi.nlm.nih.gov/books/NBK1763/
- Nucleotide-Nucleotide BLAST difference between algorithms at http://www.ncbi.nlm.nih.gov/Class/MLACourse/Modules/BLAST/nucleotide_blast.html 
- Download the latest installer or source code at http://www.ncbi.nlm.nih.gov/books/NBK21097/ and save it under /RAC/Util folder.

**MySQL**
- Current release MySQL 5.5.31
- MySQL is an open source database from Oracle.
- Download the installer at http://dev.mysql.com/downloads/mysql/.

**MySQL Connector/Python**
- Current release Connector/Python 1.0.12
- MySQL Connector/Python is a driver for Python programs to connect to MySQL databases.
- For links to tutorial and more about Connector/Pyhon, visit http://dev.mysql.com/doc/connector-python/en/index.html.
- Download the installer at http://dev.mysql.com/downloads/connector/python/.

**phpMyAdmin** (optional)
- Current release phpMyAdmin 3.5.2.2
- phpMyAdmin is a free software used for MySQL administration. This software is marked as optional as the program can run even without this program.
- For more information and the installation link of phpMyAdmin, visit http://www.phpmyadmin.net/home_page/index.php.

**CLC Genomics Workbench**
- CLC Genomics Workbench is not necessarily required to run the program but it is required to complete the final step of the project which is to be able to know if SNPs exists from the identified antibiotic resistance genes. 
- This is a licensed software that needs to be purchased.
- For more information, visit http://www.clcbio.com/products/clc-genomics-workbench/.

Illustration of the program
---------------------------
![alt tag](https://raw.github.com/gabrieledcjr/RAC/master/program_flow.png)

Suggested Improvements
----------------------
- Although the program has been created with a set goal towards cross-compatibility between Macintosh and Windows computers, but the program has only been tested on a Mac computer with an OS X 10.8.4 (Mountain Lion) operating system. It is highly suggested to test the program on a Windows computer to ensure compatibility along with the external third-party modules and softwares required to run it.
- Due to PySAM’s limited means of parsing sequences, it is highly suggested to use Python’s native code for file manipulation. This is in relation to an algorithm in assemble_hit_files.py. An array of file positions (sort of indexing) was created when searching for sequences in the unmapped sorted BAM file to retrieve the original sequence for each high scoring hit. The problem with PySAM, you can only query values in two ways; first is parsing everything at once in a dictionary which would be memory intensive especially with huge files, and the second way was by using next() built-in function that points to the next sequence. But the problem with the built-in next function is if the sequence is in the current file position, you cannot retrieve the current value but instead it only lets you get the next sequence. So what ends up happening, you go to file position before it which is 50,000 sequences before the right one, and go through all 50,000 sequences until you get back to the current file position.
- ARDB table and AR Groups table should both remain on MySQL database since the goal is to make the program accessible in the Internet but future improvements should look into parsing specimen BLAST output and saving it on a SQLITE database. If the results of the BLAST are only important locally, then there is no need to pushed it on the MySQL database, instead SQLITE is a local database that does not need any server. Also, if the specimen tables are not needed, they can also be deleted.
- Create a Graphical User Interface to make it user-friendly especially for researchers who are not familiar with the command line.

