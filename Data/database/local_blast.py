#!/usr/bin/python -tt
'''
Created on May 15, 2013
Updated on June 13, 2013

@author: Gokcen, Gabriel de la Cruz

@description: Holds a function that creates a standalone BLAST database using
BLAST command line makeblastdb
'''

# use for directory path and file manipulation
import os

# use to get system date and time
import datetime

# use to execute command line operations
import subprocess

# directory path to BLAST program
BLAST_PATH = os.getcwd() + '/../../Util/ncbi-blast-2.2.27+/bin/'


#*******************************************************************************
# @func_name: run_makeblastdb_Query()
# @desc: Builds a local blast db using the given inDBFile (in fasta format,
# nucleotide sequences) Database to be created gets the name outDBName, and can
# be referred just like that while querying. NOTE: Check out blastn command
# documentation to change output format.
# @return: None
#*******************************************************************************
def run_makeblastdb_Query(inDBFile, outDBName):
    # makeblastdb parameters
    # [-dbtype nucl] for nucleotide sequence
    # [-in <fasta_file>] input fasta file
    # [-out <db_name>] output local BLAST database name
    # [–parse_seqids] flag to enable retrieval of sequences based upon sequence identifiers
    subprocess.call([BLAST_PATH + 'makeblastdb', '-dbtype', 'nucl', '-in', inDBFile, '-out', outDBName, '-parse_seqids'])

##end run_makeblastdb_Query()
