#!/usr/bin/python -tt
'''
Created on June 14, 2013

@author: Gabriel de la Cruz

@description: Creates standalone BLAST Antibiotic Resistance Database (ARDB)
by querying data from 'ARDB' table in MySQL and converts into a FASTA file.
FASTA file is then converted into BLAST database by using BLAST command line
makeblastdb.
'''

# use for directory path and file manipulation
import os

# use for string manipulation
import string

# use to get system date and time
import datetime

# houses function run_makeblastdb_Query()
# use to create standalone BLAST database
import local_blast

# houses functions to connect to MySQL
# and run SQL statements
import mysqlc


#*******************************************************************************
# @func_name: create_fasta_file()
# @desc: Creates the FASTA file of all AR genes
# @return: FASTA filename of type <string>
#*******************************************************************************
def create_fasta_file(cursor, filename):
  
  print 'Creating FASTA file'

  # Start Timer
  startTime = datetime.datetime.now()
  print 'Start time: ', startTime.strftime('%H:%M:%S')

  # FASTA filename
  fasta_file = filename + '.fasta'
  
  # open file stream to write FASTA file
  fp = open(fasta_file, 'w')

  # traverse through each AR gene record in the SQL table 
  for (hit_id, resistance_type, phenotype, sequence) in cursor:
    # >gnl|ARDB|<hit_id> <resistance_type>|<phenotype>
    fp.write('>gnl|ARDB|%i %s|%s\n' % (hit_id, resistance_type, phenotype))
    # ATGCCCATAGCT...TAA
    fp.write('\n'.join(split_str_into_len(sequence)) + '\n')

  # closes file stream
  fp.close()

  # End Timer
  endTime = datetime.datetime.now()
  print 'End time: ', endTime.strftime('%H:%M:%S')
  print 'Time elapsed:', endTime - startTime

  # return string FASTA filename
  return fasta_file

##end create_fasta_file()


#*******************************************************************************
# @func_name: split_str_into_len()
# @desc: Split a string into chunks of specific length (default: 60)
# @return: List of string
# @reference: http://forums.devshed.com/python-programming-11/most-efficient-way
#             -of-splitting-a-string-into-fixed-size-390312.html
#*******************************************************************************
def split_str_into_len(sequence, length=60):
  
    return [sequence[i:i+length] for i in range(0, len(sequence), length)]

##end split_str_into_len()
  

#*******************************************************************************
# @func_name: create_blast_db()
# @desc: Creates the BLAST standalone database
# @return: None
#*******************************************************************************
def create_blast_db(fasta_file, blast_db_name):
  
  print 'Creating Antibiotic Resistance Database (ARDB)'

  # Start Timer
  startTime = datetime.datetime.now()
  print 'Start time: ', startTime.strftime('%H:%M:%S')

  local_blast.run_makeblastdb_Query(fasta_file, blast_db_name)

  # End Timer
  endTime = datetime.datetime.now()
  print 'End time: ', endTime.strftime('%H:%M:%S')
  print 'Time elapsed:', endTime - startTime
  
##end create_blast_db()


#*******************************************************************************
# @func_name: main()
# @desc: Program starts here
#*******************************************************************************
def main():

  # MySQL connection parameters
  db_username = 'racuser'     
  db_password = 'racuser123'
  db_host = 'localhost'
  db_name = 'RAC'

  # MySQL table name that holds all AR genes
  table_name = 'ARDB'

  # SQL statement to query all AR genes use to create FASTA file
  GET_AR_GENES = ('SELECT hit_id, resistance_type, phenotype, sequence '
                  'FROM %s ' % (table_name) +
                  'ORDER BY hit_id ASC')

  # current working directory
  cwd = os.getcwd()

  # connect to MySQL database
  db = mysqlc.connect_db(db_username, db_password, db_host, db_name)

  # query AR genes from specified table
  cursor = mysqlc.execute_query(db, True, GET_AR_GENES)

  # creates FASTA file
  fasta_filename = create_fasta_file(cursor, cwd + '/' + table_name)

  # closes cursor
  cursor.close()

  # creates BLAST standalone database
  create_blast_db(fasta_filename, cwd + '/' + table_name)

  # disconnect from MySQL database
  db.close()

##end main() 


if __name__ == '__main__':
  main()
