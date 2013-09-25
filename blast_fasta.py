#!/usr/bin/python -tt
'''
Created on May 24, 2013

@author: Gabriel de la Cruz

@description: Divides a huge FASTA file into smaller temporary FASTA files and each
temporary FASTA file is BLASTed against the standalone Antibioitic Resistance
Database (ARDB). The expected BLAST result are multiple XML files.
'''

# use for command line arguments
import sys

# use for directory path and file manipulation
import os

# use to get system date and time
import datetime

# [BioPython] third party modules use to manipulate sequence files
from Bio import SeqIO
from Bio import Seq
from Bio.Blast.Applications import NcbiblastnCommandline


# File path for BLAST software
BLAST_PATH = os.getcwd() + '/Util/ncbi-blast-2.2.27+/bin/'

# Report log List
REPORT_LOG = []


#*******************************************************************************
# @func_name: batch_iterator()
# @desc: This can be used on any iterator, for example to batch up SeqRecord
#        objects from Bio.SeqIO.parse(...), or to batch Alignment objects from
#        Bio.AlignIO.parse(...), or simply lines from a file handle.
#
#        This is a generator function, and it returns lists of the entries from
#        the supplied iterator.  Each list will have batch_size entries,
#        although the final list may be shorter.
# @return: Returns lists of length batch_size
# @reference: http://biopython.org/wiki/Split_large_file
#*******************************************************************************
def batch_iterator(iterator, batch_size) :

  entry = True #Make sure we loop once
  
  while entry:
    batch = []
    
    while len(batch) < batch_size:
      try:
        entry = iterator.next()
        
      except StopIteration:
        entry = None
        
      if entry is None:
        break  #End of file
      
      batch.append(entry)
      
    ##end while
      
    if batch:
      yield batch
      
  ##end while
      
##end batch_iterator()


#*******************************************************************************
# @func_name: delete_file()
# @desc: Deletes specified file
# @return: None
# @reference: http://www.cyberciti.biz/faq/howto-python-delete-files/
#*******************************************************************************
def delete_file(filename):
  
  try:
    os.remove(filename)
    
  except OSError, e:  ## if failed, report it back to the user ##
    print ('Error: %s - %s.' % (e.filename,e.strerror))

##end delete_file()


#*******************************************************************************
# @func_name: make_dir()
# @desc: Creates directory if it does not exist
# @return: None
#*******************************************************************************
def make_dir(folder_name):
  
  try:
    os.makedirs(folder_name)
      
  except OSError:
    if not os.path.isdir(folder_name):
      raise
        
##end make_dir()


#*******************************************************************************
# @func_name: delete_temp_folder()
# @desc: Deletes temporary file directory
# @return: None
#*******************************************************************************
def delete_temp_folder(folder):

  # ensure that folder is empty
  for the_file in os.listdir(folder):
    file_path = os.path.join(folder, the_file)
    
    try:
        if os.path.isfile(file_path):
            os.unlink(file_path)
            
    except Exception, e:
        print e

  # deletes the folder
  try:
    os.rmdir(folder)
    
  except OSError as err:
    if err.errno == errno.ENOTEMPTY:
      print 'Cannot delete %s. Folder is not empty.' % (folder)
    else:
      print err

##end delete_temp_folder()


#*******************************************************************************
# @func_name: main()
# @desc: Program starts here
#*******************************************************************************
def main(argv):

  if len(argv):
    # Specimen name
    specimen_name = argv[0]
  else:
    # missing argument
    print 'Specimen name argument missing'
    return 1

  REPORT_LOG.append(specimen_name)

  # number of records saved per file
  records_per_file = 50000

  # current working directory
  cwd = os.getcwd()

  # BLAST database directory
  db_path = cwd + '/Data/Database/'
  
  # BLAST database name
  blast_db_name = db_path + 'ARDB'

  # directory of specimen where FASTA file to BLAST is located
  unmapped_path = cwd + '/Data/Library/' + specimen_name + '/Unmapped/'
  
  # path of specimen where to save BLAST result
  xml_path = unmapped_path + 'BLAST_Output/xml/'


  #<<=== USED FOR DATABASE MAINTENANCE
  # unmapped_path = db_path
  # xml_path = cwd + '/Data/Library/' + specimen_name + '/BLAST_Output/xml/'


  # creates directory if it does not exist
  make_dir(xml_path)

  # specimen's complete path and name
  specimen_fasta = unmapped_path + specimen_name + '.fasta'

  # temporary file directory
  temp_path = xml_path + 'tmp/'
  # creates temp folder
  make_dir(temp_path)

  # Start Timer for total BLAST
  startTime = datetime.datetime.now()
  time_stamp = startTime.strftime('%Y%m%d_%H%M')
  print 'Start time: ', startTime.strftime('%H:%M:%S')
  print ''
  REPORT_LOG.append('Start time: %s' % (startTime.strftime('%H:%M:%S')))

  # parses specimen FASTA file
  record_iter = SeqIO.parse(open(specimen_fasta),'fasta')

  # iterates through the specimen FASTA file by creating separate files by specifying
  # number of records per file
  for i, batch in enumerate(batch_iterator(record_iter, records_per_file)):

    # temporary FASTA file with specific number of records
    fasta_file = temp_path + specimen_name + '_group_%i.fasta' % (i+1)

    # opens file stream of temporary FASTA file
    fasta_fp = open(fasta_file, 'w')

    # writes the record batch on temporary FASTA file and get a count
    count = SeqIO.write(batch, fasta_fp, 'fasta')

    # closes temporary FASTA file stream
    fasta_fp.close()
    # prompts user of number of records written on temp FASTA file
    print 'Wrote %i records to %s\n' % (count, fasta_file)
    REPORT_LOG.append('Wrote %i records to %s' % (count, fasta_file))

    # filename for xml BLAST result
    group_file = xml_path + 'group_%i.xml' % (i+1)
    
    # use BLAST command line operation to get result for FASTA file
    blastn_cline = NcbiblastnCommandline(cmd=BLAST_PATH + 'blastn',
                      out=group_file, outfmt=5, query=fasta_file,
                      db=blast_db_name, evalue=0.001, strand='both',
                      num_threads=2, task='megablast')

    # BLAST start timer
    startTime_blast = datetime.datetime.now()
    print 'BLAST Group %i Initiated' % (i+1)
    REPORT_LOG.append('BLAST Group %i Initiated' % (i+1))

    # BLAST fasta file
    stdout, stderr = blastn_cline()
    
    # delete temporary FASTA file
    delete_file(fasta_file)

    # BLAST end timer
    endTime_blast = datetime.datetime.now()
    print 'BLAST Group %i Completed' % (i+1)
    print 'Time elapsed:', endTime_blast - startTime_blast
    print ''
    REPORT_LOG.append('BLAST Group %i Completed' % (i+1))
    REPORT_LOG.append('Time elapsed: %s' % (endTime_blast - startTime_blast))

  ##end for

  # deletes the temporary folder
  delete_temp_folder(temp_path)

  # deletes specimen fasta file
  delete_file(specimen_fasta)

  # End Timer for total BLAST
  endTime = datetime.datetime.now()
  print 'End time: ', endTime.strftime('%H:%M:%S')
  print 'Time elapsed:', endTime - startTime
  REPORT_LOG.append('End time: %s' % (endTime.strftime('%H:%M:%S')))
  REPORT_LOG.append('Time elapsed: %s' % (endTime - startTime))

  # Creates a Report log file on Reports Folder
  report_fp = open(cwd + '/Data/Library/' + specimen_name + '/Reports/Step2.log.txt', 'w')
  report_fp.write('\n'.join(REPORT_LOG))
  report_fp.close()

##end main()

  

if __name__ == '__main__':
  main(sys.argv[1:])
