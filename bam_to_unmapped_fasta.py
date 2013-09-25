#!/usr/bin/python -tt
'''
Created on May 21, 2013

@author: de la Cruz

@description: Filtering out unmapped reads from MiSeq BAM file using third
party module Pysam. Filtered unmapped reads will be written on another BAM file.
BAM file will then be sorted, and the result will be converted into FASTA file.
'''

# use to get system date and time
import datetime

# use for directory path and file manipulation
import os

# use for command line argument
import sys

# use for string manipulation
import string

# [PySam] third party module use to manipulate SAM/BAM file format
import pysam


# location of samtools directory for command line operations
SAMTOOLS_PATH = os.getcwd() + '/Util/samtools-0.1.19/'

# Report log List
REPORT_LOG = []

#*******************************************************************************
# @func_name: filter_unmapped_reads()
# @desc: Filter unmapped reads from the aligned BAM file. Will require the
#        indexed BAM file from the same directory of the BAM file.
# @return: None
#*******************************************************************************
def filter_unmapped_reads(bam_aligned, bam_unmapped):

  # opens read and write files
  bam_aligned_fp = pysam.Samfile(bam_aligned, 'rb')
  bam_unmapped_fp = pysam.Samfile(bam_unmapped, 'wb', template=bam_aligned_fp)

  # If until_eof is given, all reads from the current file position
  # will be returned in order as they are within the file.
  # Using this option will also fetch unmapped reads.
  for read in bam_aligned_fp.fetch(until_eof = True):
    
    if read.is_unmapped:
      bam_unmapped_fp.write(read)
      
  ##end for

  # closes file streams
  bam_unmapped_fp.close()
  bam_aligned_fp.close()

##end filter_unmapped_reads()


#*******************************************************************************
# @func_name: sort_bam_file()
# @desc: Sorts specified BAM file
# @return: None
#*******************************************************************************
def sort_bam_file(unmapped, unmapped_sorted):

  # pysam.sort() unsupported in Python 2.7.5
  # pysam.sort('-n', unmapped, unmapped_sorted)
  # -n sorts it by read names
  command = SAMTOOLS_PATH + 'samtools sort -n ' + unmapped + ' ' + unmapped_sorted
  os.system(command)

##end sort_bam_file()


#*******************************************************************************
# @func_name: create_bam_index()
# @desc: Creates a BAM index file 
# @return: None
#*******************************************************************************
def create_bam_index(bam_file):

  # create index BAM file
  # pysam.index() unsupported in Python 2.7.5
  # pysam.index(bam_file, bam_file + '.bai')
  command = SAMTOOLS_PATH + 'samtools index ' + bam_file
  os.system(command)

##end create_bam_index()


#*******************************************************************************
# @func_name: bam_to_fasta()
# @desc: Converts BAM file to FASTA file
# @return: None
#*******************************************************************************
def bam_to_fasta(bam_file, fasta_file):

  # opens read and write files
  bam_fp = pysam.Samfile(bam_file, 'rb')
  fasta_fp = open(fasta_file, 'w')

  # If until_eof is given, all reads from the current file position
  # will be returned in order as they are within the file.
  # Using this option will also fetch unmapped reads.
  for read in bam_fp.fetch(until_eof = True):
    
    if read.is_unmapped:
      query_name = '>' + read.qname
      
      if read.is_read1:
        query_name += '/1'
      elif read.is_read2:
        query_name += '/2'
      
      fasta_fp.write(query_name + '\n')
      fasta_fp.write('\n'.join(split_str_into_len(read.query)) + '\n')

  bam_fp.close()
  fasta_fp.close()
  
##end bam_to_fasta()


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
# @func_name: clear_memory()
# @desc: Clears inactive memory using purge command line (PLATFORM DEPENDENT)
# @return: None
#*******************************************************************************
def clear_memory():
  os.system('purge')

##end clear_memory()

#*******************************************************************************
# @func_name: main()
# @desc: Program starts here
#*******************************************************************************
def main(argv):

  if len(argv):
    # Specimen name
    specimen_name = argv[0]  #<<<=== CHANGE SPECIMEN NAME
  else:
    # missing argument
    print 'Specimen name argument missing'
    return 1

  REPORT_LOG.append(specimen_name)

  # current working directory
  cwd = os.getcwd()


  # source path for specimen
  source_path = cwd + '/Data/Library/' + specimen_name + '/Source/'

  # unmapped path for specimen
  unmapped_path = cwd + '/Data/Library/' + specimen_name + '/Unmapped/'
  # create unmapped directory if it does not exist
  make_dir(unmapped_path)
  
  
  # aligned BAM filename
  bam_aligned = source_path + specimen_name + '.bam'
  # unmapped BAM filename
  bam_unmapped = unmapped_path + specimen_name + '.bam'


  ### Filter unmapped reads from BAM file using Pysam
  print 'Filter unmapped reads from BAM file'
  REPORT_LOG.append('Filter unmapped reads from BAM file')

  # Start Timer for BAM filtering
  startTime = datetime.datetime.now() 
  print 'Start time: ', startTime.strftime('%H:%M:%S')
  REPORT_LOG.append('Start time: %s' % (startTime.strftime('%H:%M:%S')))
  
  # filtering unmapped reads from BAM files
  filter_unmapped_reads(bam_aligned, bam_unmapped)

  # End Timer
  endTime = datetime.datetime.now()
  print 'End time: ', endTime.strftime('%H:%M:%S')
  print 'Time elapsed:', endTime - startTime
  REPORT_LOG.append('End time: %s' % (endTime.strftime('%H:%M:%S')))
  REPORT_LOG.append('Time elapsed: %s' % (endTime - startTime))
  ### End Filter unmapped reads


  ### Create BAM index file
  print '\nCreate BAM index file'
  REPORT_LOG.append('Create BAM index file')

  # Start Timer for BAM filtering
  startTime = datetime.datetime.now()
  print 'Start time: ', startTime.strftime('%H:%M:%S')
  REPORT_LOG.append('Start time: %s' % (startTime.strftime('%H:%M:%S')))
  
  # create indexed BAM file
  create_bam_index(bam_unmapped)

  # End Timer
  endTime = datetime.datetime.now()
  print 'End time: ', endTime.strftime('%H:%M:%S')
  print 'Time elapsed:', endTime - startTime
  REPORT_LOG.append('End time: %s' % (endTime.strftime('%H:%M:%S')))
  REPORT_LOG.append('Time elapsed: %s' % (endTime - startTime))
  ### End Create BAM index file


  # unmapped sorted BAM filename
  bam_unmapped_sorted = unmapped_path + specimen_name + '_sorted'  


  ### Sort unmapped BAM file
  print '\nSort unmapped BAM file'
  REPORT_LOG.append('Sort unmapped BAM file')

  # Start Timer for BAM filtering
  startTime = datetime.datetime.now()
  print 'Start time: ', startTime.strftime('%H:%M:%S')
  REPORT_LOG.append('Start time: %s' % (startTime.strftime('%H:%M:%S')))
  
  # sort unmapped BAM file
  sort_bam_file(bam_unmapped, bam_unmapped_sorted)

  # End Timer
  endTime = datetime.datetime.now()
  print 'End time: ', endTime.strftime('%H:%M:%S')
  print 'Time elapsed:', endTime - startTime
  REPORT_LOG.append('End time: %s' % (endTime.strftime('%H:%M:%S')))
  REPORT_LOG.append('Time elapsed: %s' % (endTime - startTime))
  ### End Sort unmapped BAM file


  # unmapped sorted BAM filename
  bam_unmapped_sorted += '.bam'

  # delete unsorted BAM file and index
  delete_file(bam_unmapped)
  delete_file(bam_unmapped + '.bai')

  # FASTA filename
  fasta_file = unmapped_path + specimen_name + '.fasta'

  
  ### Convert sorted BAM file to FASTA
  print '\nConvert sorted unmapped BAM file to FASTA format'
  REPORT_LOG.append('Convert sorted unmapped BAM file to FASTA format')

  # Start Timer for BAM filtering
  startTime = datetime.datetime.now() 
  print 'Start time: ', startTime.strftime('%H:%M:%S')
  REPORT_LOG.append('Start time: %s' % (startTime.strftime('%H:%M:%S')))
  
  # converts sorted SAM/BAM file to FASTA file
  bam_to_fasta(bam_unmapped_sorted, fasta_file)

  # End Timer
  endTime = datetime.datetime.now()
  print 'End time: ', endTime.strftime('%H:%M:%S')
  print 'Time elapsed:', endTime - startTime
  REPORT_LOG.append('End time: %s' % (endTime.strftime('%H:%M:%S')))
  REPORT_LOG.append('Time elapsed: %s' % (endTime - startTime))
  ### End Convert sorted BAM file to FASTA


  # report path for specimen
  report_path = cwd + '/Data/Library/' + specimen_name + '/Reports/'
  # creates folder if it does not exist
  make_dir(report_path)

  # Create a Report log file on Reports Folder
  report_fp = open(report_path + 'Step1.log.txt', 'w')
  report_fp.write('\n'.join(REPORT_LOG))
  report_fp.close()


  # clears memory using purge command line (PLATFORM DEPENDENT)
  clear_memory()
  
##end main()

  

if __name__ == '__main__':
  main(sys.argv[1:])
