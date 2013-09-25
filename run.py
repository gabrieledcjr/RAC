#!/usr/bin/python -tt
'''
Created on August 21, 2013

@author: Gabriel de la Cruz

@description: A single command that goes through all four processes
to get to the end product
'''

# use for command line argument
import sys

# use for directory path and file manipulation
import os


#******************************************************************************
# @func_name: main()
# @desc: Program starts here
#******************************************************************************
def main(argv):

  if len(argv) == 0:
    # missing argument
    print 'Specimen name argument(s) missing'
    return 1
  else:
    # specimen name
    for specimen_name in argv:
  
      step1_command = 'python bam_to_unmapped_fasta.py ' + specimen_name
      step2_command = 'python blast_fasta.py ' + specimen_name
      step3_command = 'python parse_xml_specimen.py ' + specimen_name
      step4_command = 'python assemble_hit_files.py ' + specimen_name

      print 'STRAIN: ' + specimen_name
      print ''
      
      print '> FILTER UNMAPPED READS FROM BAM FILE'
      os.system(step1_command)
      print ''
      
      print '> FASTA FILE GETS BLASTED AGAINST ARDB'
      os.system(step2_command)
      print ''
      
      print '> PARSE VALUES FROM XML FILES INTO MYSQL TABLE'
      os.system(step3_command)
      print ''
      
      print '> CREATE HIT FILES FOR ASSEMBLY'
      os.system(step4_command)
      print ''

      print '-Process Completed for ' + specimen_name + '-'
      print ''

##end main()

if __name__ == '__main__':
  main(sys.argv[1:])
