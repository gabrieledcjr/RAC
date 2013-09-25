#!/usr/bin/python -tt
'''
Version Created on June 4, 2013

@author: Gabriel de la Cruz

@description: Parses all the hits from the xml BLAST output into the
specimen table in the RAC MySQL database
'''

# use for command line argument manipulation
import sys

import gc

# use for directory path and file manipulation
import os
import fnmatch

# use for string manipulation
import string

# use to get system date and time
import datetime

# [BioPython] third party modules use to manipulate sequence files
from Bio.Blast import NCBIXML

# python file that houses functions to connect to MySQL
# and run SQL statements
import mysqlc


# Report log List
REPORT_LOG = []


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
    # specimen name
    specimen_name = argv[0]
  else:
    # missing argument
    print 'Specimen name argument missing'
    return 1

  REPORT_LOG.append(specimen_name)
  
  # name of table to be added on the database
  table_name = specimen_name


  # MySQL connection parameters
  db_username = 'racuser'     
  db_password = 'racuser123'
  db_host = 'localhost'
  db_name = 'RAC'


  # SQL statement to create specimen table
  RAC_TABLE = (
    "CREATE TABLE %s (" % (table_name) +
    " rec_no int NOT NULL AUTO_INCREMENT,"
    " hit_id int NOT NULL,"
    " alg_len int(5) NOT NULL,"
    " query varchar(100) NOT NULL,"
    " query_len int(4) NOT NULL,"
    " percent_id float(5,2) NOT NULL,"
    " hsp_exp float NOT NULL,"
    " hsp_scr float NOT NULL,"
    " hsp_bits float NOT NULL," 
    " hsp_qry_st int(5) NOT NULL,"
    " hsp_qry_en int(5) NOT NULL,"
    " hsp_sbj_st int(5) NOT NULL,"
    " hsp_sbj_en int(5) NOT NULL,"
    " hsp_alg_len int(4) NOT NULL,"
    " hsp_gaps int(11) NOT NULL,"
    " hsp_qry varchar(255) NOT NULL,"
    " hsp_mtch varchar(255) NOT NULL,"
    " hsp_sbj varchar(255) NOT NULL,"
    " PRIMARY KEY (rec_no),"
    " FOREIGN KEY (hit_id) REFERENCES ARDB(hit_id)"
    ") ENGINE=InnoDB")

  # SQL statement to add values to the specimen table
  ADD_REC = ("INSERT INTO %s " % (table_name) +
             "(hit_id, alg_len, query, query_len, percent_id, hsp_exp,"
             " hsp_scr, hsp_bits, hsp_qry_st, hsp_qry_en, hsp_sbj_st,"
             " hsp_sbj_en, hsp_alg_len, hsp_gaps, hsp_qry, hsp_mtch, hsp_sbj) "
             "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s,"
             " %s, %s, %s, %s, %s, %s, %s)")


  # current working directory
  cwd = os.getcwd()

  # set path to xml directory of specific specimen
  xml_path = cwd + '/Data/Library/' + specimen_name + '/Unmapped/BLAST_output/xml/'

  # XML group filename template
  xml_file = 'group_'

  # total number of XML files in the XML directory
  total_xml_files = len(fnmatch.filter(os.listdir(xml_path), '*.xml'))
  

  ## Start Timer
  startTime = datetime.datetime.now()
  print 'Start time: ', startTime.strftime('%H:%M:%S')
  REPORT_LOG.append('Start time: %s' % (startTime.strftime('%H:%M:%S')))

  # connect to MySQL database
  db = mysqlc.connect_db(db_username, db_password, db_host, db_name)

  # create specimen table
  mysqlc.create_table(db, table_name, RAC_TABLE)
    

  # Go through all XML files in the xml directory
  for file_num in range(1, (total_xml_files + 1)):
    
    # xml group Start Timer
    fileStartTime = datetime.datetime.now()
    print xml_file + str(file_num) + '.xml Started:', fileStartTime.strftime('%H:%M:%S')
    REPORT_LOG.append(xml_file + str(file_num) + '.xml Started: %s' % (fileStartTime.strftime('%H:%M:%S')))

    # counters
    count_que = 0 # counts # of queries
    count_alg = 0 # counts # of alignments
    count_hit = 0 # counts # of hits

    # opens XML group file stream
    result_handle = open(xml_path + xml_file + str(file_num) + '.xml')

    # uses Biopython to parse XML file
    blast_records = NCBIXML.parse(result_handle)

    # traverse through all the BLAST result in the xml file
    for blast_record in blast_records:

      # counter for queries
      count_que += 1
      
      # checks if blast record has any alignment
      if blast_record.alignments:

        # counter for alignments
        count_alg += 1

        # traverse through all alignments
        for alignment in blast_record.alignments:

          # traverse through all high-scoring segment pairs for each alignment
          for hsp in alignment.hsps:

            # compute for percent of similarity of length
            percent_len = hsp.score / float(blast_record.query_letters)
            # compute for percent identity
            percent_id = hsp.match.count('|') / float(hsp.align_length)

            # save values of the hit in a tuple
            tuple_values = (string.split(alignment.hit_id, '|')[2], alignment.length,
                            blast_record.query, blast_record.query_letters, (percent_id * 100),
                            hsp.expect, hsp.score, hsp.bits, hsp.query_start,
                            hsp.query_end, hsp.sbjct_start, hsp.sbjct_end,
                            hsp.align_length, hsp.gaps, hsp.query, hsp.match, hsp.sbjct)

            # adds hit to database
            mysqlc.insert_query(db, ADD_REC, tuple_values)

            # counts hits
            count_hit += 1

            # UNCOMMENT if you want queries prompted on screen
            # output = [blast_record.query,
            #          str(blast_record.query_letters),
            #          str(alignment.hit_id),
            #          str(alignment.length),
            #          str(hsp.align_length),
            #          '%.2f' % (percent_id) + '%']
            # print ' \t'.join(output)
            
          ##end for
                       
        ##end for
              
      ##end if
      
    ##end for

    
    # closes connections
    blast_records.close()
    result_handle.close()

    # XML group End Timer
    fileEndTime = datetime.datetime.now()
    print xml_file + str(file_num) + '.xml Completed:', fileEndTime.strftime('%H:%M:%S')
    print 'Time elapsed:', fileEndTime - fileStartTime
    REPORT_LOG.append(xml_file + str(file_num) + '.xml Completed: %s' % (fileEndTime.strftime('%H:%M:%S')))
    REPORT_LOG.append('Time elapsed: %s' % (fileEndTime - fileStartTime))
   
    # prompts user of report for each xml file
    print str(count_que) + ' queries processed [' + str(count_alg) + ' alignments, ' + str(count_hit) + ' hits]\n'
    REPORT_LOG.append(str(count_que) + ' queries processed [' + str(count_alg) + ' alignments, ' + str(count_hit) + ' hits]')

    # clears memory every 10 files using purge command line (PLATFORM DEPENDENT)
    if file_num % 10 == 0:
      clear_memory()
    
    #break  # <== THIS BREAK CAUSES THE PROGRAM TO ONLY PROCESS ONE XML FILE

  ##end for
    

  # disconnect from MySQL database
  db.close()

  # clears memory using purge command line (PLATFORM DEPENDENT)
  clear_memory()

  # End Timer
  endTime = datetime.datetime.now()
  print 'End time: ', endTime.strftime('%H:%M:%S')
  print 'Time elapsed:', endTime - startTime
  REPORT_LOG.append('End time: %s' % (endTime.strftime('%H:%M:%S')))
  REPORT_LOG.append('Time elapsed: %s' % (endTime - startTime))

  report_fp = open(cwd + '/Data/Library/' + specimen_name + '/Reports/Step3.log.txt', 'w') 
  report_fp.write('\n'.join(REPORT_LOG))
  report_fp.close()
##end main()


  
if __name__ == '__main__':
  main(sys.argv[1:])
