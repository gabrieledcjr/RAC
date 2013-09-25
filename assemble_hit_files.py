#!/usr/bin/python -tt
'''
Created on June 6, 2013

@author: Gabriel de la Cruz

@description: Generate FASTQ files to be assembled for all AR group hits.
Program will have to query all hits for every AR group filtered against
two thresholds: total number of hits per AR group and hsp_exp equals to 0.
Once filtered, program will have to retrieve the original sequence from
the sorted BAM file of unmapped reads for every hit in the AR group. These
sequences will then be saved in FASTQ files which will then be imported to
CLC Assembly to mapped against the AR Group and locate SNPs.
'''

import sys

# use for directory path and file manipulation
import os
import glob

# use to get system date and time
import datetime

# python file that houses functions to connect to MySQL
# and run SQL statements
import mysqlc

# [PySam] third party module use to manipulate SAM/BAM file format
import pysam


# path for bamUtil program
BAMUTIL_PATH = os.getcwd() + '/Util/bamUtil_1.0.9/bamUtil/bin/'


#*******************************************************************************
# @func_name: create_bam_index_list()
# @desc: Creates a list that indexed file positions of the BAM file to be used
#        to speed searching of a query
# @return: None
#*******************************************************************************
def create_bam_index_list(bam_file, bam_index):

  # counts number of queries
  ctr = 0

  # traverses through all queries in the BAM file
  while True:
    
    try:
      # reads each query in the bam file
      aligned_read = bam_file.next()
      
    except StopIteration as err:
      try:
        break  # end of file
      
      except:
        print err
    
    if (ctr % 5000) == 0:
      header_list = aligned_read.qname.split(':')
      header_list = header_list[-4:]
      header_list = map(int, header_list)
      header_list.append(bam_file.tell())
      bam_index.append(header_list)

    ctr += 1

  ##end while

##end create_bam_index_list()
    

#*******************************************************************************
# @func_name: get_file_position()
# @desc: Searches the query header using binary search against the indexed list
#        of file positions and returns the file position when located
# @return: File position of a query
#*******************************************************************************
def get_file_position(header_list, bam_index_list):
  start = 0
  end = len(bam_index_list) - 1
  while True:
    
    mid = (start + end) / 2
    #print "s: %d, e: %d, m: %d" % (start, end, mid)

    if (header_list < bam_index_list[mid][:-1]):
      end = mid
    elif (header_list > bam_index_list[mid][:-1]):
      start = mid
    else:
      # mid - 1, since it's equal, you have to start to the file position before it
      # NEEDS TO BE FIXED FOR EFFICIENCY
      return bam_index_list[mid - 1][-1]

    if start + 1 == end:
      if (header_list > bam_index_list[end][:-1]):
        return bam_index_list[end][-1]
      else:
        return bam_index_list[start][-1]
      
  ##end while

##end get_file_position()
      

#*******************************************************************************
# @func_name: bam_to_fastq()
# @desc: Using bamUtil to convert SAM/BAM file to FASTQ format. bamUtil program
#        will generate 3 FASTQ files: 1st pair read, 2nd pair read and reads
#        with no pair. NOTE: It is important that BAMUTIL_PATH has the correct
#        location of the software
# @return: None
#*******************************************************************************
def bam_to_fastq(bamFile, fastqFile):

  command = BAMUTIL_PATH + "./bam bam2FastQ --in " + bamFile + " --outBase " + fastqFile + " --readname"
  
  os.system(command)

##end bam_to_fastq()


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
# @func_name: delete_folder()
# @desc: Deletes file directory
# @return: None
#*******************************************************************************
def delete_folder(folder):

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
    print 'Deleting %s...' %(folder)
    
  except OSError as err:
    if err.errno == errno.ENOTEMPTY:
      print 'Cannot delete %s. Folder is not empty.' % (folder)
    else:
      print err

##end delete_folder()
    

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

  # MySQL table name in RAC database
  table_name = specimen_name

  # num of hits much reach threshold to be considered for assembly
  threshold_hits = 100

  # only hits with hsp_exp equal to 0 are assembled
  threshold_string = 'hsp_exp = 0'


  # MySQL connection parameters
  db_username = 'racuser'     
  db_password = 'racuser123'
  db_host = 'localhost'
  db_name = 'RAC'


  # current working directory
  cwd = os.getcwd()

  # output directory for the hit files
  output_path = cwd + '/Data/Library/' + specimen_name + '/Hit_files/'
  # creates directory if it does not exist
  make_dir(output_path)

  
  # connect to MySQL database
  db = mysqlc.connect_db(db_username, db_password, db_host, db_name)


  ### Create index list of bam file
  print '\nCreate Indexed List of BAM file'

  # Start Timer for BAM to FASTQ conversion
  startTime = datetime.datetime.now() 
  print 'Start time: ', startTime.strftime('%H:%M:%S')

  # open bam file stream to sorted & unmapped specimen bam file
  bam_in_fp = pysam.Samfile(cwd + '/Data/Library/' + specimen_name + '/Unmapped/' + specimen_name + '_sorted.bam', 'rb')

  # create empty bam index list
  bam_index_list = []
  # Create an index list of sorted & unmapped specimen bam file
  create_bam_index_list(bam_in_fp, bam_index_list)

  # End Timer
  endTime = datetime.datetime.now()
  print 'End time: ', endTime.strftime('%H:%M:%S')
  print 'Time elapsed:', endTime - startTime
  ### End Create index list of bam file
 

  # select all AR groups in MySQL specimen table
  query = ("SELECT COUNT(*) AS num_hits, AR_groups_table.group_name AS gp_name "
           "FROM %s, AR_groups_table " % (table_name) +
           "WHERE %s.hit_id = AR_groups_table.hit_id " % (table_name) +
           "GROUP BY AR_groups_table.group_name "
           "ORDER BY num_hits DESC")
  # execute select query
  cursor_ar_group = mysqlc.select_query(db, True, query)


  for (num_hits, group_name) in cursor_ar_group:

    # ensure that num of hits is over hits threshold
    if num_hits >= threshold_hits:
      
      # output directory for the AR group
      output_group_path = output_path + '/' + '_'.join(group_name.split(' ')) + '/'
      # creates directory if it does not exist
      make_dir(output_group_path)
      

      # select all AR genes in MySQL of particular AR group
      query = ("SELECT AR_groups_table.hit_id AS hit_id, AR_groups_table.alg_len AS alg_len, "
               "ARDB.resistance_type AS res_type, ARDB.phenotype AS phenotype, "
               "ARDB.description AS description, ARDB.sequence AS sequence "
               "FROM AR_groups_table, ARDB "
               "WHERE AR_groups_table.group_name LIKE '%s' " % (group_name) +
               "AND AR_groups_table.hit_id LIKE ARDB.hit_id "
               "ORDER BY AR_groups_table.alg_len DESC")
      # execute select query
      cursor_ar_gene = mysqlc.select_query(db, True, query)


      # output directory for the reference genes in the AR group
      output_group_ref_path = output_group_path + '/Reference/'
      # creates directory if it does not exist
      make_dir(output_group_ref_path)


      # creates empty list for hit ids for AR group
      hit_ids = []

      # prompt user on # of hits for AR group 
      print "\n%s \t [%d hits]" % (group_name, num_hits)

      # traverse to all genes in the AR group
      for (hit_id, alg_len, res_type, phenotype, description, sequence) in cursor_ar_gene:
        
        # prompt user on alignment length and hit id for every gene in the group
        print "> gnl|ARDB|%d\t%d" % (hit_id, alg_len)

        # appends hit id in thelist
        hit_ids.append(hit_id)

        # gene filename, convert spaces to underscores
        gene_filename = '_'.join(res_type.split(' '))

        # reference gene file stream
        ref_gene_fp = open(output_group_ref_path + gene_filename + '.fasta', 'w')
        
        # >gnl|ARDB|<hit_id> <resistance_type>|<phenotype>
        ref_gene_fp.write('>%s %s|%s|[gnl|ARDB|%i]\n' % (res_type, phenotype, description, hit_id))
        # ATGCCCATAGCT...TAA
        ref_gene_fp.write('\n'.join(split_str_into_len(sequence)) + '\n')

        # closes file stream
        ref_gene_fp.close()

      ##end for
        

      # sql 'LIKE' string
      like_statement = ""
      
      for index in range(0, len(hit_ids)):
        
        if index > 0:
          like_statement += "OR "
        like_statement += "hit_id LIKE '" + str(hit_ids[index]) + "' "
        
      ##end for


      # select all queries for selected AR genes of the selected AR group
      query = ("SELECT query "
               "FROM %s " % (table_name) +
               "WHERE (%s) AND (%s)"% (like_statement, threshold_string) +
               "GROUP BY query "
               "ORDER BY query ASC")
      # execute select query
      cursor_queries = mysqlc.select_query(db, True, query)


      # bam out file stream
      bam_out_fp = pysam.Samfile(output_group_path + '_'.join(group_name.split(' ')) + '.bam', 'wb', template=bam_in_fp)


      ### Create BAM file
      print '\nCreate BAM file'

      # Start Timer
      startTime = datetime.datetime.now()
      print 'Start time: ', startTime.strftime('%H:%M:%S')

      # initialize record counter
      count_record = 0

      # traverse to all queries
      for fasta_record in cursor_queries:

        # iterate record counter
        count_record += 1

        # takes out the /1 or /2 
        header = fasta_record[0][:-2]

        # splits the header by ':'
        header_list = header.split(':')

        # get rid of unnecessary part of the query header
        header_list = header_list[-4:]

        # converts all items in the list into int type
        header_list = map(int, header_list)

        # locates file position of the query header list
        temp = get_file_position(header_list, bam_index_list)

        if temp == 0:
          bam_in_fp.reset()
        else:
          bam_in_fp.seek(temp)


        # loops infinitely until query record is found
        while True:
          
          try:
            # sets file pointer to next record
            read = bam_in_fp.next()
            
          except StopIteration as err:
            try:
              # NOTE: should NOT print at all if program is correct
              print header_list
              print bam_index_list
              print temp
              print 'ERROR IN ALGORITHM: Query not found!'
              break
            
            except:
              print err

          # set qname to read query name
          qname = read.qname

          # identify is read is pair 1 or 2
          if read.is_read1:
            qname += '/1'
          else:
            qname += '/2'

          # check if query name is record being searched
          if qname == fasta_record[0]:
            # write on file if record located
            bam_out_fp.write(read)
            #print str(count_record) + ": " + qname + " found!"
            break
          
        ##end while
      ##end for

      # no record was saved
      if count_record == 0:
        delete_folder(output_group_ref_path)
        delete_folder(output_group_path)

      # close file stream
      bam_out_fp.close()

      # close cursor stream
      cursor_ar_gene.close()

      # End Timer
      endTime = datetime.datetime.now()
      print 'Time elapsed: ', endTime - startTime
      print 'End time: ', endTime.strftime('%H:%M:%S')      
      print "%d sequences saved" % (count_record)
      ### End Create BAM file
      
    ##end if
  ##end for

  # close file stream
  bam_in_fp.close()

  # close cursor
  cursor_ar_group.close()

  # close db connection
  db.close()
  

  # loop through all AR group folders
  for group_folder in glob.glob(output_path + '*'):
    
    # loop through all bam files in the AR group folder
    for bam_file in glob.glob(group_folder + '/' + '*.bam'):

      # prompts user of bam filename
      print '\n' + bam_file.split('/')[-1]
      
      # prompt user conversion has started
      print 'Converting to FASTQ format...'

      # Start Timer for BAM to FASTQ conversion
      startTime = datetime.datetime.now()  
      print 'Start time: ', startTime.strftime('%H:%M:%S')

      # convert bam file to fastq files
      bam_to_fastq(bam_file, bam_file[:-4])

      # deletes bam file
      delete_file(bam_file)

      # End Timer
      endTime = datetime.datetime.now()
      print 'End time: ', endTime.strftime('%H:%M:%S')
      print 'Time elapsed:', endTime - startTime
    
    ##end for
  ##end for
  os.system('purge')    
##end main()



if __name__ == '__main__':
  main(sys.argv[1:])
