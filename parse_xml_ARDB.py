#!/usr/bin/python -tt
'''
Created on May 27, 2013

@author: Gabriel de la Cruz
'''

# use for directory path and file manipulation
import os

# use for string manipulation
import string

# use to get system date and time
import datetime

# Third-party module MySQL Connector/Python
# use to connect to MySQL database
import mysql.connector
from mysql.connector import errorcode

# Third-party module BioPython
# use to parse BLAST xml file output
from Bio.Blast import NCBIXML

###############  FUNCTIONS  ################

def item_to_db(data_rec):
  cursor.execute(ADD_REC, data_rec)
  db.commit()

def counter():
  global count_seq
  count_seq += 1
  #print count_seq

############################################

################  GLOBAL  ##################

table_name = 'AR_groups_table'

db_username = 'racuser'     
db_password = 'racuser123'
db_host = 'localhost'
db_name = 'RAC'

query_name = 'ARDB'

xml_path = os.getcwd() + '/../Library/' + query_name + '/BLAST_output/xml/'
xml_file = 'group_'

#blast_dict = {}
#seq_dict = {}

count_seq = 0

start_xml_file_num = 1
end_xml_file_num = 2

threshold = 0.3

final_list = []
sets_list = []

group_name = 'AR Group'

############################################


## Start Timer
startTime = datetime.datetime.now()
print 'Start time: ', startTime.strftime('%H:%M:%S')

RAC_TABLE = (
  "CREATE TABLE %s (" % (table_name) +
  " rec_no int NOT NULL AUTO_INCREMENT,"
  " group_name varchar(12) NOT NULL,"
  " hit_id int NOT NULL,"
  " alg_len int(5) NOT NULL,"
  " PRIMARY KEY (rec_no),"
  " FOREIGN KEY (hit_id) REFERENCES ARDB(hit_id)"
  ") ENGINE=InnoDB")

ADD_REC = ("INSERT INTO %s " % (table_name) +
           "(group_name, hit_id, alg_len) "
           "VALUES (%s, %s, %s)")

DROP_TABLE = ("DROP TABLE %s " % (table_name))

try:
  # Open database connection
  db = mysql.connector.connect(user=db_username,password=db_password,host=db_host,database=db_name)
  cursor = db.cursor()

except mysql.connector.Error as err:
  if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
    print("Something is wrong with your username or password")
  elif err.errno == errorcode.ER_BAD_DB_ERROR:
    print("Database does not exists")
  else:
    print(err)


try:
  print("Creating %s Table" % (table_name))
  cursor.execute(RAC_TABLE)
except mysql.connector.Error as err:
  if err.errno == errorcode.ER_TABLE_EXISTS_ERROR:
    try:
      print("%s already exists." % (table_name))
      print("Deleting table %s!" % (table_name))
      cursor.execute(DROP_TABLE)
      print("Recreating %s Table" % (table_name))
      cursor.execute(RAC_TABLE)
    except mysql.connector.Error as err:
      print(err)
    else:
      print("OK")
  else:
    print(err)
else:
  print("OK")

# group_1 to group_##
# Go through all group xml files
for file_num in range(start_xml_file_num, end_xml_file_num):
  result_handle = open(xml_path + xml_file + str(file_num) + '.xml')
  blast_records = NCBIXML.parse(result_handle)

  for blast_record in blast_records:
    
    #print "Query: " + blast_record.query
    #print "Query Length: " + str(blast_record.query_letters)
    temp_list = []
    for alignment in blast_record.alignments:
      for hsp in alignment.hsps:
        # (length of matching sequence) * (BLAST identity score) / 
        # (length of reference gene + length of matching sequence gene)
        equation = (hsp.align_length * ((hsp.match.count('|') / float(hsp.align_length)))) / float(blast_record.query_letters + alignment.length)

        if equation >= threshold: # cut-off to filter HSPs with longer lengths
          temp_list.append((string.split(alignment.hit_id, '|')[2], alignment.length))
            
          
          #percent_id = (hsp.match.count('|') / float(hsp.align_length)) * 100
          #tuple_values = (alignment.length, blast_record.query, blast_record.query_letters, percent_id, hsp.expect, hsp.score, hsp.bits, hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end, hsp.align_length, hsp.gaps, hsp.query, hsp.match, hsp.sbjct)

          #print alignment.hit_id + ' ' + alignment.hit_def + ' ' + alignment.length
          #print str(alignment.hit_id) + '\t' + str(alignment.length) + '\t' + str(hsp.align_length)
          '''
          output = [blast_record.query[:15],
                    str(blast_record.query_letters),
                    alignment.hit_def[:15],
                    str(alignment.length),
                    str(hsp.align_length),
                    '%.2f' % (percent_id),
                    '%.1f' % (equation)]

          print ' \t'.join(output)
          '''
    if len(temp_list) == 1:
      final_list.append(temp_list)
    elif len(temp_list) > 1:
      sets_list.append(temp_list)
    
    #print "\n"
    counter()
    #if count_seq == 2:
    #  break;

  blast_records.close()
  result_handle.close()
  

while len(sets_list):
  popped_item = sets_list.pop()
  #print popped_item
  recycled = False
  
  for index in range(0,len(sets_list)):
    #print set(popped_item).isdisjoint(set(item))
    if not set(popped_item).isdisjoint(set(sets_list[index])):
      sets_list.append(list(set.union(set(popped_item),set(sets_list[index]))))
      del sets_list[index]
      recycled = True
      break
    
  if not recycled:
    final_list.append(popped_item)
  #print ','.join(sets_list.pop())

final_list.reverse()
group_num = 0

for group in final_list:
  group_num += 1
  for item in group:
    print "Adding gnl|ARDB|%s to database" % (item[0])
    item_to_db((group_name + ' ' + str(group_num), item[0], str(item[1])))


## End Timer for total BLAST
endTime = datetime.datetime.now()
print 'End time: ', endTime.strftime('%H:%M:%S')
print 'Time elapsed:', endTime - startTime


