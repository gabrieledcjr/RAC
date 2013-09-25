#!/usr/bin/python -tt
'''
Created on June 14, 2013

@author: Gabriel de la Cruz

@description: MySQL wrapper created to utilize reusability of functions
when connecting to MySQL
'''

# used for mysql connection
import mysql.connector

# use for mysql connection exception handling
from mysql.connector import errorcode


#*******************************************************************************
# @func_name: connect_db()
# @desc: Establish connection to specified MySQL database
# @return: MySQL Database connection
#*******************************************************************************
def connect_db(db_username, db_password, db_host, db_name):
  
  try:
    db = mysql.connector.connect(user=db_username,password=db_password,
                                 host=db_host,database=db_name)
    
  except mysql.connector.Error as err:
    if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
      print('MySQL_Err: Something is wrong with your username or password')
    elif err.errno == errorcode.ER_BAD_DB_ERROR:
      print('MySQL_Err: Database does not exists')
    else:
      print(err)
      
  else:
    print('MySQL: Successfully established connection to %s ' % (db_name) +
          'database.')
    return db
  
##end connect_db()
          

#*******************************************************************************
# @func_name: select_query()
# @desc: Execute a SELECT SQL statement
# @return: Cursor when SQL is successfully executed
#*******************************************************************************
def select_query(db, isBuffered, sql_statement):
  
  try:
    # setting buffered=True lets you open multiple cursors
    cursor = db.cursor(buffered=isBuffered)
    cursor.execute(sql_statement)
          
  except mysql.connector.Error as err:
    # err.errno == errorcode.ER_NO_SUCH_TABLE:
    # err.errno == errorcode.ER_TABLE_EXISTS_ERROR:
    print(err)
    
  else:
    # print('MySQL: Successfully retrieved MySQL data.') 
    return cursor
  
## end select_query()


#*******************************************************************************
# @func_name: insert_query()
# @desc: Execute an INSERT SQL statement
#*******************************************************************************
def insert_query(db, sql_statement, data):
  
  try:
    # opens a db cursor
    cursor = db.cursor()

    # adds data to database
    cursor.execute(sql_statement, data)

    # saves data to database
    db.commit()
    
  except mysql.connector.Error as err:
    print(err)

  else:
    # closes db cursor
    cursor.close()
  
## end insert_query()


#*******************************************************************************
# @func_name: create_table()
# @desc: Creates table and if table exist in the database. Existing table will
#        be deleted and a new empty table will be created
#*******************************************************************************
def create_table(db, table_name, sql_statement):
  
  try:
    # opens a db cursor
    cursor = db.cursor()
    # creates table
    cursor.execute(sql_statement)
    
  except mysql.connector.Error as err:
    # exception handling, table already exist
    if err.errno == errorcode.ER_TABLE_EXISTS_ERROR:
      try:
        print("MySQL: %s already exists." % (table_name))
        # deletes existing table
        cursor.execute(("DROP TABLE %s " % (table_name)))
        print("MySQL: %s deleted" % (table_name))

        # recreate table
        cursor.execute(sql_statement)
        
      except mysql.connector.Error as err:
        print(err)
        
      else:
        print("MySQL: Successfully Created %s Table" % (table_name))
        
    else:
      print(err)
      
  else:
    print("MySQL: Successfully Created %s Table" % (table_name))

  # closes db cursor
  cursor.close()

##end create_table()
