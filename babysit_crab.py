#!/usr/bin/env python

import sys
import os
import subprocess
from subprocess import *
import sqlite3
from time import sleep



def docrabstatus(dset):
    print 'Crab status on dataset '+dset
    os.system('crab -c %s -status' % str(dset))

def docrabget(dset):
    print 'Crab get on dataset '+dset
    os.system('crab -c %s -get' % str(dset))

def fullupdate(thisdset):
    print 'Processing full update of dataset '+thisdset
    docrabstatus(thisdset)
    docrabget(thisdset)
    docrabstatus(thisdset)

def checkandresubmit():
    for dset in dsets:
        fullupdate(dset)
        ntot=getnjobs(dset)
        jobstoresub=[]
        jobstosub=[]
        nfinished=0
        for j in range(ntot):
            n=j+1
            nfinished=nfinished+isfinished(dset,n)
            action=decideaction(dset,n)
            if (action=='OK'):
                continue
            elif (action=='Submit'):
                jobstosub.append(str(n))
            elif (action=='Resubmit'):
                jobstoresub.append(str(n))
        if (nfinished==ntot):
            dsets.remove(dset)
            continue
        if (len(jobstoresub)>0):
            os.system('echo EXECUTING crab -c %s -resubmit %s' % (dset,str(',').join(jobstoresub)))
            os.system('crab -c %s -resubmit %s' % (dset,str(',').join(jobstoresub)))
        if (len(jobstosub)>0):
            os.system('echo EXECUTING crab -c %s -submit %s' % (dset,str(',').join(jobstoresub)))
            os.system('crab -c %s -submit %s' % (dset,str(',').join(jobstoresub)))            
        sleep(20)
        docrabstatus(dset)
        jobstoresub=[]
        jobstosub=[]
        tovetoedsites={}
        for j in range(ntot):
            n=j+1
            action=decideaction(dset,n)
            if (action=='OK'):
                continue
            elif (action=='Resubmit'):
                jobstoresub.append(str(n))
            elif (action=='ResubmitDifferent'):
                temp=getdestsite(dset,n).split('.')
                if (len(temp)<3):
                    continue
                site=temp[len(temp)-2]+'.'+temp[len(temp)-1]
                if site not in tovetoedsites.keys():
                    tovetoedsites[site]=[]
                tovetoedsites[site].append(str(n))
            if (len(jobstoresub)>0):
                os.system('echo EXECUTING crab -c %s -forceResubmit %s' % (dset,str(',').join(jobstoresub)))
                os.system('crab -c %s -forceResubmit %s' % (dset,str(',').join(jobstoresub)))
        for site in tovetoedsites.keys():
            if (len(tovetoedsites[site])>0):
                os.system('echo EXECUTING crab -c %s -GRID.se_black_list=T3,%s -forceResubmit %s' % (dset,site,str(',').join(tovetoedsites[site])))
                os.system('crab -c %s -GRID.se_black_list=T3,%s -forceResubmit %s' % (dset,site,str(',').join(tovetoedsites[site])))
    for i in dsets:
        docrabstatus(i)

def isfinished(dset,n):
    if (getwrappercode(dset,n)==0):
        return True
    else:
        return False

def decideaction(dset,n):
    st = getstatus(dset,n)
    code = getwrappercode(dset,n)
    out = ''
    if (code==0):
        out= 'OK'
    elif (st=='Running' or st=='Submitted'):
        out= 'OK'
    elif (st=='Created'):
        out= 'Submit'
    elif (st=='Aborted' or st=='Cancelled' or 60000<=code<70000):
        out= 'Resubmit'
    elif (7000<=code<11000 or 50000<=code<60000):
        out= 'ResubmitDifferent'
    else:
        out= 'Printout'
    if (out!='OK'):
        print 'Found job %s nr. %s: %s %s' % (str(dset),str(n),str(st),str(code))
    return out

def babysit():
    while (len(dsets)>0):
        checkandresubmit()
        sleep(30*60)
    sys.exit(0)

def getnjobs(dset):
    curs[dset].execute('select count(*) from bl_job')
    return int(curs[dset].fetchall()[0][0])

def getwrappercode(dset,job):
    curs[dset].execute('select wrapper_return_code from bl_runningjob where job_id=? order by submission desc limit 1',(int(job),))
    out = curs[dset].fetchall()[0][0]
    if (out==None):
        return -1
    else:
        return int(out)

def getexecode(dset,job):
    curs[dset].execute('select application_return_code from bl_runningjob where job_id=? order by submission desc limit 1',(int(job),))
    return int(curs[dset].fetchall()[0][0])

def getstatus(dset,job):
    curs[dset].execute('select status_scheduler from bl_runningjob where job_id=? order by submission desc limit 1',(int(job),))
    return str(curs[dset].fetchall()[0][0])

def getdestsite(dset,job):
    curs[dset].execute('select destination from bl_runningjob where job_id=? order by submission desc limit 1',(int(job),))
    out = str(curs[dset].fetchall()[0][0])
    return str(out.partition(":")[0])

def getdetailsfail(dset,job):
    curs[dset].execute('select wrapper_return_code,application_return_code,destination from bl_runningjob where job_id=? order by submission desc limit 1',(int(job),))
    return curs[dset].fetchall()[0]




if (len(sys.argv)<2):
    print 'Wrong usage'
    sys.exit(1)

dir = str(sys.argv[1])
print 'Working on directory '+dir

os.chdir(dir)
dsets = sorted(os.listdir(dir))
dsets2=[]
for i in dsets:
    if (str(i).find('multicrab')>=0):
        continue
    if (str(i).find('history')>=0):
        continue
    dsets2.append(i)
dsets=dsets2
conn={}
curs={}
for i in dsets:
    conn[i]=sqlite3.connect(str(dir)+'/'+str(i)+'/share/crabDB')
    curs[i]=conn[i].cursor()

babysit()










    
#sqlite> .schema
#CREATE TABLE bl_job
#  (
#          id INTEGER PRIMARY KEY AUTOINCREMENT,
#              task_id INT NOT NULL,
#              job_id INT NOT NULL,
#              wmbsJob_id INT ,
#              name VARCHAR(255),
#              executable TEXT,
#              events INT,
#              arguments TEXT,
#              stdin TEXT,
#              stdout TEXT,
#              stderr TEXT,
#              input_files TEXT,
#              output_files TEXT,
#              dls_destination TEXT,
#              submission_number INT default 0,
#              closed CHAR default "N",
#              UNIQUE(job_id, task_id),
#              FOREIGN KEY(task_id) references bl_task(id) ON DELETE CASCADE
#            );
#  CREATE TABLE bl_runningjob
#    (
#            id INTEGER PRIMARY KEY AUTOINCREMENT,
#                job_id INT NOT NULL,
#                task_id INT NOT NULL,
#                submission INT NOT NULL,
#                state VARCHAR(255),
#                scheduler TEXT,
#                service TEXT,
#                sched_attr TEXT,
#                scheduler_id VARCHAR(255),
#                scheduler_parent_id VARCHAR(255),
#                status_scheduler VARCHAR(255),
#                status VARCHAR(255),
#                status_reason TEXT,
#                destination TEXT,
#                creation_timestamp TIMESTAMP,
#                lb_timestamp TIMESTAMP,
#                submission_time TIMESTAMP,
#                scheduled_at_site TIMESTAMP,
#                start_time TIMESTAMP,
#                stop_time TIMESTAMP,
#                stageout_time TIMESTAMP,
#                getoutput_time TIMESTAMP,
#                output_request_time TIMESTAMP,
#                output_enqueue_time TIMESTAMP,
#                getoutput_retry INT,
#                output_dir TEXT,
#                storage TEXT,
#                lfn TEXT,
#                application_return_code INT,
#                wrapper_return_code INT,
#                process_status TEXT default 'created',
#                closed CHAR default "N",
#                UNIQUE(submission, job_id, task_id),
#                FOREIGN KEY(job_id) references bl_job(id) ON DELETE CASCADE,
#                FOREIGN KEY(task_id) references bl_task(id) ON DELETE CASCADE
#              );
#    CREATE TABLE bl_task
#      (
#              id INTEGER PRIMARY KEY AUTOINCREMENT,
#                  name VARCHAR(255),
#                  dataset VARCHAR(255),
#                  start_dir TEXT,
#                  output_dir TEXT,
#                  global_sanbox TEXT,
#                  cfg_name TEXT,
#                  server_name TEXT,
#                  job_type TEXT,
#                  total_events INT,
#                  user_proxy TEXT,
#                  outfile_basename TEXT,
#                  common_requirements TEXT,
#                  unique(name)
#                );
      
