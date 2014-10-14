#!/usr/bin/env python

# -*- coding: utf-8 -*-

##    Description    eTOXsys-like (API 2.0) for testing web service
##                   
##    Authors:       Thomas Kleinoder and Manuel Pastor  
##

import sys
import os
import getopt

import requests
import time
from etoxwsapi.v2 import schema
import json

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

BASE_URL = 'http://localhost/etoxwsapi/v2'


def print_result(ids):
        frmt = '| %-7s | %-15s | %-15s | %-15s |'
        for (job_id,model_id) in ids:
                print "==========================================================================="
                print "Results for model '%s' (job: %s)\n"%(model_id, job_id)
                ret = requests.get('/'.join((BASE_URL, 'jobs', job_id)))
                if ret.status_code == 200:
                        results = json.loads(ret.text)
                        if results['status'] == "JOB_COMPLETED":
                                print frmt%("cmp_id", "value", "AD", "RI")
                                print '-' * len(frmt%('','','',''))
                                for r in results['results']:
                                        cid = r['cmp_id']
                                        if 'value' in r.keys():
                                                val = r['value']
                                        else:
                                                val = 'na'
                                        if r['AD']['success']:
                                                ad = r['AD']['value']
                                        else:
                                                ad = 'na'
                                        if r['RI']['success']:
                                                ri = r['RI']['value']
                                        else:
                                                ri = 'na'

                                        print frmt%(cid,val,ad,ri)
                                        #print frmt%(r['cmp_id'], r['value'], r['AD']['value'], r['RI']['value'])
                        else:
                                print "Job is not (yet) available: %s"%(results['status'])
                else:
                        print "Could not get job results: %s, %s"%(ret.status_code, ret.text)

def submit_jobs(models, filename):
        """
        preparing request for calculation
        """
        calculation_request = schema.get('calculation_request')
        req_obj = calculation_request.create_object()

        calculation_info = schema.get('calculation_info')
        
        req_obj.req_calculations = models

        fname = os.path.join(THIS_DIR, filename)
        with open(fname) as fp:
                req_obj.sdf_file = fp.read()

        req_ret = requests.post(BASE_URL+"/jobs/", data = req_obj.to_json())
        
        if req_ret.status_code == 200:
                job_ids = list()
                for stat in json.loads(req_ret.text):
                        job_ids.append(stat['job_id'])
                        print "======================================================================"
                        print "new job submitted."
                        for t in ("job_id", "status"): #, "msg"):
                                print "%s: %s"%(t, stat[t])
        else:
                raise Exception("Failed to submit jobs (%s): %s"%(req_ret.status_code, req_ret.text))

        return job_ids

def delete_job(job_id):
        ret = requests.delete(BASE_URL + '/jobs/' + job_id)
        print ret.status_code, ret.text

def observing_jobs(job_ids, duration):
        print "observing running jobs for %dsec. (or until one is completed)"%(duration)
        for i in range(duration):
                for job_id in job_ids:
                        response = requests.get(BASE_URL+"/jobs/%s"%(job_id))
                        if response.status_code == 200:
                                stat = json.loads(response.text)
                                print "status for '%s': %s"%(job_id, stat)
                                if stat['status'] != "JOB_RUNNING":
                                        print "Job not running anymore."
                                        return
                        else:
                                print response.status_code
                                print response.text
                time.sleep(1)
        
def process_request(endpoint, filename):
        try:
                ret = requests.get(BASE_URL+"/dir")
                
                if ret.status_code != 200:
                        print ret.text

                models = []
                model_ids = []
                ret = requests.get(BASE_URL+"/dir")

                if endpoint==None:
                        models = [ m for m in json.loads(ret.text)]
                        model_ids = [ model['id'] for model in models ]
                else:         
                        for m in json.loads(ret.text):
                                #print m
                                if endpoint in m['id']:
                                        models.append(m)
                                        model_ids.append(m['id'])

                print "Selected models"
                print models

                #print "All jobs:"
                #ret = requests.get(BASE_URL+"/jobs/")
                #print ret.text

                job_ids = submit_jobs(models, filename)

                observing_jobs(job_ids, 10)

                print_result(zip(job_ids, model_ids) )

                #time.sleep(5)
                #delete_job(job_ids[0])
                
        except Exception, e:
                print "Error occurred."
                print "%s"%(e)
                return 1
        return 0


def main ():

        endpoint = None
        filename = None

        try:
                opts, args = getopt.getopt(sys.argv[1:], 'e:f:')
        except getopt.GetoptError:
                print 'Error. Arguments not recognized'
                sys.exit(1)

        if args:
                print 'Error. Arguments not recognized'
                sys.exit(1)
        
        if len( opts ) > 0:
                for opt, arg in opts:
                        if opt in '-e':
                                endpoint = arg
                        elif opt in '-f':
                                filename = arg              

        if not filename:
                filename = 'input_file.sdf'

        process_request(endpoint, filename)
                
        sys.exit(0)


if __name__ == "__main__":

        main ()
