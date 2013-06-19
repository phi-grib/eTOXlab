#-*- coding utf-8 -*-

from cStringIO import StringIO
from django.http import HttpResponse
from django.utils import simplejson
from django.core.context_processors import csrf
from django.views.decorators.csrf import csrf_exempt

import os
import shutil
import subprocess
import time
import datetime
import tempfile
import random

MODEL_DIPL1 = "/Toxicity/Organ Toxicity/Phospholipidosis/DIPL/1"
MODEL_DIPL2 = "/Toxicity/Organ Toxicity/Phospholipidosis/DIPL/2"
MODEL_CACO2 = "/ADME/Absorption/Gastro Intestinal/CACO2 permeability/1"

BASEDIR = '/srv/www/webapps/etoxws_FIMIM/src/etoxws/models/'
CACODIR = '/srv/www/webapps/etoxws_FIMIM/src/etoxws/CACO2/'

def info (request):
    jsondata = { "provider": "FIMIM",
	"homepage": "http://phi.imim.es",
	"admin": "Manuel Pastor",
	"admin-email": "manuel.pastor@upf.edu"
    }
    json = simplejson.dumps(jsondata)
    return HttpResponse(json, mimetype="application/json")

def available_services (request):
    jsondata = { "predictions": [ MODEL_DIPL1, MODEL_DIPL2, MODEL_CACO2 ] } 
    json = simplejson.dumps(jsondata)
    return HttpResponse(json, mimetype="application/json")

def handle_uploaded_file(ufil, format):

    tdir = tempfile.mkdtemp(dir=BASEDIR+'temp')
    tfil = tdir + '/input_file'
    
    if 'SMILES' in format:
        tfil += '.smi'
    elif 'SDF' in format:
        tfil += '.sdf'

    fou = open (tfil, 'wb')
    for chunk in ufil.chunks():
        fou.write(chunk)
    fou.close()
    
    return tdir, tfil

@csrf_exempt
def calculate (request):
    
    results = list()
    status_code = 200;
    msg = ""
        
    try:
        
        if request.method != 'POST':
            raise Exception, "Only HTTP POST is supported"

        indata = request.POST

        property = indata['property']
        if not property in ( MODEL_DIPL1, MODEL_DIPL2, MODEL_CACO2 ):
            raise Exception, "Unknown property: '%s'"%(property)
        
        format = indata['format']
        if not format in ("SDF", "SMILES"):
            raise Exception, "Unknown file format: '%s'"%(format)

        target_dir, target_file = handle_uploaded_file(request.FILES['uploadfile'], format)
       
        tmp_project = 'project%s' % (str(random.randrange(0,100000,1)))

        os.chdir(target_dir)

        if property == MODEL_DIPL1:
            call = ['/usr/bin/python', BASEDIR+'dipl1.py','-p', tmp_project]
        elif property == MODEL_DIPL2:
            call = ['/usr/bin/python', BASEDIR+'dipl2.py']
        elif property == MODEL_CACO2:
##            call = ['/usr/bin/python', BASEDIR+'caco2.py','-p', tmp_project]
            os.environ['ETAM_OPT'] = '/opt/'
            os.environ['ETAM_HOME'] = '/srv/www/webapps/etoxws_FIMIM/src/etoxws/CACO2'
            call = ['/usr/bin/python', CACODIR+'predict.py','-f', './input_file.sdf']

        execution = subprocess.Popen(call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        result, errors = execution.communicate()

        if execution.returncode != 0:
            raise Exception, errors
        else:
            # creating the LOG file
            now = datetime.datetime.now()
            client_address = request.META['REMOTE_ADDR']

            logfile = open (BASEDIR+'LOG', 'a+')
            # print the client IP, filetype and prop request in LOG file
            logfile.write('\nREQUEST:\t%s\t%s\t %s\t %s\nINPUT_FILE:\n' % (now, client_address, property, format))
            #inputfile = open(target_file, 'rb')
            #infile = inputfile.read()
            #inputfile.close()
            #logfile.write(infile)
            logfile.write('\nRESULTS:\n')
            logfile.write(result)
            logfile.close()

            # printing the results
            results_array = result.split()

            for val in results_array:
                if 'NA' in str(val):
                    if ( property == MODEL_CACO2 ):
                        val = 0.0
                    else:
                        val = ""
                    stat = 1
                    msg = "undefined error"
                else:
                    if ( property == MODEL_CACO2 ):
                        val = float (val)
                    stat = 0
                    msg = ""

                results.append((val,stat,msg))

            # uncomment to delete the temp files and keep commented for debugging
            os.chdir("..")
            shutil.rmtree(target_dir)

    except Exception, e:
            status_code = 500
            errfile = open (BASEDIR+'ERR', 'a+')
            errfile.write('\nERROR: '+str(e)+'\n')
            errfile.close()

    jsondata = dict()
    jsondata['property'] = property
    jsondata['results'] = results
    jsondata['msg'] = msg

    json = simplejson.dumps(jsondata)
    return HttpResponse(json, mimetype="application/json", status=status_code)
