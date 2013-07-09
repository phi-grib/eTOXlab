#-*- coding utf-8 -*-

from cStringIO import StringIO
from django.http import HttpResponse
from django.utils import simplejson
from django.core.context_processors import csrf
from django.views.decorators.csrf import csrf_exempt

import os
import sys
import shutil
import subprocess
import time
import datetime
import tempfile
import cPickle as pickle

BASEDIR = '/srv/www/webapps/etoxws_FIMIM/src/etoxws/models/'

def getModels():

    services = []
    models = dict()

    mdir = BASEDIR + 'src/'
    
    for item in os.listdir (mdir):
        if os.path.isdir(mdir+item):
            label = None
            try:
                f = open (mdir+item+'/service-label.txt')
                label = f.readline()[:-1]
                f.close()
            except:
                continue
            services.append(label)
            models[label] = item
            
    return (services,models)


def info (request):
    jsondata = { "provider": "FIMIM",
	"homepage": "http://phi.imim.es",
	"admin": "Manuel Pastor",
	"admin-email": "manuel.pastor@upf.edu"
    }
    json = simplejson.dumps(jsondata)
    return HttpResponse(json, mimetype="application/json")

def available_services (request):

    services, models = getModels()
    
    jsondata = { "predictions": services } 
    json = simplejson.dumps(jsondata)
    return HttpResponse(json, mimetype="application/json")

def savefile(ufile, iformat):

    tdir  = tempfile.mkdtemp(dir=BASEDIR+'temp')
    tfile = tdir + '/input_file'
    
    if 'SMILES' in iformat:
        tfile += '.smi'
    elif 'SDF' in iformat:
        tfile += '.sdf'

    fou = open (tfile, 'wb')
    for chunk in ufile.chunks():
        fou.write(chunk)
    fou.close()
    
    return tdir, tfile

@csrf_exempt
def calculate (request):
    
    results = list()   
    msg = ""
              
    if request.method != 'POST':
        raise Exception ('Only HTTP POST is supported')

    services, models = getModels()
    
    indata = request.POST

    iproperty = indata['property']
    
    if not iproperty in services:
        raise Exception ('Unknown property: %s'%(iproperty))
    
    iformat = indata['format']
    if not iformat in ("SDF", "SMILES"):
        raise Exception ('Unknown file format: %s'%(iformat))

    tdir, tfile = savefile(request.FILES['uploadfile'], iformat)      
    
    os.chdir(tdir)
    
    call = ['/usr/bin/python', BASEDIR+'src/predict.py', '-e',  models[iproperty], '-a']

    xresults = 'undefined error'   # in case xresults is not generated
    
    try:
        #retcode = subprocess.call(call,stderr=stderrf)
        retcode = subprocess.call(call)
        #stderrf.close()
        
        if retcode != 0:
            raise Exception ('Error in call: '+str(call))

        pkl = open('./results.pkl', 'rb')
        xresults = pickle.load(pkl)
        pkl.close()

    except Exception as e:
        status_code = 500
        errfile = open (BASEDIR+'ERR', 'a+')
        msg = str(e)
        errfile.write(msg)
        errfile.write('\n')            
        errfile.close()
        xresults = msg
        
    else:
        status_code = 200;
        logfile = open (BASEDIR+'LOG', 'a+')
        logfile.write('\nREQUEST:\t%s\t%s\t %s\t %s\n' %(
            datetime.datetime.now(),
            request.META['REMOTE_ADDR'],
            iproperty,
            iformat))    
        logfile.write('RESULTS:\n')
        for val, stat, msg in xresults:
            logfile.write('val: '+str(val)+' stat: '+str(stat)+' msg: '+msg+'\n')       
        logfile.close()
        
    jsondata = dict()
    jsondata['property'] = iproperty
    jsondata['results'] = xresults
    jsondata['msg'] = msg

    os.chdir("..")
    shutil.rmtree(tdir)
    
    json = simplejson.dumps(jsondata)
    return HttpResponse(json, mimetype="application/json", status=status_code)
