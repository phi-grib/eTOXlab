#-*- coding: utf-8 -*-

import json
import sys
import shutil
import subprocess
import os
import re
import tempfile
import logging

from etoxwsapi.v2 import schema
from etoxwsapi.v2.wsbase import WebserviceImplementationBase

BASEDIR = '/srv/www/webapps/etoxws/src/etoxws/models/'

class WS2(WebserviceImplementationBase):

    def __init__(self):
        self.my_tags = dict()
        self.my_models = list()
        
        calculation_info = schema.get('calculation_info')

        mdir = BASEDIR + 'src/'

        mcount = 0
        for item in os.listdir (mdir):
            if os.path.isdir(mdir+item):
                mlabel = None
                try:
                    f = open (mdir+item+'/service-label.txt')
                    mlabel = f.readline()[:-1]
                    mtype  = f.readline()[:-1]
                    f.close()
                except:
                    continue

                # only for published models
                mcount += 1
                mid = "eTOXvault ID%d"%mcount
                rtype = schema.get("result_endpoint").schema
                if mtype == 'qualitative':
                    rtype['properties']['value'] = { "enum": ["positive", "negative"]}
                    
                new_model = calculation_info.create_object(id=mlabel, category="ENDPOINT", external_id = mid)
                new_model ['return_type_spec'] = rtype
                self.my_models.append(new_model)

                # my_program lists the endpoint tags (e.g. CACO2, ABCB1) 
                self.my_tags [mlabel] = item
                
    def info_impl(self):
        ws_info = schema.get('ws_info')
        data = { "provider": "FIMIM",
                 "homepage": "http://phi.imim.es",
                 "admin": "Manuel Pastor",
                 "admin_email": "manuel.pastor@upf.edu",
                 }
        ws = ws_info.create_object(**data)
        return ws.to_json()

    def dir_impl(self):
        return json.dumps(self.my_models)

    def calculate_impl(self, jobobserver, calc_info, sdf_file):   

        itag = self.my_tags[calc_info ['id']]
        
        tdir  = tempfile.mkdtemp(dir=BASEDIR+'temp')
        tfile = tdir + '/input_file.sdf'

        with open(tfile, "wb") as fp:
            fp.write(sdf_file)

        logging.info("calculation for %s"%(calc_info['id']))

        regex = re.compile("\*\*\* RECORD no\.:\s+(\d+)\s+read \*") 

        os.chdir(tdir)
        
        p = subprocess.Popen(['/usr/bin/python', BASEDIR+'src/predict.py','-e',itag,'-b']
                              ,stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        while True:
            retcode = p.poll() #returns None while subprocess is running
            line = p.stdout.readline()
            if (retcode is not None):
                break
            else:
                m = regex.search(line)
                if m:
                    jobobserver.report_progress(int(m.group(1)))

        jobobserver.report_status(retcode, p.stderr.read())
        if retcode == 0:
            
            with open('result.txt') as fp:
                for i, line in enumerate(fp):
                    r = line.strip().split('\t')
        
                    result = dict()
                    result['cmp_id'] = i
                    result['value'] = r[0]
                    result['AD'] = r[1]
                    result['RI'] = r[2]
                    jobobserver.report_result(i, json.dumps(result))

##        os.chdir("..")
##        shutil.rmtree(tdir)


