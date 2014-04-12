#-*- coding: utf-8 -*-

import json
import sys
import shutil
import subprocess
import os
import tempfile
import logging

from etoxwsapi.v2 import schema
from etoxwsapi.v2.wsbase import WebserviceImplementationBase

BASEDIR = '/home/modeler/soft/eTOXlab/src/'

class WS2(WebserviceImplementationBase):

    def __init__(self):
        self.my_tags = dict()
        self.my_type = dict()
        self.my_models = list()
        
        calculation_info = schema.get('calculation_info')

        mcount = 0
        for item in os.listdir (BASEDIR):
            if os.path.isdir(BASEDIR+item):
                mlabel = None
                try:
                    f = open (BASEDIR+item+'/service-label.txt')
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
                self.my_type [mlabel] = mtype
                
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

        itag  = self.my_tags[calc_info ['id']]
        itype = self.my_type[calc_info ['id']]
        
        tdir  = tempfile.mkdtemp(dir='/var/tmp')
        tfile = tdir + '/input_file.sdf'

        with open(tfile, "wb") as fp:
            fp.write(sdf_file)

        logging.info("calculation for %s"%(calc_info['id']))

        os.chdir(tdir)
        
        p = subprocess.Popen(['/usr/bin/python', BASEDIR+'predict.py','-e',itag,'-b']
                              ,stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        while True:
            retcode = p.poll() #returns None while subprocess is running
            line = p.stdout.readline()
            if (retcode is not None):
                break
            else:
                if line.startswith('completed:'):
                    jobobserver.report_progress(int(line.split()[1]))
                    logging.info(line)

        jobobserver.report_status(retcode, p.stderr.read())
        if retcode == 0:
            
            with open('result.txt') as fp:
                for i, line in enumerate(fp):
                    r = line.strip().split('\t')
        
                    result = dict()
                    result['cmp_id'] = i

                    if len (r) != 6:
                        result['success'] = False
                        result['message'] = 'unknown error'
                        result['AD'] = { "value": "", "success": False, "message": 'unknown error' }
                        result['RI'] = { "value": "", "success": False, "message": 'unknown error' }
                    
                    if r[0]=='1':
                        if itype == 'quantitative':
                            result['value'] = float(r[1])
                        else:
                            result['value'] = r[1]
                        result['success'] = True
                    else:
                        result['success'] = False
                        result['message'] = r[1]
                        
                    if r[2]=='1':
                        result['AD'] = { "value": float(r[3]), "success": True, "message": "" }
                    else:
                        result['AD'] = { "success": False, "message": r[3] }
                    
                    if r[4]=='1':
                        result['RI'] = { "value": float(r[5]), "success": True, "message": "" }
                    else:
                        result['RI'] = { "success": False, "message": r[5] }
                        
                    jobobserver.report_result(i, json.dumps(result))

##        os.chdir("..")
##        shutil.rmtree(tdir)


