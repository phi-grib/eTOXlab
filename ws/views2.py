#-*- coding: utf-8 -*-

##    Description    WS2 class. Class connecting eTOXlab with eTOXsys (API 2.0)
##                   
##    Authors:       Manuel Pastor (manuel.pastor@upf.edu)
##                   API 2.0 template provided by Thomas Kleinodel (MN)
##
##    Copyright 2014 Manuel Pastor
##
##    This file is part of eTOXlab.
##
##    eTOXlab is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation version 3.
##
##    eTOXlab is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with eTOXlab.  If not, see <http://www.gnu.org/licenses/>.

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


    ################################################################################
    ##
    ## implementation of prediction
    ##
    ## this method calls the eTOXlab predict method with:
    ##
    ##      -e itag (where itag identifies the endpoint, e.g. CACO2)
    ##      -b      (use API 2.0 formated output)
    ##
    ## the input molecule is copied to /var/tmp/TEMPDIR/input_file.sdf
    ##
    ## output to STDOUT indicates progress
    ## output to STDERR for errors
    ##
    ## prediction results are dumped to /var/tmp/TEMPDIR/resuts.txt
    ##
    ## TO DO: more flexible input/output by
    ##       1. dump extra calc info to /var/tmp/TEMPDIR/meta_in
    ##       2. other prediction results (e.g. struct) to /var/tmp/TEMPDIR/meta_out
    ################################################################################
    def calculate_impl(self, jobobserver, calc_info, sdf_file):   

        itag  = self.my_tags[calc_info ['id']]      # -e tag for predict.py
        itype = self.my_type[calc_info ['id']]      # quant/qualit endpoint
        
        tdir  = tempfile.mkdtemp(dir='/var/tmp')
        tfile = tdir + '/input_file.sdf'

        with open(tfile, "wb") as fp:
            fp.write(sdf_file)

        logging.info("calculation for %s"%(calc_info['id']))

        os.chdir(tdir)
        
        p = subprocess.Popen(['/usr/bin/python', BASEDIR+'predict.py','-e',itag,'-b']
                              ,stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        while True:
            retcode = p.poll() # returns None while subprocess is running
            line = p.stdout.readline()
            if (retcode is not None):
                break
            else:              # computation progress
                if line.startswith('completed:'):
                    jobobserver.report_progress(int(line.split()[1]))
                    logging.info(line)

        jobobserver.report_status(retcode, p.stderr.read())
        if retcode == 0:

            #################################################################
            ## results.txt is formated as: (one line per compound)
            ## 1  0.234  1  0  1  0.023
            ##
            ##  1       prediction success, True if predicted value obtained
            ##  0.234   predicted value
            ##  1       AD success, True if AD obtained
            ##  0       number of ADAN criteria broken: (0-6) for quant, 0-3 for qualit
            ##  1       RI success, True if RI obtained
            ##  0.023   95% CI for the predicted value. Only for quant
            ##
            #################################################################
            with open('result.txt') as fp:
                for i, line in enumerate(fp):
                    r = line.strip().split('\t')
        
                    result = dict()
                    result['cmp_id'] = i

                    if len (r) != 6:     # catastrophic error, unformated line
                        result['success'] = False
                        result['message'] = 'unknown error'
                        result['AD'] = { "value": "", "success": False, "message": 'unknown error' }
                        result['RI'] = { "value": "", "success": False, "message": 'unknown error' }
                    
                    if r[0]=='1':        # formated line
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


