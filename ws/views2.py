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
import traceback

from etoxwsapi.v2 import schema
from etoxwsapi.v2.wsbase import WebserviceImplementationBase

PARTNER_ID   = 'FIMIM'
PARTNER_WEB  = 'http://phi.imim.es'
ADMIN_NAME   = 'Manuel Pastor'
ADMIN_EMAIL  = 'manuel.pastor@upf.edu'

BASEDIR = '/home/modeler/soft/eTOXlab/src/'

class WS2(WebserviceImplementationBase):

    def __init__(self):
        self.my_tags = dict()
        self.my_type = dict()
        self.my_mver = dict()
        self.my_models = list()
        
        calculation_info = schema.get('calculation_info')

        for item in os.listdir (BASEDIR):
            if os.path.isdir(BASEDIR+item):

                internaldir = os.listdir (BASEDIR+item)
                if not 'version0001' in internaldir:
                    continue
                
                mlabel = None
                try:
                    f = open (BASEDIR+item+'/service-label.txt')
                    mlabel = f.readline()[:-1]
                    mtype  = f.readline()[:-1]
                    f.close()
                except:
                    continue

                # model properties common to all model versions for this endpoint
                rtype = schema.get("result_endpoint").schema
                if mtype == 'qualitative':
                    rtype['properties']['value'] = { "enum": ["positive", "negative"]}

                # only for exposed models (eTOXlab version > 0.91)
                nmodels=0
                try:
                    f = open (BASEDIR+item+'/service-version.txt')
                except:
                    continue
                
                for line in f:
                    if line[-1]=='\n': line = line[:-1]
                    
                    line_versions=line.split('\t')

                    ## old syntax ( eTOXlab < 0.95 ) for backwards compatibility only
                    if len(line_versions) == 1:          
                        try:
                            mver = int(line_versions[0]) # internal (eTOXlab) model version
                            ever = 1                     # external (eTOXsys) model version
                        except:
                            break
                        
                        if not os.path.isdir (BASEDIR+item+'/version%0.4d'%(mver)):
                            break
                        
                        mid = 'eTOXvault ID '+ mlabel + ' ' + PARTNER_ID

                        ## new API with version support
                        try:
                            new_model = calculation_info.create_object(id=mlabel, category="ENDPOINT", version=str(ever), external_id = mid)
                        except:
                            new_model = calculation_info.create_object(id=mlabel, category="ENDPOINT", external_id = mid)
                            
                        new_model ['return_type_spec'] = rtype
                        
                        self.my_models.append(new_model)
                        self.my_mver [mlabel,str(ever)] = mver
                        nmodels+=1
                        break
                    ###################################################################

                    try:
                        mver = int (line_versions[0])    # internal (dir tree    ) model version
                        ever = int (line_versions[1])    # external (user defined) model version
                    except:
                        continue
                    
                    if ever==0:                          # this model has not been exposed
                        continue
                    
                    if not os.path.isdir (BASEDIR+item+'/version%0.4d'%mver):    # make sure this version exists
                        continue
                    
                    mid = 'eTOXvault ID '+ mlabel + ' ' + PARTNER_ID
                    

                    ## new API with version support
                    try:
                        new_model = calculation_info.create_object(id=mlabel, category="ENDPOINT", version=str(ever), external_id = mid)
                    except:
                        new_model = calculation_info.create_object(id=mlabel, category="ENDPOINT", external_id = mid)
                    
                    new_model ['return_type_spec'] = rtype
                    
                    self.my_models.append(new_model)
                    self.my_mver [mlabel,str(ever)] = mver
                    nmodels+=1    
                f.close()
            
                if nmodels==0:   # do not add item and mtype unless there is any published model
                    continue
                
                # my_program lists the endpoint tags (e.g. CACO2, ABCB1) 
                self.my_tags [mlabel] = item
                self.my_type [mlabel] = mtype

                
    def info_impl(self):
        ws_info = schema.get('ws_info')
        data = { "provider": PARTNER_ID,
                 "homepage": PARTNER_WEB,
                 "admin": ADMIN_NAME,
                 "admin_email": ADMIN_EMAIL
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

        ## VERSIONING
        try:
            imver = self.my_mver[calc_info ['id'], calc_info['version']]
        except:
            imver = -1
        
        tdir  = tempfile.mkdtemp(dir='/var/tmp')
        tfile = tdir + '/input_file.sdf'

        with open(tfile, "wb") as fp:
            fp.write(sdf_file)

        logging.info("calculation for %s"%(calc_info['id']))

        os.chdir(tdir)
        
##        p = subprocess.Popen(['/usr/bin/python', BASEDIR+'predict.py','-e',itag,'-b']
##                              ,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        ## VERSIONING
        p = subprocess.Popen(['/usr/bin/python', BASEDIR+'predict.py','-e',itag,'-v',imver, '-c']
                              ,stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        jobobserver.report_progress(0) # indicate that calculation has been started
        while True:
            retcode = p.poll() # returns None while subprocess is running
            line = p.stdout.readline()
            if (retcode is not None):
                break
            else:              # computation progress
                try:
                    jobobserver.log_info(line)
                    if line.startswith('completed:'):
                        jobobserver.report_progress(int(line.split()[1]))
                        logging.info(line)
                except Exception, e:
                    jobobserver.log_warn("Failed to parse log output: %s (%s)"%(line, e))

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
                    result = dict()
                    result['cmp_id'] = i

                    result['success'] = False
                    result['message'] = 'unknown error'
                    result['AD'] = { "value": "", "success": False, "message": 'unknown error' }
                    result['RI'] = { "value": "", "success": False, "message": 'unknown error' }
                    try:
                        r = line.strip().split('\t')
                        
                        if r[0]=='1':        # value has been computed with success
                            result['success'] = True
                            result['message'] = ''
                            
                            if itype == 'quantitative':
                                result['value'] = float(r[1])
                            else:
                                result['value'] = r[1]
                        else:
                            result['success'] = False
                            result['message'] = r[1]

                        if r[2]=='1':        # AD OK
                            result['AD'] = { "value": float(r[3]), "success": True, "message": "" }
                        else:
                            result['AD'] = { "success": False, "message": r[3] }

                        if r[4]=='1':        # RI OK
                            result['RI'] = { "value": float(r[5]), "success": True, "message": "" }
                        else:
                            result['RI'] = { "success": False, "message": r[5] }
                            
                    except Exception, e:
                        result['message'] = json.dumps(traceback.format_exc().splitlines())


                    jobobserver.report_result(i, json.dumps(result))

        os.chdir("..")
        try:
            shutil.rmtree(tdir)
        except Exception, e:
            jobobserver.log_warn("Could not delete temporary dir: %s (%s)"%(tdir, e))

