#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    eTOXlab simple GUI
##                   
##    Authors:       Ines Martinez and Manuel Pastor (manuel.pastor@upf.edu) 
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

from Tkinter import *  # Importing the Tkinter (tool box) library 
import Tkconstants

import tkMessageBox
import tkFileDialog

import commands
import os
import subprocess
import shutil
import Queue

from threading import Thread

from utils import wkd
#from fabric.colors import white


class buildWorker: 

    def __init__(self, seeds, queue):
        self.seeds = seeds
        self.q = queue

    def rebuild(self):
      
        name    = self.seeds[0]
        version = self.seeds[1]
        filebut = self.seeds[2]

        # Save all in a specific workspace folder
        tempdir='/home/modeler/workspace/'+name         
                  
        try:
            if not os.path.isdir(tempdir): os.mkdir (tempdir)
        except:
            return (False,'unable to create directory '+tempdir)     
        
        try:
            shutil.copy(filebut,tempdir+'/training.sdf')
        except:
            return (False,'unable to copy training at '+tempdir)     
    
        os.chdir(tempdir)
        
        mycommand = [wkd+'/build.py','-e',name,'-f','training.sdf','-v',version]

        try:
            retcode = subprocess.call(mycommand)
        except:
            self.q.put ('building failed')

        self.q.put('building completed')


class etoxlab:

    def __init__(self, master):
        self.datos  = []
        self.seeds  = []
        self.q      = Queue.Queue()
        self.master = master
        self.myfont = 'Courier New'

        f10 = (self.myfont,10)
        f9  = (self.myfont,9)
        
        # GENERAL PANED - m0 windows that contains m1 and m2 paned
        m0 = PanedWindow()
        m0.pack(fill=BOTH, expand=1)
        
        # MODEL DETAILS - Paned m1
        m1 = PanedWindow(m0,orient=VERTICAL)        
        labelList = Label(m1, text='#endpoint   tag', anchor=W, justify=LEFT,font=f10)
        self.listbox = Listbox(root,width=70, height=20,exportselection=False,font=f10)
        m1.add(labelList)  
        m1.add(self.chargeData())
        
        # VERSION DETAILS - Paned m2              
        m2 = PanedWindow(m0, orient=VERTICAL)
        labelListv = Label(root, text="#    MD            mod    mol", anchor=W, justify=LEFT,font=f10)
        self.listboxversion = Listbox(root,width=70, font=f9)        
        m2.add(labelListv)
        m2.add(self.listboxversion)

        self.listbox.bind('<<ListboxSelect>>', self.chargeVersion)
        
        # Buttons
        self.button1 = Button(m2, text = 'model details', command = self.seeDetails, font=f10)
        self.button2 = Button(m2, text = 'rebuild'      , command = self.build     , font=f10)         
        self.button3 = Button(m2, text = 'publish'      , command = self.publish   , font=f10)
        self.button4 = Button(m2, text = 'view'         , command = self.view      , font=f10)
        self.button5 = Button(m2, text = 'Quit'         , command = quit           , font=f10)
        
        m2.add(self.button1)
        m2.add(self.button2)
        m2.add(self.button3)  
        m2.add(self.button4)
        m2.add(self.button5)

        m0.add(m1)
        m0.add(m2)    

        self.chargeVersion(0)

        self.periodicCall()

    def chargeData(self):
        #Read a file and add all the data to a file  
        output=commands.getoutput(wkd+'/manage.py --info=short')
            
        output=output.split("\n")         
       
        for line in output[1:-1]:
            if line.find('MD')==-1:
                if line.find('-')==-1: 
                    self.datos.append(line) 
                    line=line.split("[")                
                    self.listbox.insert(END, '%-12s'%(line[0])+ '%-30s'%(line[1][0:-1]))
        
        #First item selected    
        self.listbox.selection_set( first = 0 )

        #self.chargeVersion(0)
         
        return self.listbox            
    
    def chargeVersion(self,e):
        
        self.listboxversion.delete(0, END)
         
        d=self.datos[self.whichSelected(self.listbox)]    
        name= d.split()[0];
        
        # Obtain information about the model      
        output=commands.getoutput(wkd+'/manage.py -e '+name+' --info=short') 
        
        output= output.split("\n")    
       
        for line in output[2:-1]:
            l=line.strip().split()
            
            if (len(l) == 10):
    
                line =  '%-6s'%l[0] \
                        +'%-16s'%l[1].replace('MD:','') \
                        +'%-8s'%l[2].replace('mod:','') \
                        +'%-6s'%l[5] \
                        +'%-11s'%l[7].replace(':','R2:') \
                        +'%-11s'%l[9].replace(':','Q2:')
            elif (len(l) == 9 ):
                
                line =  '%-6s'%l[0] \
                        +'%-16s'%l[1].replace('MD:','') \
                        +'%-8s'%l[2].replace('mod:','') \
                        +'%-6s'%l[5] \
                        +'%-11s'%l[6] \
                        +'%-11s'%l[7] \
                        +'%-11s'%l[8]
            else:
    
                line = 'not recognized'
                
            self.listboxversion.insert(END, line)    
        
        self.listboxversion.selection_set( first = 0 )
       
        return self.listboxversion
    
    def whichSelected(self,l):
        return int(l.curselection()[0])
    
    def publish(self):   
        # Choose the model
        select = self.datos[self.whichSelected(self.listbox)]  
        name = select.split()[0];
        
        #Publish the model
        subprocess.call(wkd+'/manage.py --publish -e '+name, stdout=subprocess.PIPE, shell=True)  
    
        self.chargeVersion(0)

    def seeDetails(self):
        # Capture the name of the selected model
        d=self.datos[self.whichSelected(self.listbox)]
        name= d.split()[0];
        
        version = self.whichSelected(self.listboxversion).__str__()    
    
        # Obtain information about the model
        output=commands.getoutput(wkd+'/manage.py -e '+name+' -v' +version  +' --info=long')
    
        outputlist = output.split('\n')
        outputlist = outputlist [1:-3]
    
        output = ''
        for l in outputlist: output+= l+'\n'
        
        # Show the collected information in a new window (winDetails)
        winDetails = Tk()
        winDetails.resizable(0.5,0.5)
        win1 = PanedWindow()
        winDetails.wm_title(name+' ver '+version)
        
        lab= Label(winDetails, text = output,font=(myfont, 10),justify=LEFT)
        lab.pack()
        winDetails.mainloop()
    
    def openFile(self):
        filename = tkFileDialog.askopenfilename(parent=root,
                                                initialdir='/home/modeler/workspace',
                                                title="Select a new training set to rebuild the model",
                                                filetypes=[('image files', '.sdf')]) 
        return filename


    def view(self):   
        # Choose the model
        select = self.datos[self.whichSelected(self.listbox)]  
        name = select.split()[0]
    
        # Choose the version
        version = self.whichSelected(self.listboxversion).__str__()
        
        command = wkd+'/view.py -e '+name+' -v '+version
        subprocess.Popen(command, shell=True)


    def periodicCall(self):

        while self.q.qsize():
            try:
                msg = self.q.get(0)
                if 'building' in msg:
                    self.button2.configure(state='normal')
                    self.chargeVersion(0)
                
                tkMessageBox.showinfo("Info Message", msg)

            except Queue.Empty:
                pass

        self.master.after(500, self.periodicCall) # re-call after 500ms
    
##############################################
### Functions to build model in new thread ###
##############################################

    def buildJob(self):       
        job = buildWorker(self.seeds, self.q)                                                                       
        job.rebuild()

    def build(self):
        
        # ModelS
        select=self.datos[self.whichSelected(self.listbox)]  
        name= select.split()[0];       
    
        # Version
        version = self.whichSelected(self.listboxversion).__str__()  
        
        button_opt = {'fill': Tkconstants.BOTH, 'padx': 5, 'pady': 5}
        filebut = self.openFile()
        
        # Add argument to build list 
        self.seeds = [] 
        self.seeds.append(name)
        self.seeds.append(version)
        self.seeds.append(filebut)
        
        # Call new thread to build the model       
        self.button2.configure(state='disable')

        t = Thread(target=self.buildJob)
        t.start()   

 
if __name__ == "__main__":

    root = Tk()
##    root.geometry('1050x360+700+150')
    root.wm_title("etoxlab GUI (beta 2)")
    
    app = etoxlab(root)  
    root.mainloop()
