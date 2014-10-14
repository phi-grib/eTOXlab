#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    eTOXlab simple GUI (OOP version)
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

import Tkinter
import time
import threading
import random
import Queue
import Tkconstants
import tkMessageBox
import tkFileDialog
import commands
import os
import subprocess
import shutil
import time
from fabric.colors import white
from Tkinter import * 


class GuiPart:
    def __init__(self, master, queue, endCommand):
        global scroll,listbox,datos,datos1,listboxversion,root
        datos = []
        datos1 = []
        
        self.queue = queue
        
        # Set up the GUI                      
        
        root.geometry('1050x330+700+150')
        root.wm_title("etoxdesk")
    
        # m0 General PanedWindows that contains m1 and m2
        m0 = PanedWindow()
        m0.pack(fill=BOTH, expand=1)
        
        # Model details - Paned m1
        m1 = PanedWindow(m0,orient=VERTICAL)    
        
        labelList = Label(root, text='#endpoint   tag', anchor=W, justify=LEFT,font=("Courier New", 10))
        listbox = Listbox(root,width=70, height=20,exportselection=False,font=("Courier New", 10))
        
        m1.add(labelList)    
        m1.add(self.chargeData())
    
        m0.add(m1)
        
        # Version details - Paned m2
        labelListv = Label(root, text="#    MD            mod    mol", anchor=W, justify=LEFT,font=("Courier New", 10))
        
        #List of Versions
        listboxversion = Listbox(root,width=70, font=("Courier New", 9))        
        m2 = PanedWindow(m0, orient=VERTICAL)  
        listbox.bind('<<ListboxSelect>>', self.chargeVersion) 
        m2.add(labelListv)
        m2.add(listboxversion)
    
        self.chargeVersion(0)
            
        # Buttons
        button1 = Button(root, text = 'model details', command = self.seeDetails,font=("Courier New", 10))

        button2 = Tkinter.Button(master, text='Build', command=endCommand)         
        button3 = Button(root, text = 'publish', command = self.publish,font=("Courier New", 10))
        button4 = Button(root, text = 'Quit', command = quit,font=("Courier New", 10))
        
        m2.add(button1)
        m2.add(button2)
        m2.add(button3)  
        m2.add(button4)
        m0.add(m2)        
        

    def processIncoming(self):
        """
        Handle all the messages currently in the queue (if any).
        """
        while self.queue.qsize():
            try:
                msg = self.queue.get(0)
                # Check contents of message and do what it says
                # As a test, we simply print it
                print msg
            except Queue.Empty:
                pass
            
    def publish(self):   
        # Choose the model
        select = datos[self.whichSelected(listbox)]  
        name = select.split()[0];
        
        subprocess.call('/home/modeler/soft/eTOXlab/src/manage.py --publish -e '+name, stdout=subprocess.PIPE, shell=True)  
    
        self.chargeVersion(0)   
    
    def chargeData(self):
        #Read a file and add all the data to a file  
        output=commands.getoutput('python /home/modeler/soft/eTOXlab/src/manage.py --info=short')
            
        output=output.split("\n")         
       
        for line in output[1:-1]:
            if line.find('MD')==-1:
                if line.find('-')==-1: 
                    datos.append(line) 
                    line=line.split("[")
                    listbox.insert(END, '%-12s'%(line[0])+ '%-30s'%(line[1][0:-1]))
                    
        
        #First item selected    
        listbox.selection_set( first = 0 )
         
        return listbox
    
    def chargeVersion(self,e):
    
        listboxversion.delete(0, END)
         
        d=datos[self.whichSelected(listbox)]    
        name= d.split()[0];
        
        # Obtain information about the model      
        output=commands.getoutput('/home/modeler/soft/eTOXlab/src/manage.py -e '+name+' --info=short') 
        
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
            else :
                line =  '%-6s'%l[0] \
                        +'%-16s'%l[1].replace('MD:','') \
                        +'%-8s'%l[2].replace('mod:','') \
                        +'%-6s'%l[5] \
                        +'%-11s'%l[6] \
                        +'%-11s'%l[7] \
                        +'%-11s'%l[8]
                
            listboxversion.insert(END, line)
    
        
        listboxversion.selection_set( first = 0 )
       
        return listboxversion
                
    def whichSelected(self,l):
        return int(l.curselection()[0])
    
    def seeDetails(self):
        # Capture the name of the selected model
        d=datos[self.whichSelected(listbox)]
        name= d.split()[0];
        
        version = self.whichSelected(listboxversion).__str__()    
    
        # Obtain information about the model
        output=commands.getoutput('/home/modeler/soft/eTOXlab/src/manage.py -e '+name+' -v' +version  +' --info=long')
    
        outputlist = output.split('\n')
        outputlist = outputlist [1:-3]
    
        output = ''
        for l in outputlist: output+= l+'\n'
        
        # Show the collected information in a new window (winDetails)
        winDetails = Tk()
        winDetails.resizable(0.5,0.5)
        win1 = PanedWindow()
        winDetails.wm_title(name+' ver '+version)
        
        lab= Label(winDetails, text = output,font=("Courier New", 12),justify=LEFT)
        lab.pack()
        winDetails.mainloop()

class ThreadedClient:
    """
    Launch the main part of the GUI and the worker thread. periodicCall and
    endApplication could reside in the GUI part, but putting them here
    means that you have all the thread controls in a single place.
    """
    def __init__(self, master):
        """
        Start the GUI and the asynchronous threads. We are in the main
        (original) thread of the application, which will later be used by
        the GUI. We spawn a new thread for the worker.
        """
        self.master = master

        # Create the queue
        self.queue = Queue.Queue()

        # Set up the GUI part
        self.gui = GuiPart(master, self.queue, self.starbuild()) 
      

    def starbuild(self):
        # Set up the thread to do asynchronous I/O
        # More can be made if necessary        
        self.thread1 = threading.Thread(target=self.reBuild())
        self.thread1.start()
        
            
    def reBuild(self):
        
        # Choose the model
        select=self.gui.datos[self.gui.whichSelected(self.gui.listbox)]
        name= select.split()[0];
    
        # Choose the version
        version = self.gui.whichSelected(listboxversion).__str__()    
        
        
        button_opt = {'fill': Tkconstants.BOTH, 'padx': 5, 'pady': 5}
        filebut = self.openFile()
        # Save all in a specific workspace folder
        wkd='/home/modeler/workspace/'+name
    
        try:
            if (os.path.isdir(wkd) == 0):
                os.mkdir (wkd)
        except:
            return (False,'unable to create directory '+wkd)
        
        try:
            shutil.copy(filebut,wkd+'/training.sdf')
        except:
            return (False,'unable to copy training at '+wkd)
    
        os.chdir(wkd)
        
        mycommand = ['python','/home/modeler/soft/eTOXlab/src/build.py','-e',name,
                   '-f',filebut,'-v',version]
        
        print mycommand
#         with aLock:
        
        subprocess.Popen(mycommand)
        
        tkMessageBox.showinfo("Info Message", "Rebuilding finished")
   
    def openFile(self):
        filename = tkFileDialog.askopenfilename(parent=root,initialdir='/home/modeler',title="Select a new training set to rebuild the model" , filetypes=[('image files', '.sdf')]) ## filename not filehandle
        return filename
        
        
rand = random.Random()
root = Tkinter.Tk()

client = ThreadedClient(root)

root.mainloop()
