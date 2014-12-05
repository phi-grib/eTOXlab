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
import ttk

import tkMessageBox
import tkFileDialog
import commands
import os
import subprocess
import shutil
import Queue

from threading import Thread
from utils import wkd
import tarfile
from PIL import ImageTk, Image
import glob

'''
Creates a TreeView to shows general information about the existing models
'''
class modelViewer (ttk.Treeview):

    def __init__(self, parent):
              
        scrollbar_tree = ttk.Scrollbar(root)
        
        self.tree=ttk.Treeview.__init__(self, parent, columns = ('a','b','c','d','e'), selectmode='browse', yscrollcommand = scrollbar_tree.set)
               
        self.column("#0",minwidth=0,width=120, stretch=NO)
        self.column ('a', width=5)
        self.column ('b', width=50)
        self.column ('c', width=30)
        self.column ('d', width=30)
        self.column ('e', width=200)
        self.heading ('a', text='#')
        self.heading ('b', text='MD')
        self.heading ('c', text='mod')
        self.heading ('d', text='mol')
        self.heading ('e', text='quality')         
        
        self.datos= []

        scrollbar_tree.pack(side="left", fill=Y)
        self.pack(side="top", fill="both",expand=True,ipadx=100)

        # Move scrollbar 
        scrollbar_tree.config(command = self.yview)    
        self.bind('<<TreeviewSelect>>', self.chargeData())
                

    def clearTree(self):
        for i in self.get_children():
            self.delete(i)

    # Charges general information about models.    
    def chargeData(self):
        self.clearTree()
        version = []
        output=commands.getoutput(wkd+'/manage.py --info=short')
        output=output.split('------------------------------------------------------------------------------')
        
        for line in output[1:]:
            line=line.split("\n")
            self.datos.append(line)
            name=line[1].split("[")[0]            
            self.insert('', 'end', '%-9s'%(name), text='%-9s'%(name))
            count=0                            
            for x in line[2:-1]:                
                version=self.chargeDataDetails(x)
                if len(version)>0:
                    self.insert('%-9s'%(name), 'end', values=(version[0],version[1],version[2],version[3],version[4]), iid='%-9s'%(name)+str(count))
                    # The list of child is unfold in treeView
                    self.item('%-9s'%(name), open=True)                    
                count+=1              
                            
        # Focus in first element of the TreeView
        self.selection_set(self.get_children()[0:1])
        self.focus_set()
        self.focus(self.get_children()[0:1][0])
                       
        del version

    # Charges detailed information about each version of a given model(line).   
    def chargeDataDetails(self,line):
        
        l=line.split()
        y = []

        if (len(l) == 10):
            y.append('%-5s'%l[0])
            y.append('%-10s'%l[1].replace('MD:',''))
            y.append('%-8s'%l[2].replace('mod:',''))
            y.append('%-6s'%l[5])
            y.append('R2%-11s'%l[7]+"  "+'Q2%-11s'%l[9])

        elif (len(l) == 9):

             y.append('%-6s'%l[0])
             y.append('%-10s'%l[1].replace('MD:',''))
             y.append('%-8s'%l[2].replace('mod:',''))
             y.append('%-6s'%l[5])
             y.append('%-11s'%l[6]+'%-11s'%l[7]+"  "+'%-11s'%l[8])           

        else:
             line = 'not recognized'        
  
        return y
    
'''
Creates a new object to execute manage.py commands in new threads
'''
class Process:
    
    def __init__(self, parent, command, seeds, q, hasversion):
        
        self.model = parent
        self.command = command
        self.seeds = seeds
        self.queue = q
        self.dest =''
        self.v=hasversion
        
    def processJob(self):

        p = processWorker(self.command, self.seeds, self.queue,self.dest,self.v)                                                                       
        p.compute()                      

    def process(self):
        
        d=self.model.focus().split()
        self.seeds = []
        
        if len(d)==1:        
            self.seeds.append(d[0])
        else:
            self.seeds.append(d[0])
            self.seeds.append(d[1])
                
        try:            
            childs=self.model.get_children('%-9s'%d[0])                                
            if len(childs)==0 and (self.command==' --remove' or self.command==' --publish'):
                self.queue.put("The endpoint has no versions")
                
            elif self.command==' --remove' and len(childs)==1:
                self.queue.put("The 'sandbox' version cannot be deleted")
                
            else:
                if self.v:
                    self.dest=tkFileDialog.askdirectory(initialdir='.',title="Choose a directory...")                    

                t = Thread(target=self.processJob)
                t.start()
            
        except IndexError:
            self.queue.put("The model selected is not correct")
            pass
        
class processWorker: 

    def __init__(self, command, seeds, queue,d,hasversion):
        self.mycommand = command
        self.seeds = seeds
        self.q = queue
        self.dest = d
        self.v=hasversion

    def compute(self):

        try:
            endpoint = self.seeds[0]            

            if self.v:
                if self.dest=='':
                    self.q.put('Select a directory to save the data')  
                else:                    
                    version = self.seeds[1]
                    os.chdir(self.dest)
                    subprocess.call(wkd+'/manage.py -e '+endpoint+' -v '+version+' '+self.mycommand, stdout=subprocess.PIPE, shell=True)
                    self.q.put('Process finished')            
            else:
                subprocess.call(wkd+'/manage.py -e '+endpoint+' '+self.mycommand, stdout=subprocess.PIPE, shell=True)
                self.q.put('Process finished')
                
        except IndexError:
            self.q.put ('Process failed')
    
        
'''
Creates an object to execute view.py command in a new thread
'''
class Visualization:
    
    def __init__(self, parent, seeds, q):        
        self.model = parent
        self.seeds = seeds
        self.queue = q

    def viewJob(self):
        view = viewWorker(self.seeds, self.queue)                                                                       
        view.compute()                      

    def view(self):        
        d=self.model.focus().split()        
        try:           
            self.seeds = [] 
            self.seeds.append(d[0])
            self.seeds.append(d[1])

            # Call new thread to build the model       
            app.button4.configure(state='disable')
        
            t = Thread(target=self.viewJob)
            t.start()
            
        except IndexError:
            self.queue.put("The model selected has no version")
            pass

class viewWorker: 

    def __init__(self, seeds, queue):
        self.seeds = seeds
        self.q = queue

    def compute(self):        
        name    = self.seeds[0]
        version = self.seeds[1]
        
        mycommand=[wkd+'/view.py','-e',name,'-v',version]
        
        try:            
            subprocess.call(mycommand)                        
        except:
            self.q.put ('View process failed')

        self.q.put('View process completed')
        

'''
Creates an object to execute build.py command in a new thread
'''
class buildmodel:
    
    def __init__(self, parent, seeds, q):        
        self.model = parent
        self.seeds = seeds
        self.queue = q
          
    def buildJob(self):
        job = buildWorker(self.seeds, self.queue)
        job.rebuild()

    def build(self):
        
        d=self.model.focus().split()
        try:
            name=d[0]
            version=d[1]
            button_opt = {'fill': Tkconstants.BOTH, 'padx': 5, 'pady': 5}
            filebut = self.openFile()

            if filebut=='' or ".sdf" not in filebut:
                self.queue.put("Select a training set")
            else:                           
                # Add argument to build list 
                self.seeds = [] 
                self.seeds.append(name)
                self.seeds.append(version)
                self.seeds.append(filebut)
            
                # Call new thread to build the model       
                app.button2.configure(state='disable')
                app.pb.start(100)
            
                t = Thread(target=self.buildJob)
                t.start()
            
        except IndexError:
            self.queue.put("The model selected has no version")            
            pass

    def openFile(self):
        filename = tkFileDialog.askopenfilename(parent=root,
                                                initialdir='/home/modeler/workspace',
                                                title="Select a new training set to rebuild the model",
                                                filetypes=[('image files', '.sdf')]) 
        return filename
    

class buildWorker: 

    def __init__(self, seeds, queue):
        self.seeds = seeds
        self.q = queue

    def rebuild(self):      
        name    = self.seeds[0]
        version = self.seeds[1]
        filebut = self.seeds[2]

        # Save all in a specific workspace folder
        self.tempdir='/home/modeler/workspace/'+name         
                  
        try:
            if not os.path.isdir(self.tempdir): os.mkdir (self.tempdir)
        except:
            return (False,'unable to create directory '+self.tempdir)     
        
        try:
            shutil.copy(filebut,self.tempdir+'/training.sdf')
        except:
            return (False,'unable to copy training series to '+self.tempdir)     
    
        os.chdir(self.tempdir)
        mycommand = [wkd+'/build.py','-e',name,'-f','training.sdf','-v',version]

        # Rebuild the model and save the result in a ".txt" file
        f = open(wkd+"/output.txt", "w")
        try:
            subprocess.call(mycommand,stdout=f)
        except:
            self.q.put ('Building failed')
        f.close()

        # Check the model has been generated correctly
        f1 = open(wkd+"/output.txt", "r")

        lines = f1.readlines()
        if not "Model OK" in lines[-1]:
            self.q.put('Building failed')
        else:
            self.q.put('Building completed')    
     
    
'''
GUI class
'''
class etoxlab:

    def __init__(self, master):
        self.seeds  = []
        self.q      = Queue.Queue()
        self.master = master
        self.myfont = 'Courier New'

        self.f10 = (self.myfont,10)

        # create a toplevel menu
        menubar = Menu(root)
        filemenu = Menu(menubar)    
        filemenu.add_command(label="Exit", command= quit)
        #filemenu.add_command(label="Exit", command= root.withdraw) this doesn't work in standalone mode!
        filemenu1 = Menu(menubar)   
        filemenu1.add_command(label="About eTOXlab", command=self.showhelp)
        menubar.add_cascade(label="File", menu=filemenu)
        menubar.add_cascade(label="Help", menu=filemenu1)
        root.config(menu=menubar)
             
        i1 = Frame(root)
        i2 = Frame(root)        
        
        ## Treeview
        t1 = Frame(i1)
        self.models = modelViewer (t1)          
     
        # Main container is a Notebook
        n = ttk.Notebook (i2)
        f1 = Frame(n)
        f2 = Frame(n)
        f3 = Frame(n)        
        
        f1.pack(side="top", fill="both", expand=False)
        f2.pack(side="top", fill="both", expand=False)
        f3.pack(side="top", fill="both", expand=False)
        
        n.add (f1, text='manage')
        n.add (f2, text='build')
        n.add (f3, text='view')
        n.pack (side="top", fill="x", expand=False)
       
        # MANAGE Frame        
        ## Buttons
        f12 = ttk.Frame(f1)
    
        fnew = ttk.LabelFrame(f12, text='new endpoint')
        fnewi= Frame(fnew)
        fnew1 = ttk.Frame(fnewi)
        fnew2 = ttk.Frame(fnewi) 
       
        lnew1 = Label(fnew1, text='%-6s'%'name')        
        self.enew1 = Entry(fnew1, bd =1)
        lnew2 = Label(fnew2, text='%-9s'%'tag')
        self.enew2 = Entry(fnew2, bd =1)
        
        lnew1.pack(side="left")
        self.enew1.pack(side="right")
        lnew2.pack(side="left")
        self.enew2.pack(side="right")
        fnew1.pack(side="top")
        fnew2.pack(side="top")
        
        fnewj= Frame(fnew)
        lnew3 = Label(fnewj, text='creates a new endpoint')
        lnew3.pack(side="left", fill="y", padx=5, pady=5)
        
        button10 = Button(fnewj, text ='OK', command = self.new, height=0, width=5)
        button10.pack(side="right", expand=False, padx=5, pady=5)
        fnewi.pack(side="top", fill="x", expand=False)
        fnewj.pack(side="top", fill="x", expand=False)        
        fnew.pack(side="top", fill="x", expand=False, padx=5, pady=5)

        finfo = LabelFrame(f12, text='get information')
        linfo1 = Label(finfo, text='shows complete model info')
        linfo1.pack(side="left", fill="x", padx=5, pady=5)
        self.button11 = Button(finfo, text ='OK', command = self.seeDetails,height=0, width=5)
        self.button11.pack(side="right", fill="y",expand=False, padx=5, pady=5)        
        finfo.pack(side="top", fill="x", padx=5, pady=12)

        self.publish=Process(self.models,' --publish', self.seeds, self.q, FALSE) 

        fpublish = LabelFrame(f12, text='publish model')
        lpublish1 = Label(fpublish, text='creates a new model version')
        lpublish1.pack(side="left", fill="x",padx=5, pady=5)
        self.button12 = Button(fpublish, text ='OK', command = self.publish.process,height=0, width=5)
        self.button12.pack(side="right", fill="y", expand=False, padx=5, pady=5)
        fpublish.pack(side="top", fill="x", padx=5, pady=12)

        self.remove=Process(self.models,' --remove', self.seeds, self.q, FALSE)
        
        frem = LabelFrame(f12, text='remove model')
        lrem1 = Label(frem, text='removes a model version')
        lrem1.pack(side="left", fill="x",padx=10, pady=5)
        self.button10 = Button(frem, text ='OK', command = self.remove.process, height=0, width=5)
        self.button10.pack(side="right", fill="y", expand=False, padx=5, pady=5)
        frem.pack(side="top", fill="x", padx=5, pady=10) 

        self.gseries=Process(self.models,' --get=series', self.seeds,self.q, TRUE)
        
        fgets = LabelFrame(f12, text='get series')
        lgets1 = Label(fgets, text='retrieves the training series')
        lgets1.pack(side="left", fill="x",padx=10, pady=5)
        self.button13 = Button(fgets, text ='OK', command = self.gseries.process, height=0, width=5)
        self.button13.pack(side="right", fill="y", expand=False, padx=5, pady=5)
        fgets.pack(side="top", fill="x", padx=5, pady=10)

        self.gmodel=Process(self.models,' --get=model', self.seeds, self.q, TRUE)

        fgetm = LabelFrame(f12, text='get model')
        lgetm1 = Label(fgetm, text='retrieves the model definition')
        lgetm1.pack(side="left", fill="x",padx=10, pady=5)
        self.button14 = Button(fgetm, text ='OK', command = self.gmodel.process, height=0, width=5)
        self.button14.pack(side="right", fill="y", expand=False, padx=5, pady=5)
        fgetm.pack(side="top", fill="x", padx=5, pady=10)

        self.export=Process(self.models,' --export',self.seeds,self.q,TRUE)
        
        fexp = LabelFrame(f12, text='export')
        lexp1 = Label(fexp, text='packs current model version')
        lexp1.pack(side="left", fill="x",padx=10, pady=5)
        self.button10 = Button(fexp, text ='OK', command = self.export.process, height=0, width=5)
        self.button10.pack(side="right", fill="y", expand=False, padx=5, pady=5)
        fexp.pack(side="top", fill="x", padx=5, pady=10)        
        
        f12.pack(side="top", fill="both",expand=False)
        
        # BUILD Frame        
        self.bmodel=buildmodel(self.models, self.seeds,self.q) 
        
        f22 = Frame(f2)
        fbuild = LabelFrame(f22, text='rebuild model')

        lbuild1 = Label(fbuild, text='re-train current model')
        lbuild1.pack(side="left", fill='x', ipadx=10, expand=False)
        self.button2 = Button(fbuild, text = 'OK', command = self.bmodel.build, height=0, width=5)
        self.button2.pack(side="right",fill='y',padx=5, pady=5, expand=False)
        fbuild.pack(side="top", fill='x', padx=5, pady=5, expand=False)

        self.pb = ttk.Progressbar(f22, orient='horizontal', mode='indeterminate', value=0)
        self.pb.pack(side="top", fill='x', expand=False)       

        f22.pack(side="top", fill="x", expand=False)
 
        # VIEW Frame
        self.view=Visualization(self.models,self.seeds,self.q)

        f32 = Frame(f3)
        fview = LabelFrame(f32, text='view model')
        
        lview1 = Label(fview, text='generates graphic representations')
        lview1.pack(side="left", fill='x', ipadx=10, expand=False)
        self.button4 = Button(fview, text = 'OK', command = self.view.view, height=0, width=5)
        self.button4.pack(side="right",fill='y',padx=5, pady=5, expand=False)
        fview.pack(side="top", fill='x', padx=5, pady=5, expand=False)
        
        f32.pack(side="top",fill='x', expand=False)

        t1.pack(side="left", fill="both",expand=True)
        i1.pack(side="left", fill="both",expand=True)
        i2.pack(side="right", fill="both",expand=False)

       # Start queue listener
        self.periodicCall()

    '''
    Help window
    '''
    def showhelp(self):
        win = Tk()                       
        win.wm_title('About')

        t_help = Text(win)
        t_help.insert(INSERT, "Description    eTOXlab simple GUI\n"+
                    "Authors:       Ines Martinez and Manuel Pastor (manuel.pastor@upf.edu)\n"+
                    "Copyright 2014 Manuel Pastor"+
                    "\n\neTOXlab is free software: you can redistribute it and/or modify"+
                    "it under the terms of the GNU General Public License as published by"+
                    "the Free Software Foundation version 3.")
        t_help.config(state=DISABLED)
        t_help.pack(side="top", fill="both", expand=True)          

        win.geometry('330x180-5+40')
        win.mainloop()        
        
    '''
    Creates a new endpoint.
    '''
    def new(self):    
        endpoint = self.enew1.get()
        tag = self.enew2.get()
                
        if not endpoint:
            self.q.put ('Please provide the name of the endpoint')

        elif not tag:
           self.q.put ('Please provide the label of the eTOXsys web service')

        else:
            subprocess.call(wkd+'/manage.py -e '+endpoint+' -t '+tag+' --new', stdout=subprocess.PIPE, shell=True)
            self.enew1.delete(0,END)
            self.enew2.delete(0,END)
            self.q.put('New endpoint created')
            
        #Refresh TreeView
        self.models.chargeData()        

    '''
    Presents information about the model defined by the endpoint
    If no version is specified (general endpoint), shows all the model versions.
    '''    
    def seeDetails(self):
    
        #Get the selected model
        d=self.models.focus().split()
        name =[]
        version= []
        try:
            name=d[0]
            version=d[1]
        except IndexError:
            pass

        # Obtain information about the model
        if len(version) == 0:            
            output=commands.getoutput(wkd+'/manage.py -e '+name+' --info=long')
            outputlist = output.split('\n')
            outputlist = outputlist [1:-1]
        else:
            output=commands.getoutput(wkd+'/manage.py -e '+name+' -v' +version  +' --info=long')
            outputlist = output.split('\n')
            outputlist = outputlist [1:-3]         
        
        output = ''
        for l in outputlist: output+= l+'\n'
        
        # Show the collected information in a new window (winDetails)
        winDetails = Tk()
        winDetails.resizable(0.5,0.5)
                        
        if len(version) == 0:
            winDetails.wm_title(name+': All models')
        else:
            winDetails.wm_title(name+' ver '+version)

        scrollbar = Scrollbar(winDetails,orient=VERTICAL)
        scrollbar.pack(side=RIGHT, fill=Y)

        text = Text(winDetails, wrap=WORD, font=self.f10, yscrollcommand=scrollbar.set)
        text.insert(INSERT, output)
        text.config(state=DISABLED)
        text.pack(side="top", fill="both", expand=True)
        
        scrollbar.config(command=text.yview)     
       
        winDetails.mainloop()
       
    '''
    Handle all the messages currently in the queue (if any)
    and run events depending of its information
    '''
    def periodicCall(self):
    
        while self.q.qsize():
        
            try:
                msg = self.q.get(0)
                if 'Building' in msg:            
                    self.button2.configure(state='normal')
                    self.pb.stop()                    
                    if 'completed' in msg:
                        d=self.models.focus().split()
                        key="pls"
                        self.win=visualizewindow(d[0],key)
                        self.win.viewmodels()
                    self.models.chargeData()
                    
                if 'View' in msg:
                    self.button4.configure(state='normal')
                    key="view"
                    d=self.models.focus().split()
                    self.win=visualizewindow(d[0])
                    self.models.chargeData()

                if 'finished' in msg:
                    self.models.chargeData() 
                      
                tkMessageBox.showinfo("Info Message", msg)

            except Queue.Empty:
                pass

        self.master.after(500, self.periodicCall) # re-call after 500ms
        
'''
Creates a new window that contains the png files included in a directory(ori)
'''
class visualizewindow(Toplevel):
    
    def __init__(self,ori,key):        
        Toplevel.__init__(self)
        self.ori=ori   # in which directory search
        self.key=key   # identifier of the files
        
    def viewmodels(self):
        
        pngfiles= glob.glob("/home/modeler/workspace/"+self.ori+"/*"+self.key+"*.png")

        if len(pngfiles)==0:
            self.destroy()
        else:
            note_view = ttk.Notebook(self)
            note_view.pack()
            
            for t in pngfiles:
                f = Frame(self)
                note_view.add(f,text=t.split("/")[5])
                  
                i = ImageTk.PhotoImage(Image.open(t))
        
                l1=ttk.Label(f,image=i)        
                l1.image= i
                l1.pack()           
    

if __name__ == "__main__":

    root = Tk()
    root.wm_title("etoxlab GUI (beta 0.71)")    

    app = etoxlab(root)
    root.mainloop()
