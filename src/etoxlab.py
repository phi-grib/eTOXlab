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


################################################################
### MANAGE
################################################################
    
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
        self.seeds = []        
        self.seeds.append (self.model.selEndpoint())
        self.seeds.append (self.model.selVersion())
                            
        if self.v:
            self.dest=tkFileDialog.askdirectory(initialdir='.',title="Choose a directory...")                    

        t = Thread(target=self.processJob)
        t.start()
            
        
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
##            print "Process: "
##            print self.seeds

            if self.v:
                if self.dest=='':
                    self.q.put('Select a directory to save the data')  
                else:                    
                    version = self.seeds[1]
                  #  print version
                    os.chdir(self.dest)
                    subprocess.call(wkd+'/manage.py -e '+endpoint+' -v '+version+' '+self.mycommand, stdout=subprocess.PIPE, shell=True)
                    self.q.put('Process finished')            
            else:
                subprocess.call(wkd+'/manage.py -e '+endpoint+' '+self.mycommand, stdout=subprocess.PIPE, shell=True)
                self.q.put('Process finished')
                
        except IndexError:
            self.q.put ('Process failed')

################################################################
### VIEW
################################################################
        
'''
Creates an object to execute view.py command in a new thread
'''
class Visualization:
    
    def __init__(self, parent, seeds, q):        
        self.model = parent 
        self.seeds = seeds
        self.queue = q

    def viewJob(self):
        view = viewWorker(self.seeds, self.queue, self.vtype, self.molecules, self.background, self.refname, self.refver)                                                                       
        view.compute()                      

    def view(self):        
        self.seeds = [] 
        self.seeds.append(self.model.selEndpoint())
        self.seeds.append(self.model.selVersion())
##        print "View: "
##        print self.seeds

        # Call new thread to visualize the series       
        app.button4.configure(state='disable')

        self.vtype   = app.viewTypeCombo.get()
        self.refname = app.referEndpointCombo.get()
        self.refname = self.refname.strip()

        if self.refname=='None':
            self.refname = ''
            
        refverstr = app.referVersionCombo.get()
        refverstr = refverstr.strip()
        
        try:            
            self.refver  = int(refverstr)
        except:
            self.refver = 0     
       
        self.background = (app.viewBackground.get() == '1')
        self.molecules = ''
        
        t = Thread(target=self.viewJob)
        t.start()

    def viewQuery(self):
        self.seeds = [] 
        self.seeds.append(self.model.selEndpoint())
        self.seeds.append(self.model.selVersion())

        # Call new thread to visualize the series       
        app.buttonQuery4.configure(state='disable')

        self.vtype       = app.viewTypeComboQuery.get()
        self.refname     = ''
        self.refver      = ''
        self.molecules   = app.eviewQuery1.get()      
        self.background  = (app.viewBackgroundQuery.get() == '1')
        
        t = Thread(target=self.viewJob)
        t.start()

    def viewModel(self):
        endpointDir = wkd+'/' + app.models.selEndpoint() +'/version%0.4d'%int(app.models.selVersion())
        files = [endpointDir+'/recalculated.png', endpointDir+'/predicted.png']
        self.win=visualizewindow()
        self.win.viewMultiple(files)
            

class viewWorker: 

    def __init__(self, seeds, queue, vtype, molecules, background, refname, refver):
        self.seeds = seeds
        self.q = queue
        self.vtype = vtype
        self.molecules = molecules
        self.background = background
        self.refname = refname
        try:
            self.refver = int(refver)
        except:
            self.refver = 0
            
    def compute(self):        
        name    = self.seeds[0]
        version = self.seeds[1]

        mycommand=[wkd+'/view.py','-e',name,'-v',version,
                   '--type='+self.vtype]

        if len(self.molecules)>0:
            mycommand.append('-f')
            mycommand.append(self.molecules)

        if len(self.refname)>0:
            mycommand.append('--refname=' +self.refname)

        if self.refver != None:
            mycommand.append('--refver=%d' %self.refver)
            
        if self.background :
            mycommand.append('--background')
            
        try:            
            result = subprocess.call(mycommand)                        
        except:
            self.q.put ('View process failed')

        if result == 1 :
            self.q.put ('View process failed')
            return
        
        if self.vtype=='pca' or self.vtype=='project':
            outname = './pca-scores12.png'
        else:
            outname = './generic.png'
            
        self.q.put('View completed '+outname)


################################################################
### BUILD
################################################################
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
        name = self.model.selEndpoint()
        version = self.model.selVersion()
        filebut = app.buildSeries.get()

        if filebut=='' or ".sdf" not in filebut:
            self.queue.put("Select a training set")
            return
                        
        # Add argument to build list 
        self.seeds = [] 
        self.seeds.append(name)
        self.seeds.append(version)
        self.seeds.append(filebut)

##        print "Build: "
##        print self.seeds
    
        # Call new thread to build the model       
        app.button2.configure(state='disable')
        app.pb.start(100)
    
        t = Thread(target=self.buildJob)
        t.start()

class buildWorker: 

    def __init__(self, seeds, queue):
        self.seeds = seeds
        self.q = queue

    def rebuild(self):      
        name    = self.seeds[0]
        version = self.seeds[1]
        filebut = self.seeds[2]

        mycommand = [wkd+'/build.py','-e',name,'-v',version,'-f',filebut]

        # Rebuild the model and save the result in a ".txt" file
        f = open(wkd+"/build.log", "w")
        try:
            subprocess.call(mycommand,stdout=f)
        except:
            self.q.put ('Building failed')
        f.close()

        # Check the model has been generated correctly
        f1 = open(wkd+"/build.log", "r")

        lines = f1.readlines()
        if not "Model OK" in lines[-1]:
            self.q.put('Building failed')
        else:
            self.q.put('Building completed')    


################################################################
### TREEVIEW CLASS
################################################################

'''
Creates a TreeView to shows general information about the existing models
'''
class modelViewer (ttk.Treeview):

    def __init__(self, parent):
              
        scrollbar_tree = ttk.Scrollbar(root)
        
        self.tree=ttk.Treeview.__init__(self, parent, columns = ('a','b','c','d','e'), selectmode='browse', yscrollcommand = scrollbar_tree.set)
               
        self.column("#0",minwidth=0,width=150, stretch=NO)
        self.column ('a', width=5)
        self.column ('b', width=30)
        self.column ('c', width=80)
        self.column ('d', width=20)
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
        self.chargeData()
        
    def clearTree(self):
        for i in self.get_children():
            self.delete(i)

    def selEndpoint(self):       
        return self.focus().split()[0]

    def selVersion(self):
        d = self.set(self.focus()).get('a')
        
        if d:
            d=d.strip()
            if d=='*':
                return ('0')
            else:
                return d
        else:
            return ('0')
      
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

        self.maxver = 1
        for child in self.get_children():
            iver= len(self.get_children(child))
            if iver > self.maxver : self.maxver=iver
        
        # Focus in first element of the TreeView
        self.selection_set(self.get_children()[0:1])
        self.focus(self.get_children()[0:1][0])

    # Charges detailed information about each version of a given model(line).   
    def chargeDataDetails(self,line):

        line=line.replace(':',' ')
        l=line.split()
        #print l, len(l)
        y = []
        y.append('%-5s' %l[0])         # num 
        y.append('%-10s'%l[2])         # MD
        
        i=4                            # regression method name is often split
        imethod=''
        while l[i]!='mol':  # label of next field
            imethod+=l[i]
            imethod+=' '
            if i == len(l):
                return y
            else:
                i=i+1
               
        y.append ('%-8s'%imethod)     # regression method

        if l[i+1].isdigit():
            y.append ('%-6s'%l[i+1])  # num mol
        else:
            y.append ('na')

        if 'R2' in l:                 # quality for quantitative endpoints
            if len(l) < i+5 :
                y.append ('na')
            else:
                y.append('R2:%-11s'%l[i+3]+' Q2:%-11s'%l[i+5])
        elif 'MCC' in l:              # quality for qualitative endpoints
            if len(l) < i+7 :
                y.append ('na')
            else:
                y.append('sen:%-11s'%l[i+3]+'spe:%-11s'%l[i+5]+'MCC:%-11s'%l[i+7])        
        else:                         # fallback
             line = 'not recognized'        
  
        return y


#####################################################################################################
### VIEWER WINDOW CLASS
#####################################################################################################        
        
'''
Creates a new window that displays one or more plots given as a list of png files
'''
class visualizewindow(Toplevel):
    
    def __init__(self):        
        Toplevel.__init__(self)
        
    def viewSingle(self, fname):
        
        if fname==None or fname=='':
            self.destroy()
            return

        if not os.path.isfile (fname):
            self.destroy()
            return
        
        f = Frame(self)
        self.i = ImageTk.PhotoImage(Image.open(fname))
        ttk.Label(f,image=self.i).pack()        
        f.pack()

    def viewMultiple (self, fnames):
        
        if fnames==None or len(fnames)==0:
            self.destroy()
            return
        
        self.note_view = ttk.Notebook(self)
        self.note_view.pack()

        self.i=[]
        for t in fnames:
            if not os.path.isfile (t) : continue
            self.i.append (ImageTk.PhotoImage(Image.open(t)))
            
            f = Frame(self)
            self.note_view.add(f,text=t.split('/')[-1])
            ttk.Label(f,image=self.i[-1]).pack()

        if not len(self.i):
            self.destroy()
            return
        
        self.note_view.pack()
        
     
###################################################################################
### MAIN GUI CLASS
###################################################################################   
'''
GUI class
'''
class etoxlab:

    def __init__(self, master):
        self.seeds  = []
        self.q      = Queue.Queue()
        self.master = master
        self.myfont = 'Courier New'

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
             
        i1 = Frame(root) # frame for tree (common to all tabs)
        i2 = Frame(root) # frame for notebook
        
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
       
        ## MANAGE Frame        
        f12 = ttk.Frame(f1)
    
        fnew = ttk.LabelFrame(f12, text='new endpoint')
        fnewi = ttk.Frame(fnew)
        fnew1 = ttk.Frame(fnewi)
        fnew2 = ttk.Frame(fnewi) 
       
        lnew1 = Label(fnew1, width = 10, anchor='e', text='name')        
        self.enew1 = Entry(fnew1, bd =1)               # field containing the new endpoint name
        lnew2 = Label(fnew2, width = 10, anchor='e', text='tag')
        self.enew2 = Entry(fnew2, bd =1)               # field containing the new endpoint tag
       
        lnew1.pack(side="left")
        self.enew1.pack(side="right")
        lnew2.pack(side="left")
        self.enew2.pack(side="right")
        fnew1.pack(side="top", anchor ='w')
        fnew2.pack(side="top", anchor ='w')
        
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
        lexp1 = Label(fexp, text='packs whole model tree')
        lexp1.pack(side="left", fill="x",padx=10, pady=5)
        self.button10 = Button(fexp, text ='OK', command = self.export.process, height=0, width=5)
        self.button10.pack(side="right", fill="y", expand=False, padx=5, pady=5)
        fexp.pack(side="top", fill="x", padx=5, pady=10)        
        
        f12.pack(side="top", fill="both",expand=False)
        
        ## BUILD Frame        
        self.bmodel=buildmodel(self.models, self.seeds,self.q) 
        
        f22 = Frame(f2)

        fbuild = LabelFrame(f22, text='build model')

        fbuild0 = Frame(fbuild)
        fbuild1 = Frame(fbuild)
        
        # frame 1: entry field for selecting new training series
        lbuild0 = Label(fbuild0, width = 10, anchor='e', text='new series')        
        self.buildSeries = Entry(fbuild0, bd =1)               # field containing the new training name
        lbuild0.pack(side="left")
        self.buildSeries.pack(side="left")
        Button(fbuild0, text ='...', width=2, command = self.selectTrainingFile).pack(side="right")

        lbuild1 = Label(fbuild1, text='train selected model')
        lbuild1.pack(side="left", anchor="w", padx=10, expand=False)
        self.button2 = Button(fbuild1, text = 'OK', command = self.bmodel.build, height=0, width=5)
        self.button2.pack(side="right", anchor="e", padx=5, pady=5, expand=False)

        fbuild0.pack(side="top", anchor="w")
        fbuild1.pack(side="top", anchor="e")
        
        fbuild.pack(side="top", fill='x', padx=5, pady=5, expand=True)

        self.pb = ttk.Progressbar(f22, orient='horizontal', mode='indeterminate', value=0)
        self.pb.pack(side="top", fill='x', expand=False)       

        f22.pack(side="top", fill="x", expand=False)
 
        ## VIEW Frame
        self.view=Visualization(self.models,self.seeds,self.q)

        f32 = Frame(f3)

        fview = LabelFrame(f32, text='view series')
                
        fview0 = Frame(fview)
        fview1 = Frame(fview)
        fview2 = Frame(fview)
        fview3 = Frame(fview)
        fviewi = Frame(fview)

        # frame 0: combo-box for seletig view type
        lview0 = Label (fview0, width = 10, anchor='e', text='type')
        self.viewTypeCombo = StringVar()
        self.cboCombo = ttk.Combobox( fview0, values=('pca','property','project'), textvariable=self.viewTypeCombo, state='readonly')
        self.cboCombo.current(0)
        lview0.pack(side='left')
        self.cboCombo.pack(anchor ='w')

        # frame 1: entry field for selecting reference endpoint        
        lview1 = Label(fview1, width = 10, anchor='e', text='refername')
        self.referEndpointCombo = StringVar ()
        comboValues=("None",) + self.models.get_children()
           
        self.eview1 = ttk.Combobox( fview1, values=comboValues, textvariable=self.referEndpointCombo, state='readonly')
        self.eview1.current(0)
        lview1.pack(side="left")
        self.eview1.pack()

        # frame 2: entry field for selecting reference version

        lview2 = Label(fview2, width = 10, anchor='e', text='refver')
        self.referVersionCombo = StringVar ()
        
        comboVersions = ()
        for i in range(self.models.maxver):
            comboVersions=comboVersions+(str(i),)  # this is updated by updateGUI method
           
        self.eview2 = ttk.Combobox( fview2, values=comboVersions, textvariable=self.referVersionCombo, state='readonly')
        self.eview2.current(0)
        lview2.pack(side="left")
        self.eview2.pack()
        
        # frame 3: check button for showing background
        lview3 = Label (fview3, width = 10, anchor='e', text='   ')
        self.viewBackground = StringVar()
        self.checkBackground = ttk.Checkbutton(fview3, text='show background', variable=self.viewBackground)
        self.viewBackground.set(0)
        lview3.pack(side="left")
        self.checkBackground.pack(anchor='w')

        # frame button 
        lviewi = Label(fviewi, anchor = 'w', text='represents graphically the training series')
        lviewi.pack(side="left", fill="y", padx=5, pady=5)        
        self.button4 = Button(fviewi, text ='OK', width=5, command = self.view.view)
        self.button4.pack(side="right", anchor='e', padx=5, pady=5)

        fview0.pack(side="top", expand=YES, anchor='w')
        fview1.pack(side="top", expand=YES, anchor='w')
        fview2.pack(side="top", expand=YES, anchor='w')
        fview3.pack(side='top', expand=YES, anchor='w')
        fviewi.pack(side="top", expand=YES, anchor="e")

        fview.pack(side="top", fill="x", expand=False, padx=5, pady=5)

        fviewQuery = LabelFrame(f32, text='view query')
                
        fviewQuery0 = Frame(fviewQuery)
        fviewQuery1 = Frame(fviewQuery)
        fviewQuery2 = Frame(fviewQuery)
        fviewQueryi = Frame(fviewQuery)

        # frame 0: combo-box for seletig view type
        lviewQuery0 = Label (fviewQuery0, width = 10, anchor='e', text='type')
        self.viewTypeComboQuery = StringVar()
        self.cboComboQuery = ttk.Combobox( fviewQuery0, values=('pca','property','project'), textvariable=self.viewTypeComboQuery, state='readonly')
        self.cboComboQuery.current(0)
        lviewQuery0.pack(side='left')
        self.cboComboQuery.pack(anchor ='w')

        # frame 1: entry field for selecting reference endpoint
        lviewQuery1 = Label(fviewQuery1, width = 10, anchor='e', text='query')        
        self.eviewQuery1 = Entry(fviewQuery1, bd =1)               # field containing the new endpoint name
        lviewQuery1.pack(side="left")
        self.eviewQuery1.pack(side="left")
        Button(fviewQuery1, text ='...', width=2, command = self.selectQueryFile).pack(side="right") 

        # frame 2: check button for showing background
        lviewQuery2 = Label (fviewQuery2, width = 10, anchor='e', text='   ')
        self.viewBackgroundQuery = StringVar()
        self.checkBackgroundQuery = ttk.Checkbutton(fviewQuery2, text='show background', variable=self.viewBackgroundQuery)
        self.viewBackgroundQuery.set(0)
        lviewQuery2.pack(side="left")
        self.checkBackgroundQuery.pack(anchor='w')

        # frame button 
        lviewQueryi = Label(fviewQueryi, anchor = 'w', text='represents graphically a query series')
        lviewQueryi.pack(side="left", fill="y", padx=5, pady=5)        
        self.buttonQuery4 = Button(fviewQueryi, text ='OK', width=5, command = self.view.viewQuery)
        self.buttonQuery4.pack(side="right", anchor='e', padx=5, pady=5)

        fviewQuery0.pack(side="top", expand=YES, anchor='w')
        fviewQuery1.pack(side="top", expand=YES, anchor='w')
        fviewQuery2.pack(side="top", expand=YES, anchor='w')
        fviewQueryi.pack(side="top", expand=YES, anchor="e")

        fviewQuery.pack(side="top", fill="x", expand=False, padx=5, pady=5)

        fviewModel = LabelFrame(f32, text='view model')
                
        fviewModeli = Frame(fviewModel)

        # frame button 
        lviewModeli = Label(fviewModeli, anchor = 'w', text='represents graphically model quality')
        lviewModeli.pack(side="left", fill="y", padx=5, pady=5)        
        self.buttonModel4 = Button(fviewModeli, text ='OK', width=5, command = self.view.viewModel)
        self.buttonModel4.pack(side="right", anchor='e', padx=5, pady=5)
        
        fviewModeli.pack(side="top", expand=YES, anchor="e")
        fviewModel.pack(side="top", fill="x", expand=False, padx=5, pady=5)


        # TABS packing
        f32.pack(side="top",fill='x', expand=False)
        
        t1.pack(side="left", fill="both",expand=True)
        i1.pack(side="left", fill="both",expand=True)
        i2.pack(side="right", fill="both",expand=False)
        
        # Start queue listener
        self.periodicCall()


    def selectQueryFile(self):
        selection=tkFileDialog.askopenfilename(parent=root, filetypes=( ("Series","*.sdf"), ("All files", "*.*")) )
        if selection:
            self.eviewQuery1.delete(0, END)
            self.eviewQuery1.insert(0,selection)

    def selectTrainingFile(self):
        selection=tkFileDialog.askopenfilename(parent=root, filetypes=( ("Series","*.sdf"), ("All files", "*.*")) )
        if selection:
            self.buildSeries.delete(0, END)
            self.buildSeries.insert(0,selection)
         
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
        

    def updateGUI (self,newVersions=False):
        self.models.chargeData()

        if newVersions:
            comboVersions = ()
            for i in range(self.models.maxver):
                comboVersions=comboVersions+(str(i),)
            self.eview2['values'] = comboVersions
        

    ##################
    ### DIRECT TASKS
    ##################

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
        self.updateGUI(True)        

    '''
    Presents information about the model defined by the endpoint
    If no version is specified (general endpoint), shows all the model versions.
    '''    
    def seeDetails(self):

        name = self.models.selEndpoint()
        version = self.models.selVersion()

        # Obtain information about the model

        output=commands.getoutput(wkd+'/manage.py -e '+ name +' -v' + version  +' --info=long')
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
        
        text = Text(winDetails, wrap=WORD, font=(self.myfont,10), yscrollcommand=scrollbar.set)
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
                        endpointDir = wkd + '/' + self.models.selEndpoint() +'/version%0.4d'%int(self.models.selVersion())
                        files = glob.glob(endpointDir+"/pls-*.png")
                        files.sort()
                        self.win=visualizewindow()
                        self.win.viewMultiple(files)

                    self.updateGUI()
                    tkMessageBox.showinfo("Info Message", msg)
                    
                if 'View completed' in msg:
                    self.button4.configure(state='normal')
                    self.buttonQuery4.configure(state='normal')
                    self.win=visualizewindow()
                    self.win.viewSingle(msg[15:])

                if 'finished' in msg:
                    self.updateGUI(True)

                if msg.startswith('Please'):                         # so far, only in View
                    self.button4.configure(state='normal') # view OK
                    self.buttonQuery4.configure(state='normal') # view OK
                    tkMessageBox.showinfo("Info Message", msg)

                if msg.endswith('failed'):                           # so far, only in View
                    self.button4.configure(state='normal') # view OK
                    self.buttonQuery4.configure(state='normal') # view OK
                    tkMessageBox.showinfo("Info Message", msg)

            except Queue.Empty:
                pass

        self.master.after(500, self.periodicCall) # re-call after 500ms
 

if __name__ == "__main__":

    root = Tk()
    root.wm_title("etoxlab GUI (beta 0.82)")    

    app = etoxlab(root)
    root.mainloop()
