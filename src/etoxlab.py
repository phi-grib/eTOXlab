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
import os
import subprocess
import shutil
import Queue

from threading import Thread
from utils import wkd
from utils import VERSION
from utils import removefile
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
    
    def __init__(self, parent, command, seeds, q):       
        self.model = parent
        self.command = command
        self.seeds = seeds
        self.queue = q
        
        self.dest =''
        
    def processJob(self):
        p = processWorker(self.command, self.seeds, self.queue, self.dest)                                                                       
        p.compute()                      

    def process(self):       
        self.seeds = []        
        self.seeds.append (self.model.selEndpoint())
        self.seeds.append (self.model.selVersion())
                            
        if self.command=='--get=series' or self.command=='--export':
            self.dest=tkFileDialog.askdirectory(initialdir='.',title="Choose a directory...")                    
            
        t = Thread(target=self.processJob)
        t.start()
            
        
class processWorker: 

    def __init__(self, command, seeds, queue, dest):
        self.command = command
        self.seeds = seeds
        self.q = queue
        self.dest = dest

    def compute(self):

        endpoint = self.seeds[0]
        mycommand = [wkd+'/manage.py','-e', endpoint, self.command]
        
        if self.command=='--get=series':
            version = self.seeds[1]
            mycommand.append ('-v')
            mycommand.append (version)
            os.chdir(self.dest)
            
        elif self.command=='--export':
            os.chdir(self.dest)
        
        try:
            result = subprocess.call (mycommand)
        except:
            self.q.put ('ERROR: Unable to execute manage command')
            return
        
        if result == 1 :
            self.q.put ('ERROR: Failed to perform manage operation')
            return

        if self.command in ['--publish','--remove'] :
            self.q.put ('update')
            
        self.q.put ('Manage command finished')


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

        # Call new thread to visualize the series       
        app.viewButton1.configure(state='disable')

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
        app.viewButton2.configure(state='disable')

        self.vtype       = app.viewTypeComboQuery.get()
        self.refname     = self.seeds[0]
        self.refver      = self.seeds[1]
        self.molecules   = app.eviewQuery1.get()      
        self.background  = (app.viewBackgroundQuery.get() == '1')
        
        t = Thread(target=self.viewJob)
        t.start()

    def viewModel(self):
        #endpointDir = wkd+'/' + app.models.selEndpoint() +'/version%0.4d'%int(app.models.selVersion())
        endpointDir = app.models.selDir()
        files = [endpointDir+'/recalculated.png', endpointDir+'/predicted.png']
        for i in files:
            if not os.path.isfile(i):
                tkMessageBox.showinfo("Info Message", "No visual info available for this model")
                return
                
        self.win=visualizewindow('model: '+app.models.selEndpoint()+' ver '+app.models.selVersion())
        self.win.viewFiles(files)
            

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
            proc = subprocess.Popen(mycommand,stdout=subprocess.PIPE)
        except:
            self.q.put ('ERROR: View process failed')
            return
 
        for line in iter(proc.stdout.readline,''):
            line = line.rstrip()
            if line.startswith('ERROR:'):
                self.q.put (line)
                return

        if proc.wait() == 1 :
            self.q.put ('Unknown error')
            return
        
        if self.vtype=='pca' or self.vtype=='project':
            outname = './pca-scores12.png'
        else:
            outname = './property.png'
            
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
        series = app.buildSeries.get()
        model = app.buildModel.get()

        # If 'series' starts with '<series' then copy training.sdf, tdata.pkl, info.pkl to the sandbox
        
        if series.startswith('<series'):
            series = ''
            if version != '0' :
                try:
                    print 'copying'
                    shutil.copy(self.model.selDir()+'/training.sdf',wkd+'/'+name+'/version0000/')
                    shutil.copy(self.model.selDir()+'/tstruct.sdf',wkd+'/'+name+'/version0000/')
                    shutil.copy(self.model.selDir()+'/tdata.pkl',wkd+'/'+name+'/version0000/')
                    shutil.copy(self.model.selDir()+'/info.pkl',wkd+'/'+name+'/version0000/')
                except:
                    print 'error in copy'
                    
        if not model.startswith('<edited'):
            if version != '0' :
                try:
                    shutil.copy(self.model.selDir()+'/imodel.py',wkd+'/'+name+'/version0000/')
                except:
                    print 'error in copy'
            
        
        # Add argument to build list 
        self.seeds = [] 
        self.seeds.append(name)
        self.seeds.append('0')
        self.seeds.append(series)
    
        # Call new thread to build the model       
        app.buildButton.configure(state='disable')
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
        series  = self.seeds[2]

        mycommand = [wkd+'/build.py','-e',name,'-v',version]

        if series and series !='':
            mycommand.append ('-f')
            mycommand.append (series)

        try:
            proc = subprocess.Popen(mycommand,stdout=subprocess.PIPE)
        except:
            self.q.put ('Building process failed')
            return
 
        for line in iter(proc.stdout.readline,''):
            line = line.rstrip()
            
            if line.startswith('ERROR:'):
                self.q.put (line)
                return

            if "Model OK" in line:
                self.q.put('Building completed OK')
                return
        
        if proc.wait() == 1 :
            self.q.put ('Unknown error')
            return

        self.q.put('update')
        
        self.q.put('Building command finished')
    


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
        self.updateFocus()
        
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

    def selDir (self):
        e = self.selEndpoint()
        v = self.selVersion()
        return (wkd + '/' + e +'/version%0.4d'%int(v))
      
    # Charges general information about models.    
    def chargeData(self):
        self.clearTree()
        version = []

        process = subprocess.Popen([wkd+'/manage.py','--info=short'], stdout=subprocess.PIPE)
        output, err = process.communicate()

        output=output.split('------------------------------------------------------------------------------')
        
        for line in output[1:]:
            line=line.split("\n")
            self.datos.append(line)
            name=line[1].split("[")[0]            
            self.insert('', 'end', '%-9s'%(name), text='%-9s'%(name))
            count=0                            
            for x in line[2:-1]:
                if x.startswith ("All requested models"):
                    continue
                version=self.chargeDataDetails(x)
                if len(version)>4:
                    self.insert('%-9s'%(name), 'end', values=(version[0],version[1],version[2],version[3],version[4]), iid='%-9s'%(name)+str(count))
                    # The list of child is unfold in treeView
                    self.item('%-9s'%(name), open=True)                    
                count+=1              

        self.maxver = 1
        for child in self.get_children():
            iver= len(self.get_children(child))
            if iver > self.maxver : self.maxver=iver

    def updateFocus(self):
        # Focus in first element of the TreeView
        self.selection_set(self.get_children()[0:1])
        self.focus(self.get_children()[0:1][0])

    # Charges detailed information about each version of a given model(line).   
    def chargeDataDetails(self,line):       
        y = []

        line=line.replace(':',' ')
        l=line.split()

        y.append('%-5s' %l[0])         # num

        if 'no model info available' in line:
            y.append ('na') #MD
            y.append ('na') #mod
            y.append ('na') #mol
            y.append ('no info available') # quality
            return y
        
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
            y.append('not recognized')        
  
        return y


#####################################################################################################
### VIEWER WINDOW CLASS
#####################################################################################################        
        
'''
Creates a new window that displays one or more plots given as a list of png files
'''
class visualizewindow(Toplevel):
    
    def __init__(self, vtitle='graphic viewer'):        
        Toplevel.__init__(self)
        self.title (vtitle)

    def viewFiles (self, fnames):
        if len(fnames)<2:
            self.viewSingle (fnames[0])
        else:
            self.viewMultiple (fnames)
            
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
        self.skipUpdate = False

        # create a toplevel menu
        menubar = Menu(root)
        filemenu = Menu(menubar)    
        filemenu.add_command(label="Exit", command= quit)
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
        self.models.bind('<<TreeviewSelect>>', self.selectionChanged)
     
        # Main container is a Notebook
        n = ttk.Notebook (i2)
        f1 = Frame(n)
        f2 = Frame(n)
        f3 = Frame(n)        
        
        f1.pack(side="top", fill='x', expand=False)
        f2.pack(side="top", fill='x', expand=False)
        f3.pack(side="top", fill='x', expand=False)
        
        n.add (f1, text='manage')
        n.add (f2, text='build')
        n.add (f3, text='view')
        n.pack (side="top", fill="x", expand=False)
       
        ## MANAGE Frame        
        f12 = Frame(f1)
    
        fnew = LabelFrame(f12, text='new endpoint')
        
        fnewi = Frame(fnew)
        fnew1 = Frame(fnewi)
        fnew2 = Frame(fnewi)
        fnewj = Frame(fnew)

        vcmd1 = (root.register(self.validateEndpoint), '%S', '%P')
        
        Label(fnew1, width = 10, anchor='e', text='name').pack(side="left")       
        self.enew1 = Entry(fnew1, bd =1, validate = 'key', validatecommand = vcmd1 )
        self.enew1.pack(side="left")               # field containing the new endpoint name

        vcmd2 = (root.register(self.validateTag), '%S', '%P')
        
        Label(fnew2, width = 10, anchor='e', text='tag').pack(side="left")
        self.enew2 = Entry(fnew2, bd =1, validate = 'key', validatecommand = vcmd2)
        self.enew2.pack(side="left")               # field containing the new endpoint tag
       

        Label(fnewj, text='creates a new endpoint').pack(side="left", padx=5, pady=5)       
        Button(fnewj, text ='OK', command = self.new, width=5).pack(side="right", padx=5, pady=5)

        fnew1.pack(fill='x')
        fnew2.pack(fill='x')

        fnewi.pack(fill="x" )
        fnewj.pack(fill="x" )        
        fnew.pack(fill="x", padx=5, pady=5)
        
        finfo = LabelFrame(f12, text='get information')
        Label(finfo, text='shows complete model info').pack(side="left", padx=5, pady=5)
        Button(finfo, text ='OK', command = self.seeDetails, width=5).pack(side="right", padx=5, pady=5)        
        finfo.pack(fill='x', padx=5, pady=5)
        
        self.publish=Process(self.models,'--publish', self.seeds, self.q) 

        fpublish = LabelFrame(f12, text='publish model')
        Label(fpublish, text='creates a new model version').pack(side="left",padx=5, pady=5)
        Button(fpublish, text ='OK', command = self.publish.process, width=5).pack(side="right", padx=5, pady=5)
        fpublish.pack(fill='x', padx=5, pady=5)

        self.remove=Process(self.models,'--remove', self.seeds, self.q)
        
        frem = LabelFrame(f12, text='remove model')
        Label(frem, text='removes a model version').pack(side="left",padx=5, pady=5)
        Button(frem, text ='OK', command = self.remove.process, width=5).pack(side="right", padx=5, pady=5)
        frem.pack(fill='x', padx=5, pady=5) 

        self.gseries=Process(self.models,'--get=series', self.seeds, self.q)
        
        fgets = LabelFrame(f12, text='get series')
        Label(fgets, text='saves the training series').pack(side="left", padx=5, pady=5)
        Button(fgets, text ='OK', command = self.gseries.process, width=5).pack(side="right", padx=5, pady=5)
        fgets.pack(fill='x', padx=5, pady=5)

        fgetm = LabelFrame(f12, text='get model')
        Label(fgetm, text='edits the model definition').pack(side="left", padx=5, pady=5)
        Button(fgetm, text ='OK', command = self.modelEdit, width=5).pack(side="right", padx=5, pady=5)
        fgetm.pack(fill='x', padx=5, pady=5)

        self.export=Process(self.models,'--export',self.seeds,self.q)
        
        fexp = LabelFrame(f12, text='export')
        Label(fexp, text='packs whole model tree').pack(side="left",padx=5, pady=5)
        Button(fexp, text ='OK', command = self.export.process, width=5).pack(side="right", padx=5, pady=5)
        fexp.pack(fill="x", padx=5, pady=5)        
       
        fimp = LabelFrame(f12, text='import')
        fimp0 = Frame(fimp)
        fimp1 = Frame(fimp)
        
        Label(fimp0, width = 10, anchor='e', text='import tar').pack(side='left')        
        self.importTar = Entry(fimp0, bd =1)
        self.importTar.pack(side='left')
        Button(fimp0, text ='...', width=2, command = self.selectImportFile).pack(side='left')
        
        Label(fimp1, text='imports packed model tree').pack(side="left", padx=5, pady=5)
        Button(fimp1, text ='OK', command = self.modImport, width=5).pack(side="right", padx=5, pady=5)

        fimp0.pack(fill='x')
        fimp1.pack(fill='x')
        
        fimp.pack(fill='x', padx=5, pady=5)        
        
        f12.pack(fill='x')
        
        ## BUILD Frame        
        self.bmodel=buildmodel(self.models, self.seeds,self.q) 
        
        f22 = Frame(f2)

        fbuild = LabelFrame(f22, text='build model')

        fbuild0 = Frame(fbuild)
        fbuild1 = Frame(fbuild)
        fbuild2 = Frame(fbuild)
        
        Label(fbuild0, width = 10, anchor='e', text='series').pack(side='left')       
        self.buildSeries = Entry(fbuild0, bd =1)
        #self.buildSeries = Entry(fbuild0, bd =1, state='disable')
        self.buildSeries.pack(side='left')      
        Button(fbuild0, text ='...', width=2, command = self.selectTrainingFile).pack(side='left')

        Label(fbuild1, width = 10, anchor='e', text='model').pack(side='left')       
        self.buildModel = Entry(fbuild1, bd =1)
        self.buildModel.pack(side='left')      
        Button(fbuild1, text ='...', width=2, command = self.editModelFile).pack(side='left')

        Label(fbuild2, text='build model  in sandbox with above components').pack(side="left", padx=5, pady=5)
        self.buildButton = Button(fbuild2, text = 'OK', command = self.bmodel.build, width=5)
        self.buildButton.pack(side="right", padx=5, pady=5)

        fbuild0.pack(fill='x')
        fbuild1.pack(fill='x')
        fbuild2.pack(fill='x')
        
        fbuild.pack(fill='x', padx=5, pady=5)

        self.pb = ttk.Progressbar(f22, orient='horizontal', mode='indeterminate', value=0)
        self.pb.pack(fill='x')       

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
        Label (fview0, width = 10, anchor='e', text='type').pack(side='left')
        self.viewTypeCombo = StringVar()
        self.cboCombo = ttk.Combobox(fview0, values=('pca','property','project'), textvariable=self.viewTypeCombo, state='readonly')
        
        self.cboCombo.current(0)
        self.cboCombo.pack()

        # frame 1: entry field for selecting reference endpoint        
        Label (fview1, width = 10, anchor='e', text='refername').pack(side='left')
        self.referEndpointCombo = StringVar ()
        comboValues=("None",) + self.models.get_children()
           
        self.eview1 = ttk.Combobox(fview1, values=comboValues, textvariable=self.referEndpointCombo, state='readonly')
        self.eview1.current(0)
        self.eview1.pack()

        # frame 2: entry field for selecting reference version
        Label(fview2, width = 10, anchor='e', text='refver').pack(side='left')
        self.referVersionCombo = StringVar ()
        
        comboVersions = ()
        for i in range(self.models.maxver):
            comboVersions=comboVersions+(str(i),)  # this is updated by updateGUI method
           
        self.eview2 = ttk.Combobox(fview2, values=comboVersions, textvariable=self.referVersionCombo, state='readonly')
        self.eview2.current(0)
        self.eview2.pack()

        self.eview1.configure(state="disable")
        self.eview2.configure(state="disable")

        # frame 3: check button for showing background
        Label(fview3, width = 10, anchor='e', text='   ').pack(side='left')
        self.viewBackground = StringVar()
        self.checkBackground = ttk.Checkbutton(fview3, text='show background', variable=self.viewBackground, command =lambda: self.updateBack(True))
        
        self.viewBackground.set(0)
        self.checkBackground.pack()
        
        self.cboCombo.bind('<<ComboboxSelected>>', self.updateBack)
        
        # frame button 
        Label(fviewi, text='represents graphically the training series').pack(side="left", padx=5, pady=5)        
        self.viewButton1 = Button(fviewi, text ='OK', width=5, command = self.view.view)
        self.viewButton1.pack(side="right", padx=5, pady=5)

        fview0.pack(anchor='w')

        fview1.pack(anchor='w')
        fview2.pack(anchor='w')
        fview3.pack(anchor='w')
        fviewi.pack(fill='x')

        fview.pack(fill="x", padx=5, pady=5)

        fviewQuery = LabelFrame(f32, text='view query')
                
        fviewQuery0 = Frame(fviewQuery)
        fviewQuery1 = Frame(fviewQuery)
        fviewQuery2 = Frame(fviewQuery)
        fviewQueryi = Frame(fviewQuery)

        # frame 0: combo-box for seletig view type
        Label (fviewQuery0, width = 10, anchor='e', text='type').pack(side='left')
        self.viewTypeComboQuery = StringVar()
        self.cboComboQuery = ttk.Combobox( fviewQuery0, values=('pca','property','project'), textvariable=self.viewTypeComboQuery, state='readonly')
        self.cboComboQuery.current(0)
        self.cboComboQuery.pack(anchor ='w')

        # frame 1: entry field for selecting reference endpoint
        Label(fviewQuery1, width = 10, anchor='e', text='query').pack(side='left')      
        self.eviewQuery1 = Entry(fviewQuery1, bd =1)               # field containing the new endpoint name
        self.eviewQuery1.pack(side="left")
        Button(fviewQuery1, text ='...', width=2, command = self.selectQueryFile).pack(side="right") 

        # frame 2: check button for showing background
        Label (fviewQuery2, width = 10, anchor='e', text='   ').pack(side='left')
        self.viewBackgroundQuery = StringVar()
        self.checkBackgroundQuery = ttk.Checkbutton(fviewQuery2, text='show background', variable=self.viewBackgroundQuery)
        self.viewBackgroundQuery.set(0)
        self.checkBackgroundQuery.pack(anchor='w')

        # frame button 
        Label(fviewQueryi, anchor = 'w', text='represents graphically a query series').pack(side="left", padx=5, pady=5)        
        self.viewButton2 = Button(fviewQueryi, text ='OK', width=5, command = self.view.viewQuery)
        self.viewButton2.pack(side="right", padx=5, pady=5)

        fviewQuery0.pack(anchor='w')
        fviewQuery1.pack(anchor='w')
        fviewQuery2.pack(anchor='w')
        fviewQueryi.pack(fill='x')

        fviewQuery.pack(fill="x", padx=5, pady=5)

        fviewModel = LabelFrame(f32, text='view model')
                
        fviewModeli = Frame(fviewModel)

        # frame button 
        Label(fviewModeli, anchor='w', text='represents graphically model quality').pack(side="left", padx=5, pady=5)        
        Button(fviewModeli, text ='OK', width=5, command = self.view.viewModel).pack(side="right", padx=5, pady=5)
        
        fviewModeli.pack(fill='x')
        fviewModel.pack(fill="x", padx=5, pady=5)


        # TABS packing
        f32.pack(side="top",fill='x', expand=False)
        
        t1.pack(side="left", fill="both",expand=True)
        i1.pack(side="left", fill="both",expand=True)
        i2.pack(side="right", fill="both",expand=False)
        
        # Start queue listener
        self.periodicCall()

    def selectionChanged (self, event):
        # copy this to build series and models
        # print self.models.selEndpoint (), ' ver ', self.models.selVersion()

        if self.skipUpdate:
            self.skipUpdate = False
            return
        
        v = self.models.selVersion()
        self.buildSeries.delete(0, END)
        self.buildSeries.insert(0, '<series ver '+v+'>')

        self.buildModel.delete (0, END)
        self.buildModel.insert (0, '<model ver '+v+'>')

    def updateBack(self, event):
        enableType = (self.viewTypeCombo.get() == 'project' )
        enableBack = (self.viewBackground.get() == '1')

        if enableType or enableBack:
            self.eview1.configure(state="enable")
            self.eview2.configure(state="enable")
        else:
            self.eview1.configure(state="disable")
            self.eview1.current(0)
            self.eview2.configure(state="disable")
            self.eview2.current(0)
            
    def selectImportFile(self):
        selection=tkFileDialog.askopenfilename(parent=root, filetypes=( ("Pack","*.tgz"), ("All files", "*.*")) )
        if selection:
            self.importTar.delete(0, END)
            self.importTar.insert(0,selection)

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

    def editModelFile(self):
        # copy imodel.py of the selected version to the sandbox
        e = self.models.selEndpoint()
        v = self.models.selVersion()
        vdir = wkd + '/' + e + '/version%0.4d/'%int(v)
        zdir = wkd + '/' + e + '/version0000/'

        if (vdir!=zdir):
            try:
                shutil.copy(vdir+'imodel.py', zdir)   # copy imodel.py to the sandbox. This will be the base version for build 
            except:
                tkMessageBox.showerror("Error Message", "Unable to access source imodel.py")
                return

            removefile (zdir+'info.pkl')  # remove imodel.py from the sandbox to force model rebuilding
            
            self.skipUpdate=True
            self.models.chargeData()
            iid = '%-9s'%e + v
            self.models.selection_set((iid,))
            self.models.focus(iid)

        # launch idle with the imodel.py of the sandbox
        try:
            subprocess.Popen(['/usr/bin/idle',zdir+'imodel.py'])
        except:
            tkMessageBox.showerror("Error Message", "Unable to edit imodel.py")
            pass
            
        self.buildModel.delete(0, END)
        self.buildModel.insert(0, '<edited model (save first)>')
            
    def validateEndpoint(self, char, entry_value):
        return (char in 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890_-')
    
    def validateTag(self, char, entry_value):
        return (char in 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890_-/')
         
         
    '''
    Help window
    '''
    def showhelp(self):
        win = Tk()
        win.title('About...')
        msg = Message (win,text="An eTOXlab simple GUI\n\n"+
                    "Ines Martinez and Manuel Pastor (manuel.pastor@upf.edu)\n"+
                    "Copyright 2014, 2015 Manuel Pastor", width=600)
        msg.config(bg='white', justify=CENTER, font=("sans",14))
        msg.pack(fill='x', expand=True)
        ops = Message (win,text="\n\neTOXlab is free software: you can redistribute it and/or modify"+
                    "it under the terms of the GNU General Public License as published by"+
                    "the Free Software Foundation version 3.", width=600)
        ops.config(bg='white', justify=LEFT, font=("sans",10))
        ops.pack(fill='x', expand=True)
        win.mainloop()  
        

    def updateGUI (self,newVersions=False):
        self.models.chargeData()
        self.models.updateFocus()

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
            tkMessageBox.showerror("Error Message", "Please enter the name of the endpoint")
            return

        elif not tag:
            tkMessageBox.showerror("Error Message", "Please enter the name of the tag")
            return

        for line in self.models.get_children():
            labels = line.split()
            
            if endpoint == labels[0]:
                tkMessageBox.showerror("Error Message", "This endpoint already exists!")
                return

        try:
            mycommand = [wkd+'/manage.py', '--new', '-e', endpoint, '-t', tag]
            result = subprocess.call(mycommand)
        except:
            self.q.put ('ERROR: Unable to execute manage command')
            return
        
        if result == 1 :
            self.q.put ('ERROR: Failed to create new endpoint')
            return

        self.enew1.delete(0,END)
        self.enew2.delete(0,END)

        self.updateGUI(True)
        
        tkMessageBox.showinfo("Info Message",'New endpoint created')    


    def modelEdit(self):    
        endpointDir = self.models.selDir()
        modelName = endpointDir + '/imodel.py'

        if not os.path.isfile (modelName):
            self.q.put ('ERROR: no model found')
            return
        
        try:
            subprocess.Popen(['/usr/bin/idle',modelName])
        except:
            self.q.put ('ERROR: no editor found')

    '''
    Presents information about the model defined by the endpoint
    If no version is specified (general endpoint), shows all the model versions.
    '''    
    def seeDetails(self):

        name = self.models.selEndpoint()
        version = self.models.selVersion()

        # Obtain information about the model
        mycommand = [wkd+'/manage.py', '-e', name, '-v', version, '--info=long']
        try:
            process = subprocess.Popen(mycommand, stdout=subprocess.PIPE)
            output, err = process.communicate()
        except:
            tkMessageBox.showerror("Error Message", "Unable to obtain information")
            return
            
        #print output

        outputlist = output.split('\n')
        outputlist = outputlist [1:-2]         
        
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

    def modImport (self):

        importfile = self.importTar.get()

        endpoint = importfile.split('/')[-1]
        endpoint = endpoint [:-4]
        
        if os.path.isdir (wkd+'/'+endpoint):
            tkMessageBox.showerror("Error Message", "This endpoint already exists!")
            return
        
        shutil.copy (importfile,wkd)
        os.chdir(wkd)
       
        tar = tarfile.open(endpoint+'.tgz', 'r:gz')
        tar.extractall()
        tar.close()

        os.remove (endpoint+'.tgz')
        self.importTar.delete(0, END)

        self.updateGUI(True)

        tkMessageBox.showinfo("Info Message",'New endpoint imported')
        

    '''
    Handle all the messages currently in the queue (if any)
    and run events depending of its information
    '''
    def periodicCall(self):
    
        while self.q.qsize():
        
            try:
                msg = self.q.get(0)
                if 'Building' in msg:            
                    self.buildButton.configure(state='normal')
                    self.pb.stop()
                    if 'completed' in msg:
                        endpointDir = wkd + '/' + self.models.selEndpoint() + '/version0000'
                        files = glob.glob(endpointDir+"/pls-*.png")
                        files.sort()
                        self.win=visualizewindow('model: '+self.models.selEndpoint()+' ver 0')
                        self.win.viewFiles(files)
                        
                    #self.updateGUI()
                    iid = '%-9s0'%self.models.selEndpoint()
                    
                    self.models.chargeData()
                    self.models.selection_set((iid,))
                    self.models.focus(iid)
                    
                    tkMessageBox.showinfo("Info Message", msg)
                    
                elif 'View completed' in msg:
                    self.viewButton1.configure(state='normal')
                    self.viewButton2.configure(state='normal')
                    self.win=visualizewindow('series: '+self.models.selEndpoint()+' ver '+self.models.selVersion())
                    files=[msg[15:]]
                    self.win.viewFiles(files)

                elif msg.endswith('failed') or msg.startswith ('ERROR:'):
                    self.viewButton1.configure(state='normal') # view OK
                    self.viewButton2.configure(state='normal') # view OK
                    self.buildButton.configure(state='normal') # building
                    self.pb.stop()
                    tkMessageBox.showerror("Error Message", msg)

                elif 'finished' in msg:
                    tkMessageBox.showinfo("Info Message", msg)
                    
                elif 'update' in msg:
                    self.updateGUI(True)

            except Queue.Empty:
                pass

        self.master.after(500, self.periodicCall) # re-call after 500ms
 

if __name__ == "__main__":

    root = Tk()
    root.title("etoxlab GUI ("+VERSION+")")    

    app = etoxlab(root)
    root.mainloop()
