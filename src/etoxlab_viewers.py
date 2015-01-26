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
from utils import cleanSandbox
import tarfile
from PIL import ImageTk, Image
import glob


#####################################################################################################
### visualizeHelp
#####################################################################################################        

class visualizeHelp (Toplevel):
    def __init__(self):        
        Toplevel.__init__(self)
        self.title ('About..')

    def showAbout (self):
        f = Frame(self)
        msg = Message (f,text="An eTOXlab simple GUI\n\n"+
                    "Ines Martinez and Manuel Pastor (manuel.pastor@upf.edu)\n"+
                    "Copyright 2014, 2015 Manuel Pastor", width=600)
        msg.config(bg='white', justify=CENTER, font=("sans",14))
        msg.pack(fill='x', expand=True)

        if os.path.isfile(wkd+'/logoeTOX.png'):
            self.image = ImageTk.PhotoImage(Image.open(wkd+'/logoeTOX.png'))
            self.logo = Label (f, image=self.image,bg='white' )
            self.logo.pack(fill='x', expand=True)

        ops = Message (f,text="\n\neTOXlab is free software: you can redistribute it and/or modify"+
                    "it under the terms of the GNU General Public License as published by"+
                    "the Free Software Foundation version 3.", width=600)
        ops.config(bg='white', justify=LEFT, font=("sans",10))
        ops.pack(fill='x', expand=True)
        f.pack()

#####################################################################################################
### visualizeDetails
#####################################################################################################   

class visualizeDetails (Toplevel):
    def __init__(self):        
        Toplevel.__init__(self)

    def showDetails (self, model, output):
        self.title (model)

        scrollbar = Scrollbar(self,orient=VERTICAL)
        scrollbar.pack(side=RIGHT, fill=Y)
        
        text = Text(self, wrap=WORD, font=('Courier New',10), yscrollcommand=scrollbar.set)
        text.insert(INSERT, output)
        text.config(state=DISABLED)
        text.pack(side="top", fill="both", expand=True)
        
        scrollbar.config(command=text.yview)      
        

#####################################################################################################
### visualize graphics 
##################################################################################################### 

'''
Creates a new window that displays one or more plots given as a list of png files
'''
class visualizewindow(Toplevel):
    
    def __init__(self, vtitle='graphic viewer'):        
        Toplevel.__init__(self)
        self.title (vtitle)

    def viewFiles (self, fnames):
        #if not fnames : return
        
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
            self.note_view.add(f,text=os.path.splitext(os.path.basename(t))[0])
            ttk.Label(f,image=self.i[-1]).pack()

        if not len(self.i):
            self.destroy()
            return
        
        self.note_view.pack()

#####################################################################################################
### visualizePrediction
#####################################################################################################        
        
'''
Creates a new window that displays one or more plots given as a list of png files
'''
class visualizePrediction (Toplevel):
    
    def __init__(self):   
        Toplevel.__init__(self)
        self.title ('Prediction results')

        f0 = Frame (self)
        scrollbar_tree = ttk.Scrollbar(f0)
        self.tree = ttk.Treeview (f0, columns = ('#','a','b','c'), selectmode='browse',yscrollcommand = scrollbar_tree.set)
        self.tree.column ("#", width=50, anchor='center' )
        self.tree.column ('a', width=120, anchor='e')
        self.tree.column ('b', width=50, anchor='center')
        self.tree.column ('c', width=120, anchor='e')
        self.tree.heading ('#', text='mol#')
        self.tree.heading ('a', text='value')
        self.tree.heading ('b', text='AD')
        self.tree.heading ('c', text='CI')

        scrollbar_tree.pack(side="left", fill=Y)
        scrollbar_tree.config(command = self.tree.yview)        
        
        self.tree.pack(side='top', expand=True, fill='both')
        f0.pack(side="top", expand=True, fill='both')

        
    def show (self, endpoint, version):

        ## check if endpoint+version already exists
        if endpoint+version in self.tree.get_children():
            self.tree.delete(endpoint+version)
        
        self.tree.insert ('','end', endpoint+version, text=endpoint+' ver '+ version, open=True )
        f = open ('/var/tmp/results.txt','r')
        count = 0
        for line in f:
            result = line.split()

            if len(result) < 6:
                continue

            value = 'na'
            AD = 'na'
            CI = 'na'

            if not "failed" in line:
                
                if result[0]!='0':
                    try:
                        v = float(result[1])
                        value = '%10.3f'%v
                    except:
                        value = result[1]

                if result[2]!='0':        
                    AD = result[3]

                if result[4]!='0':
                    try:
                        c = float(result[5])
                        CI = '%10.3f'%c
                    except:
                        CI = result[5]

            self.tree.insert(endpoint+version, 'end', values=(str(count),value,AD,CI), iid=endpoint+version+str(count))
            count+=1
            
        f.close()
