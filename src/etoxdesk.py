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

from Tkinter import *           # Importing the Tkinter (tool box) library 
import Tkconstants

import tkMessageBox
import tkFileDialog

import commands
import os
import subprocess
import shutil

from fabric.colors import white


def chargeData():
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


def chargeVersion(e):
    
    listboxversion.delete(0, END)
     
    d=datos[whichSelected(listbox)]    
    name= d.split()[0];
    
    # Obtain information about the model      
    output=commands.getoutput('/home/modeler/soft/eTOXlab/src/manage.py -e '+name+' --info=short') 
    
    output= output.split("\n")    
   
    for line in output[2:-1]:
        l=line.strip().split()
        
        if (len(l) == 10):

            line= '%-6s'%l[0]+'%-16s'%l[1]+'%-14s'%l[2].replace('mod:',' ')+'%-10s'%l[5]+'%-12s'%l[7].replace(':',' ')+'%-12s'%l[9].replace(':',' ')
            listboxversion.insert(END, line)
    
    listboxversion.selection_set( first = 0 )
   
    return listboxversion


def whichSelected(l):
    return int(l.curselection()[0])


def seeDetails():
    # Capture the name of the selected model
    d=datos[whichSelected(listbox)]
    name= d.split()[0];
    
    version = whichSelected(listboxversion).__str__()    

    # Obtain information about the model
    output=commands.getoutput('/home/modeler/soft/eTOXlab/src/manage.py -e '+name+' -v' +version  +' --info=long')
      
    # Show the collected information in a new window (winDetails)
    winDetails = Tk()
    winDetails.resizable(0.5,0.5)
    win1 = PanedWindow()
    winDetails.wm_title(name+' ver '+version)
    
    lab= Label(winDetails, text = output,font=("Courier New", 12),justify=LEFT)
    lab.pack()
    winDetails.mainloop()

    
def reBuild():
    # Choose the model
    select=datos[whichSelected(listbox)]  
    name= select.split()[0];
    
    # Choose the version
    version = whichSelected(listboxversion).__str__()    

    # Choose new training set
    button_opt = {'fill': Tkconstants.BOTH, 'padx': 5, 'pady': 5}
    filebut = openFile()
    
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

    subprocess.Popen('python /home/modeler/soft/eTOXlab/src/build.py -e '+name+' -f training.sdf -v'+version, stdout=subprocess.PIPE, shell=True)
   
    # Show a info message
    tkMessageBox.showinfo("Info Message", "ReBuilding model "+name+'with dataset '+filebut)


def openFile():
    filename = tkFileDialog.askopenfilename(parent=root,initialdir='/home/modeler',title="Select a new training set to rebuild the model" , filetypes=[('image files', '.sdf')]) ## filename not filehandle
    return filename


def publish():   
    # Choose the model
    select = datos[whichSelected(listbox)]  
    name = select.split()[0];
    
    subprocess.call('/home/modeler/soft/eTOXlab/src/manage.py --publish -e '+name, stdout=subprocess.PIPE, shell=True)  

    chargeVersion(0)   


def Visualize():
    lab= Label(root, text = '')
    lab.pack()


def makeWindow():
    global scroll,listbox,datos,datos1,listboxversion,root
    datos = []
    datos1 = []
    
    root = Tk()
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
    m1.add(chargeData())

    m0.add(m1)
    
    # Version details - Paned m2
    labelListv = Label(root, text="#    MD             mod        mol       R2        Q2", anchor=W, justify=LEFT,font=("Courier New", 10))
    
    #List of Versions
    listboxversion = Listbox(root,width=70, font=("Courier New", 9))        
    m2 = PanedWindow(m0, orient=VERTICAL)  
    listbox.bind('<<ListboxSelect>>', chargeVersion) 
    root.unbind('<Leave>')
    m2.add(labelListv)
    m2.add(listboxversion)

    chargeVersion(0)
        
    # Buttons
    button1 = Button(root, text = 'model details', command = seeDetails,font=("Courier New", 10))
    button2 = Button(root, text = 'rebuild', command = reBuild,font=("Courier New", 10))
    button3 = Button(root, text = 'publish', command = publish,font=("Courier New", 10))
    button4 = Button(root, text = 'view charts', command = Visualize,font=("Courier New", 10))
    
    m2.add(button1)
    m2.add(button2)
    m2.add(button3)  
    m2.add(button4)
    m0.add(m2)
    
    return root


def main ():
    
    root = makeWindow()
    root.mainloop()
    sys.exit(0)

if __name__ == '__main__':
    
    main()
