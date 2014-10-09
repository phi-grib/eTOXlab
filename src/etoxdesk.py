from Tkinter import *           # Importing the Tkinter (tool box) library 
import tkMessageBox, tkFileDialog, Tkconstants
import commands, os, subprocess
import thread

def chargeData():
    #Read a file and add all the data to a file    
    file = open('manage.txt', 'r')
    
    for line in file:
        datos.append(line)
        
    #Insert the list on a Listbox    
    for dato in datos:
        listbox.insert(END, dato)
    
    #First item selected    
    listbox.selection_set( first = 0 )
    
    return listbox

def chargeVersion(e):
    
    listboxversion.delete(0, END)
     
    
    d=datos[whichSelected(listbox)]
    name= d.split()[0];
       
    # Obtain information about the model
    output=commands.getoutput('manage.py -e '+name+' --info=short')
    output= output.split("\n")

    for line in output[2:-1]:
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
    output=commands.getoutput('manage.py -e '+name+' -v' +version  +' --info=long')
    
    # Show the collected information in a new window (winDetails)
    winDetails = Tk()
    winDetails.resizable(0.5,0.5)
    win1 = PanedWindow()
    winDetails.wm_title(name+'_version'+version)
    lab= Label(winDetails, text = output,font=("Helvetica", 12),justify=LEFT)
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
    os.system('mkdir -p /home/modeler/workspace/'+name+ ' && cp -rf '+filebut+' /home/modeler/workspace/'+name)
    os.chdir('/home/modeler/workspace/'+name)
 
             
    os.spawnlp(os.P_NOWAIT, 'python /home/modeler/soft/eTOXlab/src/build.py -e '+name+' -f '+filebut+' -v'+version)
    
    # Show a info message
    tkMessageBox.showinfo("Info Message", "ReBuilding model "+name+'with dataset '+filebut)

def openFile():
    filename = tkFileDialog.askopenfilename(parent=root,initialdir='/home/modeler',title="Select a new training set to rebuild the model" , filetypes=[('image files', '.sdf')]) ## filename not filehandle
    return filename

def Visualize():
    lab= Label(root, text = 'See results')
    lab.pack()

def makeWindow():
    global select,scroll, listbox,datos,datos1, listboxversion,root
    datos = []
    datos1 = []
    
    root = Tk()
    root.geometry('600x270+600+100')
    root.wm_title("eTOXLAB Gui")
    m1 = PanedWindow()
    m1.pack(fill=BOTH, expand=1)
    
    # Model details
    listbox = Listbox(root,exportselection=False)
    m1.add(chargeData())
    
    listboxversion = Listbox(root)        
    m2 = PanedWindow(m1, orient=VERTICAL)  
    listbox.bind('<<ListboxSelect>>', chargeVersion) 
    m2.add(listboxversion)
    m1.add(m2)
    
    # Buttons
    button1 = Button(root, text = 'Model details', command = seeDetails)
    button2 = Button(root, text = 'ReBuild', command = reBuild)
    button3 = Button(root, text = 'View charts', command = Visualize)
    top = Label(m2, text="top pane")

    m2.add(button1)
    m2.add(button2)
    m2.add(button3)  
    
    return root
     
root = makeWindow()
root.mainloop() 


