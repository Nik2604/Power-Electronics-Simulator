import numpy as np
import tkinter as tk
import matplotlib.pyplot as plt
import math



from pylab import *
from scipy.optimize import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure



np.seterr(divide='ignore', invalid='ignore')

win =tk.Tk()
win.title("Simulation of Power Converters")
win.minsize(600, 570)

a_title= tk.Label(win,text="SIMULATION OF POWER CONVERTERS",fg = "blue",bg = "yellow",font = "Helvetica 16 bold italic")

a_title.grid(column=0, row=0, columnspan=3)

mighty = tk.LabelFrame(win, text='Simulation Parameters')

mighty.grid(column=0, row=1, sticky=tk.W+tk.N, padx=10, pady=10)

contypes = tk.LabelFrame(win, text='Types of Converters (Select One, Simulate and Plot)')

contypes.grid(column=1, row=1, sticky=tk.W+tk.N, padx=10, pady=10)

res = tk.LabelFrame(win, text='Plot Area')

res.grid(column=2, row=1,  rowspan=3)

loadtype = tk.LabelFrame(win, text='Type of Load (Select One)')

#loadtype.grid(column=0, row=2, sticky=tk.W, padx=10, pady=10)

loadtype.place(x=10, y=250, anchor="w")

noofphase = tk.LabelFrame(win, text='Number of Phases (Select One)')

#noofphase.grid(column=1, row=2, sticky=tk.W, padx=10, pady=10)

noofphase.place(x=10, y=350, anchor="w")


a_input = tk.Label(mighty, text="RMS Value of Sine Wave/DC Voltage",font = "Verdana 10 bold",fg = "red")

a_input.grid(column=0, row=0, sticky=tk.W)

peamp = tk.Entry(mighty,width=10)

peamp.insert(0,"230")

peamp.grid(column=1, row=0)


a_supply = tk.Label(mighty, text="Supply Frequency in Hz",font = "Verdana 10 bold",fg = "red",justify=tk.LEFT)

a_supply.grid(column=0, row=1, sticky=tk.W)

freq = tk.Entry(mighty,width=10)

freq.insert(0,"50")

freq.grid(column=1, row=1)



a_time = tk.Label(mighty, text="Simulation Time in secs",font = "Verdana 10 bold",fg = "red")

a_time.grid(column=0, row=2, sticky=tk.W)

simtime = tk.Entry(mighty,width=10)

simtime.insert(0,"0.06")

simtime.grid(column=1, row=2)


a_freq = tk.Label(mighty, text="Switching Frequency in Hz",font = "Verdana 10 bold",fg = "magenta")

a_freq .grid(column=0, row=3, sticky=tk.W)

sf = tk.Entry(mighty,width=10)

sf.insert(0,"1000")

sf.grid(column=1, row=3)



a_mi = tk.Label(mighty, text="Modulation Index/Duty Cycle",font = "Verdana 10 bold",fg = "magenta",justify=tk.LEFT)

a_mi.grid(column=0, row=4, sticky=tk.W)

m = tk.Entry(mighty,width=10)

m.insert(0,"0.8")

m.grid(column=1, row=4)



a_fire = tk.Label(mighty, text="Firing Angle in Degrees",font = "Verdana 10 bold",fg = "magenta")

a_fire.grid(column=0, row=5, sticky=tk.W)

alpha = tk.Entry(mighty,width=10)

alpha.insert(0,"45")

alpha.grid(column=1, row=5)



selected = tk.IntVar()

rad1 = tk.Radiobutton(contypes,text='Half Wave Controlled Rectifier', value=1, variable=selected,font = "Verdana 10 bold",fg = "blue",)

rad2 = tk.Radiobutton(contypes,text='Full Wave Controlled Rectifier', value=2, variable=selected,font = "Verdana 10 bold",fg = "blue")

rad3 = tk.Radiobutton(contypes,text='Pulse Width Modulation Inverter', value=3, variable=selected,font = "Verdana 10 bold",fg = "blue")

rad4 = tk.Radiobutton(contypes,text='AC Voltage Regulator', value=4, variable=selected,font = "Verdana 10 bold",fg = "blue")

rad5 = tk.Radiobutton(contypes,text='DC Chopper', value=5, variable=selected,font = "Verdana 10 bold",fg = "blue")

rad1.grid(column=0, row=1, sticky=tk.W)

rad2.grid(column=0, row=2, sticky=tk.W)

rad3.grid(column=0, row=3, sticky=tk.W)

rad4.grid(column=0, row=4, sticky=tk.W)

rad5.grid(column=0, row=5, sticky=tk.W)

selected.set(1)

def sel():
    global l
    l = int(loadchoice.get())
    if l==1:
        lvalue.config(state='disabled')
    else:
        lvalue.config(state='normal')
        

loadchoice = tk.IntVar()

load1 = tk.Radiobutton(loadtype,text='R Load',  value=1, variable=loadchoice,font = "Verdana 10 bold",fg = "blue",command=sel)

load2 = tk.Radiobutton(loadtype,text='RL Load',  value=2, variable=loadchoice,font = "Verdana 10 bold",fg = "blue",command=sel)

load1.grid(column=0, row=1, sticky=tk.W)

load2.grid(column=0, row=2, sticky=tk.W)

loadchoice.set(2)


    

phasechoice = tk.IntVar()

phase1 = tk.Radiobutton(noofphase,text='Single Phase',  value=1, variable=phasechoice,font = "Verdana 10 bold",fg = "blue")

phase3 = tk.Radiobutton(noofphase,text='Three Phase',  value=2, variable=phasechoice,font = "Verdana 10 bold",fg = "blue")

phase1.grid(column=0, row=1, sticky=tk.W)

phase3.grid(column=0, row=2, sticky=tk.W)

phasechoice.set(1)


r_input = tk.Label(loadtype, text="R Value in Ohms",font = "Verdana 10 bold",fg = "red")

r_input.grid(column=1, row=1, sticky=tk.W)

rvalue = tk.Entry(loadtype,width=10)

rvalue.insert(0,"100")

rvalue.grid(column=3, row=1)

l_input = tk.Label(loadtype, text="L Value (H) ",font = "Verdana 10 bold",fg = "red")

l_input.grid(column=1, row=2, sticky=tk.W)

lvalue = tk.Entry(loadtype,width=10)

lvalue.insert(0,"0.1")

lvalue.grid(column=3, row=2)

def func(b):
    k = int(loadchoice.get())   
    a1=math.radians(int(alpha.get()))
    fr=int(freq.get())
    L1=float(lvalue.get())
    R1=float(rvalue.get())
    if k==1:
        langle = (math.atan(2*pi*fr*0/R1))
    else:
        langle = (math.atan(2*pi*fr*L1/R1))
            
    r1=(a1-b)/math.tan(langle)
    return math.sin(langle-b)-(math.sin(langle-a1)*math.exp(r1))




def simulate():
    global s, l, time, amplitude1, amplitude2, amplitude3, amplitude4, amplitude5, lcurrent
    s= int(selected.get())
    l = int(loadchoice.get()) 
    tt = float(simtime.get())
    V= float(peamp.get())
    f1=int(freq.get())
    a=int(alpha.get())
    L=float(lvalue.get())
    R=float(rvalue.get())
    m1=float(m.get())
    f2 = int(sf.get())
    time        = np.arange(0, tt, 0.00001);
    g=len(time)
    beta = (fsolve(func, 3.14))
    print(math.degrees(beta))
    if l==1:
        langle = (math.atan(2*pi*f1*0/R))
        L=0
    else:
        langle = (math.atan(2*pi*f1*L/R))

    Z = sqrt((R*R)+(2*pi*f1*2*pi*f1*L*L))
    ad = radians(a)
    x1 = sin((2*pi*f1*time)-langle)
    x2 = sin(langle-ad)
    e1 = -((2*pi*f1*time)-ad)/tan(langle)
    x3 = exp(e1)
    lcurrent = (V/Z)*((x2*x3)+x1)
    
    if s==1:
        P=2/f1
        amplitude1=(360*(P-abs((time*2%P)-P)))/P
        amplitude2= V*sqrt(2)*(sin(2*pi*f1*time))
        amplitude3=V*sqrt(2)*(sin(2*pi*f1*time))
        i=0
        while i < g:
            if amplitude1[i] >= a and amplitude1[i] <=  math.degrees(beta):
                 amplitude3[i]=amplitude3[i]
                 lcurrent[i] = lcurrent[i]
            else:
                 amplitude3[i]=0
                 lcurrent[i]=0
            i +=1
            
        tk.messagebox.showinfo("Simulation", "Simulation Completed")
        
    elif s==2:
        beta1 = math.degrees(beta)
        beta2 = beta1-180
        if beta1 >= (180 + a):
            beta2 = a
        else:
            beta2 = beta2
            
        print(beta2)
        P=2/f1
        amplitude1=(360*(P-abs((time*2%P)-P)))/P
        amplitude2= V*sqrt(2)*sin(2*pi*f1*time)
        amplitude3=V*sqrt(2)*(sin(2*pi*f1*time))
        i=0
        while i < g:
            if amplitude1[i] >= 0 and amplitude1[i] <=  beta2:
                amplitude3[i]=-amplitude3[i]
                lcurrent[i] = -lcurrent[i]
            elif amplitude1[i] >= a and amplitude1[i] <=  180+beta2:
                 amplitude3[i]=amplitude3[i]
                 lcurrent[i] = lcurrent[i]
            elif amplitude1[i] >= (180+a) and amplitude1[i] <=  360+beta2:
                 amplitude3[i]=-amplitude3[i]
                 lcurrent[i] =-lcurrent[i]
            else:
                 amplitude3[i]=0
                 lcurrent[i] = 0
            i +=1
            
        tk.messagebox.showinfo("Simulation", "Simulation Completed")
        
    elif s==3:
        P=1/f2
        amplitude1=(1*(P-abs((time*2%P)-P)))/P
        amplitude2   = m1*sin(2*pi*f1*time)
        amplitude3   = m1*sin(2*pi*f1*time-pi)
        amplitude4   = sin(2*pi*f1*time-pi)
        amplitude5   = sin(2*pi*f1*time-pi)
        i=0
        while i < g:
             if amplitude1[i] > amplitude2[i]:
                 amplitude4[i]=0
             else:
                 amplitude4[i]=1
             if amplitude1[i] > amplitude3[i]:
                amplitude5[i]=0
             else:
                amplitude5[i]=1
             i +=1
        tk.messagebox.showinfo("Simulation", "Simulation Completed")
    elif s==4:
        beta1 = math.degrees(beta)
        beta2 = beta1-180
        if beta1 >= (180 + a):
            beta2 = a
        else:
            beta2 = beta2
            
        print(math.degrees(beta))
        print(beta2)
        P=2/f1
        amplitude1=(360*(P-abs((time*2%P)-P)))/P
        amplitude2= V*sqrt(2)*sin(2*pi*f1*time)
        amplitude3=V*sqrt(2)*(sin(2*pi*f1*time))
        i=0
        while i < g:
            if amplitude1[i] >= 0 and amplitude1[i] <=  beta2:
                amplitude3[i]=amplitude3[i]
                lcurrent[i] = lcurrent[i]
            elif amplitude1[i] >= a and amplitude1[i] <=  180+beta2:
                 amplitude3[i]=amplitude3[i]
                 lcurrent[i] = lcurrent[i]
            elif amplitude1[i] >= (180+a) and amplitude1[i] <=  360+beta2:
                 amplitude3[i]=amplitude3[i]
                 lcurrent[i] =lcurrent[i]
            else:
                 amplitude3[i]=0
                 lcurrent[i] = 0
            i +=1
            
        tk.messagebox.showinfo("Simulation", "Simulation Completed")
        
    elif s==5:
        P=1/f2
        amplitude1=(1*(P-abs((time*2%P)-P)))/P
        amplitude2   = m1*sin(2*pi*f1*time)
        amplitude3   = m1*sin(2*pi*f1*time-pi)
        amplitude4=sin(2*pi*f1*time-pi)
        amplitude5=sin(2*pi*f1*time-pi)
        i=0
        while i < g:
             if amplitude1[i] > m1:
                 amplitude4[i]=0
             else:
                 amplitude4[i]=V
             
             i +=1
        tk.messagebox.showinfo("Simulation", "Simulation Completed")
   
def plt1():
    if s==1:
        
        fig = Figure(figsize=(6, 5), dpi=100)
        a1 = fig.add_subplot(211)
        a1.plot(time, amplitude2,label="INPUT VOLTAGE")
        a1.plot(time , amplitude3, label="LOAD VOLTAGE")
        a1.set_title('SINGLE PHASE CONTROLLED RECTIFIER-HALF CONTROLLER ')
        a1.set_ylabel('Volts')
        a1.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        a1.grid(b=True, which='major', color='#666666', linestyle=':')
        
        a3 = fig.add_subplot(212)
        #a3.plot(time, amplitude1,label="ANGLE")
        a3.plot(time, lcurrent,label="Load Current")
        a3.set_xlabel('Time in ms')
        a3.set_ylabel('Current in Amps')
        a3.grid(b=True, which='major', color='#666666', linestyle=':')
        
        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig, master=win)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().grid(row=1,column=2,rowspan = 3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5,column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas,toolbarFrame)
        toolbar.update()      
    elif s==2:
        fig1 = Figure(figsize=(6, 5), dpi=100)
        #t = np.arange(0, 3, .01)
        a1 = fig1.add_subplot(211)
        a1.plot(time, amplitude2,label="INPUT VOLTAGE")
        a1.plot(time , amplitude3, label="LOAD VOLTAGE")
        a1.set_title('SINGLE PHASE CONTROLLED RECTIFIER-FULL CONTROLLER ')
        a1.set_ylabel('Volts')
        a1.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        a1.grid(b=True, which='major', color='#666666', linestyle=':')
        

        
        a3 = fig1.add_subplot(212)
        a3.plot(time, lcurrent,label="Load Current")
        a3.set_xlabel('Time in ms')
        a3.set_ylabel('Current in Amps')
        a3.grid(b=True, which='major', color='#666666', linestyle=':')
        
        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig1, master=win)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().grid(row=1,column=2,rowspan = 3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5,column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas,toolbarFrame)
        toolbar.update()
    elif s==3:
        fig2 = Figure(figsize=(6, 5), dpi=100)
        a1 = fig2.add_subplot(311)
        a1.plot(time, amplitude1,label="Carrier")
        a1.plot(time, amplitude2,label="Reference 1")
        a1.plot(time, amplitude3,label="Reference 2")
        a1.set_title('PWM GENERATION FOR SINGLE PHASE INVERTER')
        a1.set_ylabel('Volts')
        a1.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        a1.grid(b=True, which='major', color='#666666', linestyle=':')
        

        
        a3 = fig2.add_subplot(312)
        a3.plot(time , amplitude4, label="PWM FOR SWITCHES S1 AND S2")
        a3.set_ylabel('Volts')
        a3.grid(b=True, which='major', color='#666666', linestyle=':')
        a3.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)

        a4=fig2.add_subplot(313)
        a4.plot(time , amplitude5, label="PWM FOR SWITCHES S3 AND S4")
        a4.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        a4.grid(b=True, which='major', color='#666666', linestyle=':')
        a4.set_ylabel('Volts')
        a4.set_xlabel('Time in ms')
        
        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig2, master=win)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().grid(row=1,column=2,rowspan = 3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5,column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas,toolbarFrame)
        toolbar.update()
    elif s==4:
        fig3 = Figure(figsize=(6, 5), dpi=100)
        #t = np.arange(0, 3, .01)
        a1 = fig3.add_subplot(211)
        a1.plot(time, amplitude2,label="INPUT VOLTAGE")
        a1.plot(time , amplitude3, label="LOAD VOLTAGE")
        a1.set_title('SINGLE PHASE AC VOLTAGE REGULATOR')
        a1.set_ylabel('Volts')
        a1.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        a1.grid(b=True, which='major', color='#666666', linestyle=':')
        

        
        a3 = fig3.add_subplot(212)
        a3.plot(time, lcurrent,label="Load Current")
        a3.set_xlabel('Time in ms')
        a3.set_ylabel('Current in Amps')
        a3.grid(b=True, which='major', color='#666666', linestyle=':')
        
        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig3, master=win)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().grid(row=1,column=2,rowspan = 3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5,column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas,toolbarFrame)
        toolbar.update()
    elif s==5:
        fig4 = Figure(figsize=(6, 5), dpi=100)
        a1 = fig4.add_subplot(211)
        a1.plot(time, amplitude4,label="Load Voltage")
        a1.set_title('LOAD VOLTAGE DC-DC BUCK CONVERTER')
        a1.set_ylabel('Volts')
        a1.grid(b=True, which='major', color='#666666', linestyle=':')
        

        
        a3 = fig4.add_subplot(212)
        a3.plot(time , amplitude1, label="CARRIER SIGNAL")
        a3.set_ylabel('Volts')
        a3.grid(b=True, which='major', color='#666666', linestyle=':')
        a3.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)

       
        
        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig4, master=win)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().grid(row=1,column=2,rowspan = 3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5,column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas,toolbarFrame)
        toolbar.update()

btn = tk.Button(win, text="SIMULATE", command=simulate,font = "Verdana 10 bold")

#btn.grid(column=0, row=3,padx=25, pady=15, sticky=tk.W+tk.N)

btn.place(x=120, y=540, anchor="w")

btn1 = tk.Button(win, text="PLOT", command=plt1,font = "Verdana 10 bold")

#btn1.grid(column=1, row=3,padx=25, pady=15, sticky=tk.W+tk.N)

btn1.place(x=490, y=540, anchor="w")

win.mainloop()
