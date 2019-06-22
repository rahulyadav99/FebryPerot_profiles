
#!/usr/bin/env python 
from sys import version_info
if version_info[0] == 2:
    import Tkinter as Tk
else:
    import tkinter as Tk
    
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from math import pi
from scipy.signal import fftconvolve

######################################################################

class fpgui:
    #
    # Compute and display FP profiles 
    #
    #
    # plot window sizes 
    window_width = 8
    window_height= 4

    # entry field width
    fwidth=4

    ######################################################################

    def __init__(self, master):

        # load data only once--

        self.al, self. at = self.get_atmostrans()
               
        self.master = master
        self.master.title("Febrey-Perot profiles")

        self.init_parameters()
        rowcounter = self.init_buttons()
        self.init_plot(rowcounter)

        self.redraw()

    ############################a##########################################
    def init_scale(self):
        rowcounter = 0
        w = Tk.Scale(self.master, from_=0, to=1, orient=HORIZONTAL)
        w.pack()
        rowcounter+=1


    def init_buttons(self):
 
        #start adding buttons and text
        rowcounter = 0

        #window label
        input_label = Tk.Label(self.master, text='Set FP parameters:')
        input_label.grid(row=rowcounter, column=0, columnspan=4,sticky='n')
        rowcounter+=1
        
        #w = Tk.Scale(self.master, from_=0, to=1, resolution = 0.1,orient = Tk.HORIZONTAL,variable = self.R)
        
        #w.bind('<Return>',self.redraw_from_event)
        #w.grid(row=rowcounter, column=1)
        #rowcounter+=1   


        # reflectivity
        Tk.Label(self.master, text='Reflectivity [0-1]').grid(row=rowcounter, column=0,sticky='w')
        R_entry = Tk.Entry(self.master, width= self.fwidth, textvariable = self.R)
        R_entry.bind('<Return>',self.redraw_from_event)
        R_entry.grid(row=rowcounter, column=1)
        rowcounter+=1


        # thickness of crystal
        Tk.Label(self.master, text='Seperation between plates [mm]').grid(row=rowcounter, column=0,sticky='w')
        T_entry = Tk.Entry(self.master, width= self.fwidth, textvariable = self.T)
        T_entry.bind('<Return>',self.redraw_from_event)
        T_entry.grid(row=rowcounter, column=1)
        rowcounter+=1
        
        # referactive index of crystal
        Tk.Label(self.master, text='Refractive index').grid(row=rowcounter, column=0,sticky='w')
        T_entry = Tk.Entry(self.master, width= self.fwidth, textvariable = self.SN)
        T_entry.bind('<Return>',self.redraw_from_event)
        T_entry.grid(row=rowcounter, column=1)
        rowcounter+=1

        # lmin,
        Tk.Label(self.master, text='lambda_min/max [A]').grid(row=rowcounter, column=0,sticky='w')
        lmin_entry = Tk.Entry(self.master, width=self.fwidth, textvariable = self.lmin)
        lmin_entry.bind('<Return>',self.redraw_from_event)
        lmin_entry.grid(row=rowcounter, column=1)

        #lmax
        Tk.Label(self.master, text='--').grid(row=rowcounter, column=2)
        lmax_entry = Tk.Entry(self.master, width=self.fwidth, textvariable = self.lmax)
        lmax_entry.bind('<Return>',self.redraw_from_event)
        lmax_entry.grid(row=rowcounter, column=3)
        rowcounter+=1
        
        # quit button
        quit_button = Tk.Button(self.master, text="Quit", command=self.quit)
        quit_button.grid(row=rowcounter, column=0, columnspan=4,sticky='ew')

        return rowcounter

    ######################################################################

    def init_plot(self, rowspan):

        # create Figure for the plot 
        self.fig = Figure(figsize=(self.window_width, self.window_height), dpi=100)
        self.axes=[]

        nrow = 1
        ncol = 1
        for i in range(nrow*ncol):
            self.axes.append(self.fig.add_subplot(nrow, ncol, i+1))

        # link figure to gui
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        self.canvas.get_tk_widget().grid(row=0,column=4,rowspan=rowspan)
        self.canvas.draw()

    ######################################################################
        
    def init_parameters(self):
        R_init = 0.9
        self.R = Tk.DoubleVar()
        self.R.set(R_init)

        T_init = 0.26
        self.T  = Tk.DoubleVar()
        self.T.set(T_init)

        SN_init = 2.28
        self.SN = Tk.DoubleVar()
        self.SN.set(SN_init)

        lmin_init=6160
        self.lmin = Tk.DoubleVar()
        self.lmin.set(lmin_init)

        lmax_init=6180
        self.lmax = Tk.DoubleVar()
        self.lmax.set(lmax_init)

        v_init = 0.7
        self.v = Tk.DoubleVar()
        self.v.set(v_init)
 
    ######################################################################

    def quit(self):
        self.master.destroy()

    ######################################################################


    ######################################################################

    def redraw_from_event(self, event):

        self.redraw()

    ######################################################################

    def redraw(self):

        # read the input parameters from buttons and fields
   
        R  = self.R.get()               # spectral resolution
        print('RR',R) 
        d = self.T.get()                # thickness of material/ spacing between plates
        n0 = self.SN.get()              # refractive index
        lmin = self.lmin.get()          # wavelength range
        lmax = self.lmax.get()      
        telescopetrans = self.T.get()   # telescope transmission

        #v = self.v.get() * self.km_to_m # m/s


   
        #d1 = 0.226    # thickness (mm)
        #R = .87       # reflectivity

        theta = 0.0   # angle of incidence
        vt = 0
        #n0 = 2.28496

        lamb = np.linspace(lmin,lmax,5001)
        trans, shift = self.fp_profile( lamb, d, R, n0, vt, theta)

        ll = self.al
        tt =  self.at
        self.plot(lamb,trans,ll,tt,lmin,lmax)


    ######################################################################

    def plot(self, xx,yy,ll,tt,lmin,lmax):


        for ax in self.axes:
            ax.clear()


        ll = self.al
        tt = self.at
        ax = self.axes[0]
        ax.plot(xx-6173.34, yy )

        ax.plot(ll-6173.34,tt,color='gray',alpha=0.3)

        #ax.legend(loc='best')
        #ax.get_xaxis().set_ticklabels([])
        ax.set_ylabel('Intensity')
        ax.set_xlabel('$\lambda$-6173 [$\AA$]')
        ax.set_xlim(lmin-6173,lmax-6173)
        ax.set_title('spectrum')
        #ax.grid(True)

        self.fig.tight_layout()


        self.canvas.draw()
      
    ######################################################################  

    def fp_profile(self, lamb, d, R, n0, V, theta):
        print(d, R, n0, V, theta)
        lam = lamb
        lam0 = 6173.34
        F = (4.*R)/(1.0-R)**2.      
        finess = np.pi*np.sqrt(F)*0.5

        r13 = 5.2                   #pm/v
        E = V/(d)                     #electric field
        n = n0+0.5*n0**3.*E*r13*1e-6   #change in n0

        tmp = 2*n*d*np.cos(theta)          #phase

        fsr = lam0**2/(tmp)                #free spectral range (FSR)    
        fsr = fsr*1.e-7                   #FSR (A)

        fwhm = fsr/finess                 #fullwidth half maxima

        delt = (2.*np.pi)*tmp

        delt = (delt*1.e7)/(lam)

        denom = (1.+ F*(np.sin(delt/2.))**2.)
        trans = 1./denom
        print('FSR, FWHM, Finesse :', fsr,fwhm,finess)

        shift=0.5*(n0**2)*E*r13*1e-6*lam0
        return trans, shift


######################################################################
    def get_atmostrans(self):
        #
        #read Bass2000-solar spectrum across 6173
        # 

        a = np.loadtxt('bass_6173.dat')
        ll = a[:,0]
        yin = a[:,1]
        TT = yin/np.max(yin)
        return ll, TT

#######################################################################
        
root = Tk.Tk()
my_gui = fpgui(root)
root.mainloop()
