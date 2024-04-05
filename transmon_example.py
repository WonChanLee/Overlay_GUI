# PyQt5 packages
from PyQt5.QtWidgets import*
from PyQt5.uic import loadUi
from PyQt5 import QtGui

# matplotlib packages
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure

#fitting method
from scipy.optimize import curve_fit

# numpy, os packages
import numpy as np
import random
import os
import pickle

# custom libraries
from replotter_module import replotter

from overlay_helpers import stitch_data_files
from overlay_helpers import overlay_GUI_transmon

# window taskbar icon
import ctypes
myappid = 'mycompany.myproduct.subproduct.version' # arbitrary string
ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

##############################
# Define Gui
##############################
class start_gui(QMainWindow):
    '''
    The start_gui class runs the Overlay_gui.ui file which can be made and modified through Pyqt5 designer. 
    The buttons/toggles/text/graphs on the ui file is designed through designer in the sense that python script and ui file is decoulped each other.
    This means position/size/properties of the object in the ui file does not affect the python script.
    The other important aspect of the Pyqt5 is objects are callable with their object name. 
    There are some preexist functions that we can use. Ex. button.clicked.connect(function).

    The basic design of the overlay gui is almost similar to the original overlay gui with matplotlib, 
    but now fitting feature and some other minor details are added.
    1. Fitting bare transition by clicking 10 reference points by a user.
    2. Load files automatically. 
    Files that includes term "_avoided_" will be stitched at the very end so that avoided crossing feature can be displayed on top.
    3. Fit the transition line using current slider value as an initial guess.
    So use reasonable values for EC, EJ, and EJ_diff that can be found from other methods instead of relying 100% on fitting function. 
    Otherwise, fitting result can converge to multiple local minimums depending on initial guess or fitting method.


    User Manual in steps
    * The only variables that need to be entered through code are volts_at_half, volts_per_flux.
    
    1. run the gui
    2. load pkl files and plot graph
    3. adjust parameters 
    - vmin, vmax - to make more clear plot
    - Ecav, EC, EJ, EJ_diff - make reasonable initial guess for fitting
    4. If Ecav, EC, EJ, EJ_diff need to be changed, write values and click "Set Params".
    5. choose the lines to see
    6. click "Bare Transition Ref". 
    10 reference points on the plot.
    Upto 11.08.23, only bare transition can be fitted.
    It will be saved as coordinates that is used as the fitting data
    7. click "Bare Transition Fit" 
    8. Fitted result will be displayed in the box.
    9. To reset everything, click "Reset".
    This will initialize every variable and plot.
    '''
    
    def __init__(self):
        QMainWindow.__init__(self)
        loadUi("Overlay_Gui.ui",self)
        self.setWindowTitle("Overlay GUI")
        self.setWindowIcon(QtGui.QIcon('hyperbolic.png'))
        self.addToolBar(NavigationToolbar(self.MplWidget.canvas, self))
        

        ### Default params
        self.graph_on = False
        self.freq_on_x_axis = False
        self.bare = False
        self.two_photon = False
        self.cav_assist = False
        self.Ecav = 7.484
        self.EC = 0.3
        self.EJ_sum = 43.5
        self.EJ_diff = 5.5
        self.volts_at_half = 0.985                                      # Not user defined parameter yet. Should enter in the script 
        self.volts_per_flux = 3.66                                      # Not user defined parameter yet. Should enter in the script
        self.offset = self.volts_at_half + self.volts_per_flux/2        # Not user defined parameter yet. Should enter in the script
        self.Ecav_slider.setMaximum(int((float(self.Ecav)+0.5)*1000))
        self.Ecav_slider.setMinimum(int((float(self.Ecav)-0.5)*1000))
        self.EC_slider.setMaximum(int((float(self.EC)+0.2)*1000))
        self.EC_slider.setMinimum(int((float(self.EC)-0.2)*1000))
        self.EJ_sum_slider.setMaximum(int((float(self.EJ_sum)+10)*1000))
        self.EJ_sum_slider.setMinimum(int((float(self.EJ_sum)-10)*1000))
        self.EJ_diff_slider.setMaximum(int((float(self.EJ_diff)+3)*1000))
        self.EJ_diff_slider.setMinimum(int((float(self.EJ_diff)-3)*1000))
        self.Ecav_1.setText(str(round(float(self.Ecav)-0.5, 3)))
        self.Ecav_2.setText(str(round(float(self.Ecav)+0.5, 3)))
        self.EC_1.setText(str(round(float(self.EC)-0.2, 3)))
        self.EC_2.setText(str(round(float(self.EC)+0.2, 3)))
        self.EJ_sum_1.setText(str(round(float(self.EJ_sum)-10, 3)))
        self.EJ_sum_2.setText(str(round(float(self.EJ_sum)+10, 3)))
        self.EJ_diff_1.setText(str(round(float(self.EJ_diff)-3, 3)))
        self.EJ_diff_2.setText(str(round(float(self.EJ_diff)+3, 3)))


        ### MplWidget
        self.guifig = self.MplWidget.canvas.figure     # define figure
        self.guiScans = []                             # define the list of loaded files
        self.files_str = ''                            # define the list of names of loaded files
        self.files_dir = []                            # define files directory in the list with sorted order
        self.lines = []                                # list of lines 
        self.coords = []                               # list of coords that clicked
        self.x_fit = np.zeros(10)                      # array of x coords for fitting
        self.y_fit = np.zeros(10)                      # array of y coords for fitting 
        self.rep_click = 0                             # number for counting click 
        self.line_fit = []                             # list of lines for fitting
        self.flux_range = []                           # define the 
        self.xlim = []                                 # define the xlim of loaded files; default: freq
        self.ylim_flux = []                            # define the ylim_flux of loaded files
        self.ylim_voltage = []                         # define the ylim_voltage of loaded files
        self.vmin = -7                                 # define the vmin
        self.vmax = 2                                  # define the vmax
        self.vmin_slider.setMaximum(self.vmin + 4)
        self.vmin_slider.setMinimum(self.vmin - 4)
        self.vmax_slider.setMaximum(self.vmax + 4)
        self.vmax_slider.setMinimum(self.vmax - 4)
        

        ### load file 
        self.load_file_button.clicked.connect(self.loadfiles)


        ### set transmon params, set slider range
        self.set_params.clicked.connect(self.set_transmon_params)


        ### plot graph
        self.plot_graph.clicked.connect(self.plot_graph_function)


        ### set vmin, vamx and replot graph
        self.vmin_slider.sliderMoved.connect(self.vminmax_slider_change)
        self.vmax_slider.sliderMoved.connect(self.vminmax_slider_change)


        ### check toggles
        self.freq_on_x_axis_chkbox.stateChanged.connect(self.freq_on_x_axis_function)
        self.bare_transition_chkbox.stateChanged.connect(self.bare_transitions_plot)
        self.two_photon_chkbox.stateChanged.connect(self.two_photon_plot)
        self.cav_assist_chkbox.stateChanged.connect(self.cav_assist_plot)


        ### slider change
        self.Ecav_slider.sliderMoved.connect(self.Ecav_slider_change)
        self.EC_slider.sliderMoved.connect(self.EC_slider_change)
        self.EJ_sum_slider.sliderMoved.connect(self.EJ_sum_slider_change)
        self.EJ_diff_slider.sliderMoved.connect(self.EJ_diff_slider_change)


        ### select coords
        self.fit_ref_button.clicked.connect(self.prev_onclick)


        ### Fitting
        self.fit_button.clicked.connect(self.bare_transitions_fit)


        ### Reset
        self.reset_button.clicked.connect(self.reset)





    def loadfiles(self):
        '''
        Opens file browser to select files that user want.
        Opened file directory is saved to the list.
        '''
        self.files_dir_unsorted = QFileDialog.getOpenFileNames(self)
        xlim_search = []                #xlim in multiple files
        ylim_search_voltage = []        #ylim in multiple files
        ylim_search_flux = []
        
        for i in range(len(self.files_dir_unsorted[0])):
            if '_avoided_' in self.files_dir_unsorted[0][i]:
                self.files_dir.append(self.files_dir_unsorted[0][i])
            else:
                self.files_dir.insert(0,self.files_dir_unsorted[0][i])
        
        ### make the set of metadata files that will be included in the GUI
        for i in range(len(self.files_dir)):
            pickledict = pickle.load(open(self.files_dir[i], "rb" ) )
            #xlims: specdata xaxis
            #ylims: specdata yaxis
            xlims = [np.min(pickledict['Data']['specdata']['xaxis']), np.max(pickledict['Data']['specdata']['xaxis'])]
            xlim_search.append(xlims[0])
            xlim_search.append(xlims[1])
            ylims = [np.min(pickledict['Data']['voltages']), np.max(pickledict['Data']['voltages'])]
            ylim_search_voltage.append(ylims[0])
            ylim_search_voltage.append(ylims[1])
            ylim_search_flux.append(round((ylims[0]-self.offset)/self.volts_per_flux, 6))
            ylim_search_flux.append(round((ylims[1]-self.offset)/self.volts_per_flux, 6))
            xSize = 6.5
            ySize = 4.0
            cmap = 'YlGnBu'

            data = replotter(self.files_dir[i],
                                vmin = self.vmin,
                                vmax = self.vmax,
                                cmap = cmap,
                                xlims = xlims,
                                ylims = ylims,
                                fig_xSize = xSize,
                                fig_ySize = ySize,
                                smoothe_data = False,
                                verbose = False,
                                GHzVHz_correction = True)
            self.guiScans.append(data)
            self.files_str += str(self.files_dir[i] + '\n')

        self.load_file_txtbox.setText(self.files_str)
        self.xlim = [np.min(xlim_search), np.max(xlim_search)]
        self.ylim_flux = [np.min(ylim_search_flux), np.max(ylim_search_flux)]
        self.ylim_voltage = [np.min(ylim_search_voltage), np.max(ylim_search_voltage)]
        print(self.xlim)
        print(self.ylim_flux)
        print(self.ylim_voltage)
        return


    def set_transmon_params(self):
        ### read written values 
        self.Ecav = self.Ecav_ref_txt.toPlainText()
        self.EC = self.EC_ref_txt.toPlainText()
        self.EJ_sum = self.EJ_sum_ref_txt.toPlainText()
        self.EJ_diff = self.EJ_diff_ref_txt.toPlainText()

        ### set range
        self.Ecav_1.setText(str(round(float(self.Ecav)-0.5, 3)))
        self.Ecav_2.setText(str(round(float(self.Ecav)+0.5, 3)))
        self.EC_1.setText(str(round(float(self.EC)-0.2, 3)))
        self.EC_2.setText(str(round(float(self.EC)+0.2, 3)))
        self.EJ_sum_1.setText(str(round(float(self.EJ_sum)-10, 3)))
        self.EJ_sum_2.setText(str(round(float(self.EJ_sum)+10, 3)))
        self.EJ_diff_1.setText(str(round(float(self.EJ_diff)-3, 3)))
        self.EJ_diff_2.setText(str(round(float(self.EJ_diff)+3, 3)))

        ### set min, max times 1000 (to int because Qslider works only with an integer)
        self.Ecav_slider.setMaximum(int((float(self.Ecav)+0.5)*1000))
        self.Ecav_slider.setMinimum(int((float(self.Ecav)-0.5)*1000))
        self.EC_slider.setMaximum(int((float(self.EC)+0.2)*1000))
        self.EC_slider.setMinimum(int((float(self.EC)-0.2)*1000))
        self.EJ_sum_slider.setMaximum(int((float(self.EJ_sum)+10)*1000))
        self.EJ_sum_slider.setMinimum(int((float(self.EJ_sum)-10)*1000))
        self.EJ_diff_slider.setMaximum(int((float(self.EJ_diff)+3)*1000))
        self.EJ_diff_slider.setMinimum(int((float(self.EJ_diff)-3)*1000))
        return
    
    
    def Ecav_slider_change(self):
        ### moving slider Ecav
        self.Ecav = self.Ecav_slider.value()/1000
        self.Ecav_value.setText(str(self.Ecav))
        self.plot_line_function()
        return

    def EC_slider_change(self):
        ### moving slider EC
        self.EC = self.EC_slider.value()/1000
        self.EC_value.setText(str(self.EC))
        self.plot_line_function()
        return
    
    def EJ_sum_slider_change(self):
        ### moving slider EJ_sum
        self.EJ_sum = self.EJ_sum_slider.value()/1000
        self.EJ_sum_value.setText(str(self.EJ_sum))
        self.plot_line_function()
        return
    
    def EJ_diff_slider_change(self):
        ### moving slider EJ_diff
        self.EJ_diff = self.EJ_diff_slider.value()/1000
        self.EJ_diff_value.setText(str(self.EJ_diff))
        self.plot_line_function()
        return


    def plot_graph_function(self):
        ### stitch data file
        stitch_data_files(self.MplWidget.canvas.axes, self.guiScans,
                  convert_to_flux=convert_to_flux,
                  offset=self.offset, 
                  volts_per_flux=self.volts_per_flux,
                  use_phase_data = use_phase_data,
                  freq_on_x_axis = self.freq_on_x_axis)
        self.MplWidget.canvas.axes.set_title('FM 16 spec_Overlay')
        self.MplWidget.canvas.axes.set_xlim(self.ylim_flux[0], self.ylim_flux[1])
        self.MplWidget.canvas.axes.set_ylim(self.xlim[0], self.xlim[1])
        self.MplWidget.canvas.draw_idle()
        return 


    def plot_line_function(self):
        EC_current = self.EC
        EJ_sum_current = self.EJ_sum
        EJ_diff_current = self.EJ_diff
        cav_freq_current = self.Ecav
        min_flux_current = np.min(self.ylim_flux)
        max_flux_current = np.max(self.ylim_flux)
        bare_current = self.bare
        two_photon_current = self.two_photon
        cav_assist_current = self.cav_assist
        freq_on_x_axis_current = self.freq_on_x_axis

        
        for line in self.lines:
            try:
                line.remove()
            except:
                continue
        self.lines = []
        
        overlay_GUI_transmon(self.guifig, self.MplWidget.canvas.axes, self.lines,
                          EC = EC_current, 
                          EJ_sum = EJ_sum_current,
                          EJ_diff = EJ_diff_current,
                          cav_freq = cav_freq_current,
                          min_flux= min_flux_current,
                          max_flux= max_flux_current,
                          flux_points=flux_points, 
                          bare = bare_current, 
                          twoPhoton = two_photon_current, 
                          cavAssist = cav_assist_current,
                          freq_on_x_axis = freq_on_x_axis_current)
        
        if self.freq_on_x_axis:
            self.MplWidget.canvas.axes.set_xlim(self.xlim[0], self.xlim[1])
            self.MplWidget.canvas.axes.set_ylim(self.ylim_flux[0], self.ylim_flux[1])
            pass
        else:
            self.MplWidget.canvas.axes.set_xlim(self.ylim_flux[0], self.ylim_flux[1])
            self.MplWidget.canvas.axes.set_ylim(self.xlim[0], self.xlim[1])
            pass
        return
    

    def vminmax_slider_change(self):
        self.guiScans = []
        for i in range(len(self.files_dir)):
            self.vmin = self.vmin_slider.value()
            self.vmax = self.vmax_slider.value()
            #xlims: specdata xaxis
            #ylims: specdata yaxis
            xlims = self.xlim
            ylims = self.ylim_voltage
            xSize = 6.5
            ySize = 4.0
            cmap = 'YlGnBu'

            data = replotter(self.files_dir[i],
                                vmin = self.vmin,
                                vmax = self.vmax,
                                cmap = cmap,
                                xlims = xlims,
                                ylims = ylims,
                                fig_xSize = xSize,
                                fig_ySize = ySize,
                                smoothe_data = False,
                                verbose = False,
                                GHzVHz_correction = True)
            self.guiScans.append(data)
        self.plot_graph_function()
        self.vmin_label.setText(str(self.vmin))
        self.vmax_label.setText(str(self.vmax))
        return


    def bare_transitions_plot(self):
        if self.bare_transition_chkbox.isChecked():
            self.bare = True
        else:
            self.bare = False
        self.plot_line_function()
        return

    def two_photon_plot(self):
        if self.two_photon_chkbox.isChecked():
            self.two_photon = True  
        else:
            self.two_photon = False
        self.plot_line_function()
        return

    def cav_assist_plot(self):
        if self.cav_assist_chkbox.isChecked():
            self.cav_assist = True
        else:
            self.cav_assist = False
        self.plot_line_function()
        return

    def freq_on_x_axis_function(self):
        if self.freq_on_x_axis_chkbox.isChecked():
            self.freq_on_x_axis = True
        else:
            self.freq_on_x_axis = False
        self.plot_graph_function()
        self.plot_line_function()
        return


    def prev_onclick(self):
        ### connects onclick function when plot is clicked
        self.cid = self.guifig.canvas.mpl_connect('button_press_event', self.onclick)
        return
    
    def onclick(self, event):
        ### reads 10 x, y coordinate and saves in the x_fit, y_fit array
        ix, iy = event.xdata, event.ydata
        print (f'x = {ix}, y = {iy}')
        self.coords.append((ix, iy))

        self.x_fit[self.rep_click] = ix
        self.y_fit[self.rep_click] = iy
        self.MplWidget.canvas.axes.scatter(ix, iy, s = 200, c = 'red')
        self.MplWidget.canvas.draw()
            
        if len(self.coords) == 10:
            self.guifig.canvas.mpl_disconnect(self.cid)
        self.rep_click += 1
        return


    def bare_transitions_function(self, phi, EJ_sum, EJ_diff, EC):
        d = EJ_diff/EJ_sum
        EJ_eff = EJ_sum *np.abs(np.cos(np.pi*phi)) * np.sqrt(np.abs(1 + d**2 * (np.tan(np.pi*phi))**2))
        return (np.sqrt(8*EC*EJ_eff))
    
    def bare_transitions_fit(self):
        ### fitting calculation
        popt, pcov = curve_fit(self.bare_transitions_function, self.x_fit, self.y_fit, method='lm', p0 = np.array([self.EJ_sum, self.EJ_diff, self.EC]), maxfev = 10000)
        print(popt)
        
        ### plot line_fit
        min_flux = np.min(self.ylim_flux)
        max_flux = np.max(self.ylim_flux)
        flux_points = 51
        flux = np.linspace(min_flux, max_flux, flux_points)

        for line in self.line_fit:
            try:
                line.remove()
            except:
                continue
        self.line_fit = []
        
        self.line_fit.append(self.MplWidget.canvas.axes.plot(flux, self.bare_transitions_function(flux, popt[0], popt[1], popt[2]), color='k', marker="None", linestyle='solid', linewidth = 3, label = 'Utility Function'))
        self.MplWidget.canvas.axes.legend(loc='lower right')
        
        ### reset the values in the slider
        self.EC_fit_value.setText(str(round(popt[2],3)))
        self.EJ_sum_fit_value.setText(str(round(popt[0],3)))
        self.EJ_diff_fit_value.setText(str(round(popt[1],3)))
        self.EC_value.setText(str(round(popt[2],3)))
        self.EJ_sum_value.setText(str(round(popt[0],3)))
        self.EJ_diff_value.setText(str(round(popt[1],3)))
        self.EC = popt[2]
        self.EJ_sum = popt[0]
        self.EJ_diff = popt[1]
        self.EC_slider.setValue((int(round(popt[2],3)*1000)))
        self.EJ_sum_slider.setValue(int(round(popt[0],3)*1000))
        self.EJ_diff_slider.setValue(int(round(popt[1],3)*1000))


        ### redraw the plot with the fitted curve and lines
        self.MplWidget.canvas.draw()
        self.plot_line_function()
        return
    

    def reset(self):
        ### reset params
        self.graph_on = False
        self.freq_on_x_axis = False
        self.bare = False
        self.two_photon = False
        self.cav_assist = False
        self.Ecav = 7.484
        self.EC = 0.3
        self.EJ_sum = 43.5
        self.EJ_diff = 5.5
        self.Ecav_1.setText(str(round(float(self.Ecav)-0.5, 3)))
        self.Ecav_2.setText(str(round(float(self.Ecav)+0.5, 3)))
        self.EC_1.setText(str(round(float(self.EC)-0.2, 3)))
        self.EC_2.setText(str(round(float(self.EC)+0.2, 3)))
        self.EJ_sum_1.setText(str(round(float(self.EJ_sum)-10, 3)))
        self.EJ_sum_2.setText(str(round(float(self.EJ_sum)+10, 3)))
        self.EJ_diff_1.setText(str(round(float(self.EJ_diff)-5, 3)))
        self.EJ_diff_2.setText(str(round(float(self.EJ_diff)+5, 3)))
        self.Ecav_slider.setMaximum(int((float(self.Ecav)+0.5)*1000))
        self.Ecav_slider.setMinimum(int((float(self.Ecav)-0.5)*1000))
        self.EC_slider.setMaximum(int((float(self.EC)+0.2)*1000))
        self.EC_slider.setMinimum(int((float(self.EC)-0.2)*1000))
        self.EJ_sum_slider.setMaximum(int((float(self.EJ_sum)+10)*1000))
        self.EJ_sum_slider.setMinimum(int((float(self.EJ_sum)-10)*1000))
        self.EJ_diff_slider.setMaximum(int((float(self.EJ_diff)+5)*1000))
        self.EJ_diff_slider.setMinimum(int((float(self.EJ_diff)-5)*1000))
        self.Ecav_slider.setValue(int((float(self.Ecav)-0.5)*1000))
        self.EC_slider.setValue(int((float(self.EC)-0.5)*1000))
        self.EJ_diff_slider.setValue(int((float(self.EJ_diff)-0.5)*1000))
        self.EJ_sum_slider.setValue(int((float(self.EJ_sum)-0.5)*1000))

        ### reset objects
        self.guiScans = []
        self.files_str = ''
        self.load_file_txtbox.setText(self.files_str)
        self.files_dir = []
        self.lines = []
        self.coords = []
        self.x_fit = np.zeros(10)
        self.y_fit = np.zeros(10)
        self.rep_click = 0
        self.line_fit = []
        print(self.lines, self.guiScans)
        self.flux_range = []    
        self.xlim = []         
        self.ylim_flux = []            
        self.ylim_voltage = []    
        self.vmin = -7                                 
        self.vmax = 2
        self.vmin_label.setText(str(self.vmin))
        self.vmax_label.setText(str(self.vmax))
        self.vmin_slider.setMaximum(self.vmin + 4)
        self.vmin_slider.setMinimum(self.vmin - 4)
        self.vmax_slider.setMaximum(self.vmax + 4)
        self.vmax_slider.setMinimum(self.vmax - 4)
        self.vmin_slider.setValue(self.vmin)
        self.vmax_slider.setValue(self.vmax)


        ### redraw plots
        self.MplWidget.canvas.axes.cla()
        self.MplWidget.canvas.draw()
        return




##############################
# transmon params
##############################
'''
EC = 0.208
EJ_sum = 42.40
EJ_diff = 12.93
Ecav = 7.484
'''


min_freq = 2.0
max_freq = 9.3

min_flux = -0.75
max_flux = 0.5
flux_points = 51



###########
#decide on view orientation
###########
convert_to_flux = True
use_phase_data = False





app = QApplication([])
window = start_gui()
window.show()
app.exec_()