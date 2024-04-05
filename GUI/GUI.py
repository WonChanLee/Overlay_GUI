# PyQt5 packages
from PyQt5.QtWidgets import QMainWindow, QFileDialog
from PyQt5.uic import loadUi
from PyQt5 import QtGui

# matplotlib packages
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)

# numpy, os packages
import numpy as np
import pickle
import os

# custom libraries
from replotter.replotter_module import replotter
from fitting.scqubits_fit import fit_function, qubit_energies

from overlay_helpers import stitch_data_files
from overlay_helpers import overlay_GUI_transmon

# window taskbar icon
import ctypes
myappid = 'mycompany.myproduct.subproduct.version' # arbitrary string
ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

convert_to_flux = True
use_phase_data = False

min_freq = 2.0
max_freq = 9.3

min_flux = -0.75
max_flux = 0.5
flux_points = 51

##############################
# Define Gui
##############################
class GUI_template(QMainWindow):
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
    Upto date, only bare transition can be fitted.
    It will be saved as coordinates that is used as the fitting data
    7. click "Bare Transition Fit" 
    8. Fitted result will be displayed in the box.
    9. To reset everything, click "Reset".
    This will initialize every variable and plot.
    '''
    
    def __init__(self, params, starting_vals, lines, line_vals, qubit_type):
        QMainWindow.__init__(self)
        file_loc = os.path.dirname(os.path.realpath(__file__))
        gui_path = os.path.join(file_loc,'Overlay_Gui.ui')
        loadUi(gui_path,self)
        self.setWindowTitle("Overlay GUI")
        self.setWindowIcon(QtGui.QIcon('hyperbolic.png'))
        self.addToolBar(NavigationToolbar(self.MplWidget.canvas, self))
        

        ### Default params
        self.graph_on = False
        self.freq_on_x_axis = False
        self.params = params
        self.starting_vals = starting_vals
        self.overlay_lines = lines
        self.overlay_bools = line_vals
        self.qubit_type = qubit_type
        self.set_starting_params()

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
        self.xlim = [min_freq, max_freq]                                 # define the xlim of loaded files; default: freq
        self.ylim_flux = [min_flux, max_flux]                            # define the ylim_flux of loaded files
        self.ylim_voltage = []                         # define the ylim_voltage of loaded files
        self.vmin = -7                                 # define the vmin
        self.vmax = 2                                  # define the vmax
        self.vmin_slider.setMaximum(self.vmin + 4)
        self.vmin_slider.setMinimum(self.vmin - 4)
        self.vmax_slider.setMaximum(self.vmax + 4)
        self.vmax_slider.setMinimum(self.vmax - 4)

        self.volts_per_flux = 1
        self.offset = 0
        self.volts_per_flux_slider.setMaximum(self.volts_per_flux*2*1000)
        self.volts_per_flux_slider.setMinimum(self.volts_per_flux/2*1000)
        self.offset_slider.setMaximum(self.volts_per_flux*1000)
        self.offset_slider.setMinimum(-self.volts_per_flux*1000)
        self.volts_per_flux_label.setText(str(self.volts_per_flux))
        self.offset_label.setText(str(self.offset))
        ### load file 
        self.load_file_button.clicked.connect(self.loadfiles)

        ### set transmon params, set slider range
        #self.set_params.clicked.connect(self.set_transmon_params)

        ### plot graph
        self.plot_graph.clicked.connect(self.plot_graph_function)

        ### set vmin, vamx and replot graph
        self.vmin_slider.sliderMoved.connect(self.vminmax_slider_change)
        self.vmax_slider.sliderMoved.connect(self.vminmax_slider_change)

        self.volts_per_flux_slider.sliderMoved.connect(self.volts_per_flux_slider_change)
        self.offset_slider.sliderMoved.connect(self.volts_per_flux_slider_change)

        ### check toggles
        self.freq_on_x_axis_chkbox.stateChanged.connect(self.freq_on_x_axis_function)
        # Connect checkbox buttons for the overlay script
        self.connect_checkboxes()

        # Connect all the relevant sliders (from list of params)
        self.format_sliders()
        self.connect_sliders()

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
        print(self.files_dir_unsorted)
        for file in self.files_dir_unsorted[0]:
            if file in self.files_dir:
                print('file already in the list')
                continue
            else:
                self.files_dir.append(file)
        self.files_str = ''
        ### make the set of metadata files that will be included in the GUI
        for file in self.files_dir:
            with open(file,'rb') as f:
                pickledict = pickle.load(f)
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

            data = replotter(file,
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
            self.files_str += str(file + '\n')
        if self.files_str:
            self.load_file_txtbox.setText(self.files_str)
            self.xlim = [np.min(xlim_search), np.max(xlim_search)]
            self.ylim_flux = [np.min(ylim_search_flux), np.max(ylim_search_flux)]
            self.ylim_voltage = [np.min(ylim_search_voltage), np.max(ylim_search_voltage)]
        # print(self.xlim)
        # print(self.ylim_flux)
        # print(self.ylim_voltage)
        return
            
    def set_params_auto(self):
        # Updating the behavior of this function: we check if the user
        # has provided reference/ center values and if not, we read off 
        # the current values from the sliders
        ### read written values 
        for param in self.params:
            slider_name = param
            try:
                value = self.get_ref_value(slider_name)
            except:
                print('No ref val found, using current slider value')
                value = self.get_slider_value(slider_name)
            
            self.__setattr__(param, value)
        self.plot_line_function()
        return

    def get_ref_value(self, slider_name):
        # Read the value store in the textbox under the slider
        # this sets the middle point of the sliders after calling
        # "format_sliders"
        reference_label = slider_name+'_ref_txt'
        return_val = self.__getattribute__(reference_label)
        return float(return_val.toPlainText())
    
    def get_slider_value(self, slider_name):
        # The sliders are defined in integer values so we scale the value
        # by 1000 when setting it/ reading it back
        slider_label = slider_name+'_slider'
        return_val = self.__getattribute__(slider_label)
        return return_val.value()/1000
    
    def get_params_auto(self):
        # Create a dictionary of the fit params for the relevant qubit type
        # The fit and overlay functions should take this in as an argument to
        # make things more consistent
        param_dict = {}
        for param in self.params:
            param_dict[param] = self.get_slider_value(param)
        for param in self.overlay_lines:
            param_dict[param] = self.get_checkbox_value(param)
        return param_dict
    
    def set_starting_params(self):
        # Initialize starting values for the overlay for the relevant qubit params
        for param, val in zip(self.params, self.starting_vals):
            self.__setattr__(param, val)
        # Initialize the starting booleans for the line overlay code
        for param, val in zip(self.overlay_lines, self.overlay_bools):
            self.__setattr__(param, val)

    ### Configuring the slider backend connection (GUI to code)
    def connect_sliders(self):
        # This function loops through all the sweepable parameters for the qubit
        # and connects the GUI slider object to the update function that actually 
        # sweeps the parameter of interest
        for param in self.params:
            slider_name = param+'_slider'
            slider = self.__getattribute__(slider_name)
            slider.sliderMoved.connect(self.change_func(param))

    def change_func(self, param):
        # This function returns the "update" function needed to connect the GUI
        # slider to the code. It reads the slider value, updates the store value
        # then calls the plot_line_function to update the overlay 
        def slider_func():
            slider_name = param+'_slider'
            slider = self.__getattribute__(slider_name)
            value = slider.value()/1000
            self.__setattr__(param, value)
            # Update the display value on the given slider
            self.__getattribute__(param+'_value').setText(str(value))
            self.plot_line_function()
        return slider_func

    ### Configuring the checkboxes
    def connect_checkboxes(self):
        # Same idea as with the sliders, here we generate and connect
        # the checkboxes relating to the relevant lines proceduraly 
        for line in self.overlay_lines:
            checkbox_name = line+'_chkbox'
            checkbox = self.__getattribute__(checkbox_name)
            checkbox.stateChanged.connect(self.checkbox_func(line))

    def checkbox_func(self, name):
        # Same magic function as above: we return a function that PyQT
        # then uses to conenct the GUI button to the actual checkbox boolean
        def check_box():
            checkbox_name = name+'_chkbox'
            box = self.__getattribute__(checkbox_name)
            self.__setattr__(name, box.isChecked())
            self.plot_line_function()
        return check_box
    
    def get_checkbox_value(self, name):
        checkbox = self.__getattribute__(name+'_chkbox')
        return checkbox.isChecked()
    
    def plot_graph_function(self):
        ### stitch data file
        self.MplWidget.canvas.axes.cla()
        stitch_data_files(self.MplWidget.canvas.axes, self.guiScans,
                  convert_to_flux=convert_to_flux,
                  offset=self.offset, 
                  volts_per_flux=self.volts_per_flux,
                  use_phase_data = use_phase_data,
                  freq_on_x_axis = self.freq_on_x_axis)
        self.MplWidget.canvas.axes.set_title('FM 16 spec_Overlay')
        #self.MplWidget.canvas.axes.set_xlim(self.ylim_flux[0], self.ylim_flux[1])
        self.MplWidget.canvas.axes.set_ylim(self.xlim[0], self.xlim[1])
        self.MplWidget.canvas.draw_idle()
        return 

    def plot_line_function(self):
        params_dict = self.get_params_auto()
        EC_current = params_dict['EC']
        EJ_sum_current = params_dict['EJ_sum']
        EJ_diff_current = params_dict['EJ_diff']
        cav_freq_current = params_dict['Ecav']
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
    
    def volts_per_flux_slider_change(self):
        self.volts_per_flux = self.volts_per_flux_slider.value()/1000
        self.offset = self.offset_slider.value()/1000
        self.plot_graph_function()
        self.plot_line_function()
        self.volts_per_flux_label.setText(str(self.volts_per_flux))
        self.offset_label.setText(str(self.offset))


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
        self.MplWidget.canvas.axes.scatter(ix, iy, s = 150, c = 'red')
        self.MplWidget.canvas.draw()
            
        if len(self.coords) == 10:
            self.guifig.canvas.mpl_disconnect(self.cid)
        self.rep_click += 1
        return

    def bare_transitions_fit(self):
        fit_params = fit_function(self.get_params_auto(), self.qubit_type, self.x_fit, self.y_fit)
        print(fit_params)
        ### plot line_fit
        min_flux = np.min(self.ylim_flux)
        max_flux = np.max(self.ylim_flux)
        flux_points = 51
        flux = np.linspace(min_flux, max_flux, flux_points)

        qubit_lines = qubit_energies(fit_params, self.qubit_type, flux)
        
        for line in self.line_fit:
            try:
                line.remove()
            except:
                continue
        self.line_fit = []
        
        self.line_fit.append(self.MplWidget.canvas.axes.plot(flux, qubit_lines, color='k', 
                                                             marker="None", linestyle='solid', linewidth = 3, 
                                                             label = 'Fit result'))
        
        self.MplWidget.canvas.axes.legend(loc='lower right')
        
        ### reset the values in the slider
        self.set_fit_params(fit_params)

        ### redraw the plot with the fitted curve and lines
        self.MplWidget.canvas.draw()
        self.plot_line_function()
        return

    def set_fit_params(self, fit_params):
        for param in fit_params.keys():
            value = fit_params[param]
            str_val = str(round(value, 3))
            self.__setattr__(param, value)
            self.__getattribute__(param+'_fit_value').setText(str_val)
            self.__getattribute__(param+'_value').setText(str_val)
            self.__getattribute__(param+'_slider').setValue(int(value*1000))

    def configure_slider(self, name, val, val_range):
        '''
        This function does the dirty work behind configuring a slider:
        it computes the min/max vals and then updates the text as well 
        as the slider range. We probably should also just update the 
        value of the slider so we end up "centered" in our range 
        '''
        min_val = val-val_range
        max_val = val+val_range
        str_value = str(round(val,3))
        slider_name = name+'_slider'
        min_label = name+'_1'
        max_label = name+'_2'
        value_label = name+'_value'
        min_string = str(round(min_val,3))
        max_string = str(round(max_val,3))
        #Configures labels
        self.__getattribute__(min_label).setText(min_string) 
        self.__getattribute__(max_label).setText(max_string) 
        self.__getattribute__(value_label).setText(str_value)
        #Configures actual slider vals
        slider = self.__getattribute__(slider_name)
        slider.setMinimum(int(min_val*1000))
        slider.setMaximum(int(max_val*1000))
        slider.setValue(int(val*1000))

    def format_sliders(self):
        # Base function to create sliders. This gets overwritten in the child class
        # with whatever range the user would like to use. Here I defined a unit 
        # range so that the slider object is correctly initialized so we can query
        # it without receiving a "0" value
        for param in self.params:
            self.configure_slider(param, self.__getattribute__(param),1)

    def reset(self):
        ### reset params
        self.graph_on = False
        self.freq_on_x_axis = False
        self.set_starting_params()

        self.format_sliders()
        
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
