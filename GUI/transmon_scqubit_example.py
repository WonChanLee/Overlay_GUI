# PyQt5 packages
from PyQt5.QtWidgets import QApplication

from GUI import GUI_template

# window taskbar icon
import ctypes
myappid = 'mycompany.myproduct.subproduct.version' # arbitrary string
ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

##############################
# Define Gui
##############################
class start_gui(GUI_template):
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
    
    def __init__(self):
        params = ['Ecav','EC','EJ_sum','EJ_diff']
        starting_vals = [7.5,0.2,45,10]
        lines = ['bare','two_photon','cav_assist']
        line_bools = [False,False,False]
        qubit_type = "Transmon"

        super().__init__(params, starting_vals, lines, line_bools, qubit_type)
        
        ### Default params
        self.graph_on = False
        self.freq_on_x_axis = False
        
        self.volts_at_half = -0.182                                      # Not user defined parameter yet. Should enter in the script 
        self.volts_per_flux = 1.169718                                     # Not user defined parameter yet. Should enter in the script
        self.offset = self.volts_at_half + self.volts_per_flux/2        # Not user defined parameter yet. Should enter in the script
        self.format_sliders()

        ### set transmon params, set slider range
        self.set_params.clicked.connect(self.set_transmon_params)

        ### Fitting
        self.fit_button.clicked.connect(self.bare_transitions_fit)

        ### Reset
        self.reset_button.clicked.connect(self.reset)
              
    def set_transmon_params(self):
        # Updating the behavior of this function: we check if the user
        # has provided reference/ center values and if not, we read off 
        # the current values from the sliders
        ### read written values 
        self.set_params_auto()
        self.format_sliders()
        return

    def format_sliders(self):
        self.configure_slider('Ecav', self.Ecav, 0.25)
        self.configure_slider('EC', self.EC, 0.1)
        self.configure_slider('EJ_sum', self.EJ_sum, 10)
        self.configure_slider('EJ_diff', self.EJ_diff, 5)


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