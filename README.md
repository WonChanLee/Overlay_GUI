# Overlay-GUI
Stand alone replotting tool and fitting overlay for spec data
 
## General description
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

### Steps
1. run the gui
2. load pkl files and plot graph
3. adjust parameters 
- vmin, vmax - to make more clear plot
- Ecav, EC, EJ, EJ_diff - make reasonable initial guess for fitting
4. If Ecav, EC, EJ, EJ_diff need to be changed, write values and click "Set Params".
5. choose the lines to see
6. click "Bare Transition Ref". 
10 reference points on the plot.
- Upto 11.08.23, only bare transition can be fitted.
- It will be saved as coordinates that is used as the fitting data
7. click "Bare Transition Fit" 
8. Fitted result will be displayed in the box.
9. To reset everything, click "Reset".
This will initialize every variable and plot.
 

## To do list:
- Overlay/ fitting tool:
	- [ ] Generate models for qubits of interest (fluxonium, transmon, general avoided crossing)
	- [ ] Implement SC Qubit package