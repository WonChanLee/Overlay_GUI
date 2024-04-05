import pickle
import os

folder = os.path.dirname(__file__)
file_name = r'spec_scan_readout_avoided_test_5_20230726_233342.pkl'
filepath = os.path.join(folder, file_name)
pickledict = pickle.load(open(filepath, "rb" ) )

print(pickledict)


'''
a = objectt['Data']
b = a['specdata']
print(b)
'''