import os
from distutils.sysconfig import get_python_lib

site_directory = get_python_lib()

package_directory = os.path.dirname(os.path.realpath(__file__))

filepath = os.path.join(site_directory,"Overlay_GUI.pth")

print(site_directory)
print(package_directory)
print(filepath)

f = open(filepath,"w")
f.write(package_directory)