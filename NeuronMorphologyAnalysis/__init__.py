import imagej

#local_fiji_dir =  r"C:\Users\judith.baka\OneDrive - Allen Institute\Desktop\MouseLight\Fiji.app" #getpath().decode('utf-8')
local_fiji_dir =  r'C:\Fiji.app' #getpath().decode('utf-8')
#local_fiji_dir =  r"C:\Users\Judith\OneDrive - Allen Institute\Desktop\MouseLight\Fiji.app" #getpath().decode('utf-8')
try:
    ij = imagej.init(local_fiji_dir, headless=False)
except:
    print('fiji already initiated?')
from .analysis import analysis
from .utils import io, parameters, ccf, online_notebook, snt, plot

