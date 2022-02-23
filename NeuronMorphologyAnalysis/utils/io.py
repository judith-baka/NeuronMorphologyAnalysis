import tkinter as tk
from tkinter import filedialog
import numpy as np
import os
def load_raw_data_matrices(parameters_dict):    
    root = tk.Tk()
    root.withdraw()
    file_to_load = filedialog.askopenfilename(initialdir = parameters_dict['save_dir'])
    print('loading {}'.format(file_to_load))
    original_data = np.load(os.path.join(parameters_dict['save_dir'],file_to_load),allow_pickle = True)
    original_data = original_data.tolist()
    parameters_dict = original_data['parameters_dict']
    parameters_dict['loaded_data_file'] = file_to_load
    allen_df = original_data['allen_df']
    return original_data, parameters_dict,allen_df