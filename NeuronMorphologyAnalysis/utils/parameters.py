import os
def set_parameters(parameters_dict):
    project = parameters_dict['project']
    base_dir = r"C:\Users\Judith\OneDrive - Allen Institute\Desktop"
    base_dir = r"C:\Users\judith.baka\OneDrive - Allen Institute\Desktop"
    if project =='medulla':
        parameters_dict['metadata_dir'] = os.path.join(base_dir,r'MouseLight\MO cells MouseLight\metadata')
        parameters_dict['save_dir'] = os.path.join(base_dir,r'MouseLight\MO cells MouseLight')
        parameters_dict['online_notebook_ID'] = 'MO cells'
        #juci_cluster_colors = ['r','b','g','y']
    elif project =='cortex':
        parameters_dict['metadata_dir'] = os.path.join(base_dir,r'MouseLight\Cortex cells MouseLight\metadata')
        parameters_dict['save_dir'] = os.path.join(base_dir,r'MouseLight\Cortex cells MouseLight')
        parameters_dict['online_notebook_ID'] = 'CORTEX cells'
    elif project =='thalamus':
        parameters_dict['metadata_dir'] = os.path.join(base_dir,r'MouseLight\Thalamus cells MouseLight\metadata')
        parameters_dict['save_dir'] = os.path.join(base_dir,r'MouseLight\Thalamus cells MouseLight')
        parameters_dict['online_notebook_ID'] = 'Thalamus cells'
        #juci_cluster_colors = ['r','b','g','y']
    else:
        print('unknown project: {}'.format(project))
        parameters_dict['metadata_dir'] = None
        parameters_dict['save_dir'] = None
        parameters_dict['online_notebook_ID'] = None
    
    if float(parameters_dict['ccf_version']) == 3.:
        parameters_dict['json_dir'] = os.path.join(base_dir,r'MouseLight\all_cells_json\unsorted-ccf30')
        parameters_dict['path_allen_df_with_volumes'] = os.path.join(base_dir,r'MouseLight\Allen_p56_mouse_annotation\allen_df_with_volumes_ccf3.csv')
    elif float(parameters_dict['ccf_version']) == 2.5:
        parameters_dict['json_dir'] = os.path.join(base_dir,r'MouseLight\all_cells_json\unsorted')
        parameters_dict['path_allen_df_with_volumes'] = os.path.join(base_dir,r'MouseLight\Allen_p56_mouse_annotation\v2.5\allen_df_with_volumes.csv')
    else:
        print('unknown ccf version :{}'.format(parameters_dict['ccf_version']))
        parameters_dict['json_dir']=None
        parameters_dict['path_allen_df_with_volumes'] = None
        
    
    return parameters_dict