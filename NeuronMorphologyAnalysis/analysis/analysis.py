import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import scipy 
import scipy.ndimage as ndimage
import datetime
from ..utils import ccf, snt
import os
import miniball

#%%

def generate_sensory_motor_groups(original_data):
    """
    creates clusters based on being motor or sensory related

    Parameters
    ----------
    original_data : dict
        data extracted from json files with the analyze_json_files() function

    Returns
    -------
    output_dict : dict

    """
    cluster_names = {'Sensory_related':'Sensory related', 
                     'Motor_related':'Motor related'}
    cluster_colors = {'Sensory_related':'green', 
                     'Motor_related':'red'}
    cluster_indices = {'Sensory_related':np.asarray(original_data['soma_area_in_the_medulla']) == 'Sensory related', 
                       'Motor_related':np.asarray(original_data['soma_area_in_the_medulla']) == 'Motor related'}
    
        
    output_dict = {'cluster_indices':cluster_indices,
                   'cluster_colors':cluster_colors,
                   'cluster_names':cluster_names,
                   'clustering_description': 'Cells were assigned to sensory or motor related nuclei.'}
    
    return output_dict 


def generate_sensory_motor_groups_with_projections(original_data): # NINCSENKESZEN!!# TODO
    """
    creates clusters based on being motor or sensory related, but separates
    cerebellum projecting cells from motor related nuclei and
    midbrain projecting cells from sensory nuclei

    Parameters
    ----------
    original_data : dict
        data extracted from json files with the analyze_json_files() function

    Returns
    -------
    output_dict : dict

    """
    cluster_names = {'Sensory_related_midbrain_projecting':'Sensory related midbrain projecting',
                     'Sensory_related_cerebellum_projecting':'Sensory related cerebellum projecting',
                     'Sensory_related_rest':'Rest of sensory related', 
                     'Motor_related_cerebellum':'Motor related cerebellum projecting',
                     'Motor_related_rest':'Motor related rest'}
    cluster_colors = {'Sensory_related_midbrain_projecting':'yellowgreen',
                      'Sensory_related_cerebellum_projecting':'orange',
                     'Sensory_related_rest':'darkgreen', 
                     'Motor_related_cerebellum':'red',
                     'Motor_related_rest':'brown'}
    
    cluster_indices = {'Sensory_related_midbrain_projecting':np.asarray(np.ones(len(original_data['cell_names']))*False,bool),
                       'Sensory_related_cerebellum_projecting':np.asarray(np.ones(len(original_data['cell_names']))*False,bool),
                     'Sensory_related_rest':np.asarray(np.ones(len(original_data['cell_names']))*False,bool), 
                     'Motor_related_cerebellum':np.asarray(np.ones(len(original_data['cell_names']))*False,bool),
                     'Motor_related_rest':np.asarray(np.ones(len(original_data['cell_names']))*False,bool)}
    
    for i,(motor_vs_sensory, projection_class) in enumerate(zip(original_data['soma_area_in_the_medulla'],original_data['cell_juci_clusters'])):
        if motor_vs_sensory == 'Motor related':
            if projection_class == 'Cerebellum':
                cluster_indices['Motor_related_cerebellum'][i] = True
            else:
                cluster_indices['Motor_related_rest'][i] = True
        elif motor_vs_sensory == 'Sensory related':
            if projection_class == 'Midbrain':
                cluster_indices['Sensory_related_midbrain_projecting'][i] = True
            elif projection_class == 'Cerebellum':
                cluster_indices['Sensory_related_cerebellum_projecting'][i] = True
            else:
                cluster_indices['Sensory_related_rest'][i] = True
        else:
            print('nor motor nor sensory related cell.. wtf?')
    
        
    output_dict = {'cluster_indices':cluster_indices,
                   'cluster_colors':cluster_colors,
                   'cluster_names':cluster_names,
                   'clustering_description': 'Cells were assigned to sensory or motor related nuclei. Some projection groups were removed from these big clusters.'}
    
    return output_dict 

def generate_soma_location_groups(original_data):
    """
    creates clusters based on ground truth soma location - not raw ccf 
    .. the so called "Lauren cluters"

    Parameters
    ----------
    original_data : dict
        data extracted from json files with the analyze_json_files() function

    Returns
    -------
    output_dict : dict

    """
    soma_locations_header= np.unique(original_data['soma_locations'])
    
    cluster_names = {'trigeminal_cplx':'Trigeminal complex', 
                     'reticular_formation':'Reticular formation',
                     'vestibular_cplx':'Vestibular complex',
                     'dorsal_column':'Dorsal column nucleus',
                     'lateral_reticular_nucl': 'Lateral reticular nucleus',
                     'motor_nucl': 'Motor nuclei'}
    cluster_colors = {'trigeminal_cplx':'darkorange', 
                     'reticular_formation':'c',
                     'vestibular_cplx':'darkviolet',
                     'dorsal_column':'orchid',
                     'lateral_reticular_nucl':'black',
                     'motor_nucl':'yellow'}
    
    lauren_cluster_nuclei = {'trigeminal_cplx':['Spinal nucleus of the trigeminal, caudal part',
                                      'Spinal nucleus of the trigeminal, interpolar part', 
                                      'Spinal nucleus of the trigeminal, oral part'],
                             'reticular_formation':['Intermediate reticular nucleus',
                                           'Gigantocellular reticular nucleus',
                                           'Paragigantocellular reticular nucleus',
                                           'Paragigantocellular reticular nucleus, lateral part',
                                           'Paragigantocellular reticular nucleus, dorsal part',
                                           'Parvicellular reticular nucleus',
                                           'Magnocellular reticular nucleus',
                                           'Medullary reticular nucleus',
                                           'Medullary reticular nucleus, dorsal part'],
                             'vestibular_cplx':['Lateral vestibular nucleus',
                                           'Medial vestibular nucleus',
                                           'Spinal vestibular nucleus',
                                           'Superior vestibular nucleus'],
                             'dorsal_column':['Cuneate nucleus',
                                     'Gracile nucleus',
                                     'External cuneate nucleus'],
                             'lateral_reticular_nucl':['Lateral reticular nucleus'],
                             'motor_nucl':['Facial motor nucleus',
                                     'Hypoglossal nucleus']}
    
    lauren_cluster_indxes = {'trigeminal_cplx':np.ones(len(original_data['cell_names']))*False,
                             'reticular_formation':np.ones(len(original_data['cell_names']))*False,
                             'vestibular_cplx':np.ones(len(original_data['cell_names']))*False,
                             'dorsal_column':np.ones(len(original_data['cell_names']))*False,
                             'lateral_reticular_nucl':np.ones(len(original_data['cell_names']))*False,
                             'motor_nucl':np.ones(len(original_data['cell_names']))*False}
    
    

    
    for cluster_key_now in lauren_cluster_nuclei.keys():
        for soma_now in lauren_cluster_nuclei[cluster_key_now]:
            nucl_idx = (soma_locations_header == soma_now)
            if sum(nucl_idx)==0:
                continue
            lauren_cluster_indxes[cluster_key_now] = (lauren_cluster_indxes[cluster_key_now] + np.asarray(original_data['soma_locations_list'])[:,nucl_idx].flatten())>0
    
    output_dict = {'cluster_indices':lauren_cluster_indxes,
                   'cluster_colors':cluster_colors,
                   'cluster_names':cluster_names,
                   'cluster_nuclei':lauren_cluster_nuclei,
                   'clustering_description': 'Somas were assigned to nuclei by Judith and Lauren.'}
    
    return output_dict    

def generate_main_projection_groups(original_data):
    cluster_names = {'Cerebellum':'Cerebellum', 
                     'Medulla_Pons':'Medulla/Pons',
                     'Midbrain':'Midbrain', 
                     'Thalamus_Hypothalamus':'Thalamus/Hypothalamus',
                     'Spinal_cord': 'Spinal cord'}
    cluster_colors = {'Cerebellum':'red', 
                     'Medulla_Pons':'dodgerblue',
                     'Midbrain':'limegreen', 
                     'Thalamus_Hypothalamus':'gray',
                     'Spinal_cord':'yellow'}
    cluster_indices = {}
    
    for juci_cluster in np.unique(original_data['cell_juci_clusters']):
        cluster_indices[juci_cluster] = np.asarray(original_data['cell_juci_clusters'])==juci_cluster
        
        
    output_dict = {'cluster_indices':cluster_indices,
                   'cluster_colors':cluster_colors,
                   'cluster_names':cluster_names,
                   'clustering_description': 'Cells were assigned to main projection groups by Judith.'}
    
    return output_dict  

def generate_projection_groups(original_data):
    """
    creates clusters based on ccf projection patterns - a cell can belong to multiple clusters

    Parameters
    ----------
    original_data : dict
        data extracted from json files with the analyze_json_files() function

    Returns
    -------
    output_dict : dict

    """
    allen_df = original_data['allen_df']
    minimum_axon_end_point_num_forprojections = 3
    cluster_names = {'prereticular':'Pre-reticular cells', 
                     'premotor':'Pre-motor cells'}
    cluster_colors = {'prereticular':'olive', 
                     'premotor':'gold'}
    cluster_indxes = {'prereticular':np.ones(len(original_data['cell_names']))*False,
                      'premotor':np.ones(len(original_data['cell_names']))*False}
    
    cluster_target_nuclei = {'prereticular':['Intermediate reticular nucleus',
                                             'Gigantocellular reticular nucleus',
                                             'Paragigantocellular reticular nucleus, lateral part',
                                             'Paragigantocellular reticular nucleus, dorsal part',
                                             'Parvicellular reticular nucleus',
                                             'Magnocellular reticular nucleus',
                                             'Medullary reticular nucleus'],
                             'premotor':['Facial motor nucleus',
                                         'Hypoglossal nucleus',
                                         'Nucleus ambiguus',
                                         'Motor nucleus of trigeminal']}
    
    for cluster_key in cluster_target_nuclei.keys():
        lista = list()
        for target_now in cluster_target_nuclei[cluster_key]:
            idx = np.where(allen_df['name']==target_now)[0][0]
            enpoints_ipsi = original_data['allen_axon_end_points_matrix'][:,idx]*original_data['axon_end_point_numbers']
            enpoints_contra = original_data['allen_axon_end_points_matrix'][:,idx+len(allen_df)]*original_data['axon_end_point_numbers']
            endpoints = enpoints_ipsi + enpoints_contra
            lista.append(endpoints>=minimum_axon_end_point_num_forprojections)
            #break
        lista = np.asarray(lista)
        cluster_indxes[cluster_key] = np.sum(lista,0)>0
    
    output_dict = {'cluster_indices':cluster_indxes,
                   'cluster_colors':cluster_colors,
                   'cluster_names':cluster_names,
                   'cluster_target_nuclei':cluster_target_nuclei,
                   'clustering_description': 'A cell belongs to a cluster if it has at least {} endpoints in the corresponding nuclei'.format(minimum_axon_end_point_num_forprojections)}
    
    return output_dict    
        

def generate_cell_list(allen_df,original_data,basic_group,soma_locations_needed=None,axon_projections_needed=None):
   
    if type(basic_group) == list:
        needed_cells =basic_group
    elif basic_group == 'light red':
        needed_cells = ['AA1082',
                        'AA1338',
                        'AA1330',
                        'AA1084', 
                        'AA1345', 
                        'AA1358', 
                        'AA1064', 
                        'AA1353', 
                        'AA1068', 
                        'AA1347', 
                        'AA1083', 
                        'AA1357', 
                        'AA1336', 
                        'AA1342', 
                        'AA1065', 
                        'AA1360', 
                        'AA0950', 
                        'AA1340', 
                        'AA1093', 
                        'AA1351', 
                        'AA1085', 
                        'AA1344', 
                        'AA1346', 
                        'AA1070', 
                        'AA1341', 
                        'AA1325', 
                        'AA1339',
                        'AA1335'] #vilagos pirosak  - 'all' #(1355 talan sotetkek)
    elif basic_group == 'dark red':
        needed_cells = ['AA1077',
                        'AA1062',
                        'AA0951',
                        'AA0431', 
                        'AA0922', 
                        'AA0953', 
                        'AA1063',
                        'AA1404'] #sotet piros csopi #(1404 talan sotetkek)
        
    elif basic_group == 'dark green':
        needed_cells = ['AA0989',
                        'AA1426',
                        'AA1431',
                        'AA1535'] #sotet zold csopi (talan mind sotet kek)
    elif basic_group == 'light green':
        needed_cells = ['AA1327',
                        'AA1333',
                        'AA1337',
                        'AA1402', 
                        'AA1328', 
                        'AA1343', 
                        'AA1334',
                        'AA1363',
                        'AA1516'] #vilagos zold csopi
        
    elif basic_group == 'light blue':
        needed_cells = ['AA0948',
                        'AA1356',
                        'AA1355', # (1355 talan sotetkek)
                        'AA1354',
                        'AA1352',
                        'AA0952',
                        'AA1075',
                        'AA1454',
                        'AA1422',
                        'AA1424',
                        'AA1427',
                        'AA1520'] #vilagos kek kaosz csopi 
    elif basic_group == 'dark blue':
        needed_cells = ['AA1196',
                        'AA1313',
                        'AA0947',
                        'AA1331',
                        'AA1405',
                        'AA1329',
                        'AA1332',
                        'AA1348',
                        'AA1349',
                        'AA1359',
                        'AA1361',
                        'AA1362',
                        'AA0996',
                        'AA1350',
                        'AA1403',
                        'AA0504',
                        'AA0508',
                        'AA0516',
                        'AA0517',
                        'AA0434',
                        'AA0430',
                        'AA1489',
                        'AA1429',
                        'AA1453',
                        'AA1460',
                        'AA1428',
                        'AA1514'] #sotet kek kaosz csopi
        
    elif basic_group in ['gray','grey']:
        needed_cells = ['AA1310',
                        'AA1069',
                        'AA1030',
                        'AA0515',
                        'AA0503',
                        'AA1511',
                        'AA1521'] # szurke csopi (talan 1310 es 1069 sotetkek)
        
    elif basic_group == 'yellow':
          needed_cells = ['AA1326',
                          'AA1452',
                          'AA1455',
                          'AA1456',
                          'AA1458',
                          'AA1461',
                          'AA1536',
                          'AA1525',
                          'AA1524',
                          'AA1523',
                          'AA1522',
                          'AA1519',
                          'AA1517',
                          'AA1513']
                          
          
    elif basic_group == 'motor':
        needed_cells = ['AA1348',
                        'AA1362',
                        'AA1361',
                        'AA1359',
                        'AA0996',
                        'AA1405',
                        'AA1331',
                        'AA1196',
                        'AA0434'] # motoros magba projektalok
        
    elif basic_group == 'all':
        needed_cells = original_data['cell_names']
    elif basic_group in np.unique(original_data['cell_juci_clusters']):
        needed_cells = np.asarray(original_data['cell_names'])[np.where(np.asarray(original_data['cell_juci_clusters']) ==basic_group)[0]]
    elif basic_group == 'soma location':
        #soma_locations_needed = ['Intermediate reticular nucleus', 'Medullary reticular nucleus', 'Medullary reticular nucleus, dorsal part', 'Parvicellular reticular nucleus']
        cell_idxs = list()
        for soma_location_now in soma_locations_needed:
            cell_idxs.append(np.where(np.asarray(original_data['soma_locations'])==soma_location_now)[0])
        cell_idxs = np.concatenate(cell_idxs)
        needed_cells = list(np.asarray(original_data['cell_names'])[cell_idxs])
    elif basic_group == 'axon projection':
        metric_to_use = 'allen_axon_matrix'#'allen_axon_end_points_matrix' #'allen_axon_matrix'
        #axon_projections_needed = ['VII','XII','AMB','V']
        cell_idxs = list()
        for axon_projection_now in axon_projections_needed:
            if (allen_df['acronym']==axon_projection_now).sum()>0:
                projection_idx_ipsi = (allen_df['acronym']==axon_projection_now).argmax()
                projection_idx_contra = projection_idx_ipsi + len(allen_df)
            elif (allen_df['name']==axon_projection_now).sum()>0:
                projection_idx_ipsi = (allen_df['name']==axon_projection_now).argmax()
                projection_idx_contra = projection_idx_ipsi + len(allen_df)
            else:
                print('unknown brain region {}'.format(axon_projection_now))
            cell_idxs.append(np.where(original_data[metric_to_use][:,projection_idx_ipsi]>0)[0])
            cell_idxs.append(np.where(original_data[metric_to_use][:,projection_idx_contra]>0)[0])
        cell_idxs = np.unique(np.concatenate(cell_idxs))
        needed_cells = list(np.asarray(original_data['cell_names'])[cell_idxs])
    else:
        print('unrecognized group')
        needed_cells = None
    
    return needed_cells

def subselect_cells(allen_df,original_data,basic_group,soma_locations_needed=None,axon_projections_needed=None):

    needed_cells = generate_cell_list(allen_df,original_data,basic_group,soma_locations_needed=None,axon_projections_needed=None)
    
    cell_idx_needed = list()
    for cell_id in needed_cells:
        cell_idx_needed.append(np.argmax(np.asarray(original_data['cell_names'])==cell_id))
    cell_idx_needed = np.asarray(cell_idx_needed)
    data_now = {}
    for variable_name in original_data.keys():
        if not 'header' in variable_name:
            try:
                #exec('{} = np.asarray(original_data['"'{}'"'])[cell_idx_needed].copy()'.format(variable_name,variable_name))
                data_now[variable_name] = np.asarray(original_data[variable_name])[cell_idx_needed].copy()
            except:
                if len(original_data[variable_name])>0:
                    print('could not subselect {}'.format(variable_name))
        else:
            #exec('{} = np.asarray(original_data['"'{}'"']).copy()'.format(variable_name,variable_name))
            data_now[variable_name] = np.asarray(original_data[variable_name]).copy()
            print('{} skipped'.format(variable_name))
    return data_now

def get_allen_annotation(x,y,z,mirror_cell = False,ccf_version = '3.0',allen_annotation_volume=None):#,ccf_version = '3.0'
    midline_coordinate = 5700
    offset = -.5
    if ccf_version =='2.5':
        if mirror_cell:# and x<midline_coordinate:
            x=2*midline_coordinate -x
        allen_idxs= np.array((np.array([z,y,x])/10+offset).round(),int)
    else:
        if mirror_cell:# and z<midline_coordinate:
            z=2*midline_coordinate -z
        allen_idxs= np.array((np.array([x,y,z])/10+offset).round(),int)
    allen_id = allen_annotation_volume[allen_idxs[0],allen_idxs[1],allen_idxs[2]]            
    return allen_id

        
def remove_nan_from_spreadsheet(spreadsheet):
    #%
    for col in spreadsheet.items(): 
        head = col[0]
        col = col[1]
        try:
            for item in col.iteritems():
                if type(item[1])== float and np.isnan(item[1]):
                    spreadsheet.loc[item[0],head] = spreadsheet[head][item[0]-1]
        except:
            print('could not remove all nans from {}'.format(head))
            
                #%
    return spreadsheet   



    


#missing_cells = list()
#%% generate table for clustering
# ez csinalja a munka nagy reszet, kibanyassza a sejtek adatait a json fileokbol es a google tablazatbol.. lassu (ctrl+enter)
# itt a kockak meretet allithatod az axon es dendrit eloszlasos voxel analizishez
#columns = ['Neuron ID']
#big_spreadsheet = pd.DataFrame(columns = ) 
# ----------------------------- PARAMETERS -----------------------------♥
#ccf_disagreements_list = list()
#node_num_list = list()

def analyze_json_files(allen_df,parameters_dict):
    #% clean up allen_df database

    allen_annotation_volume = ccf.load_ccf_volume(parameters_dict)
    #%
    if parameters_dict['assign_spinal_cord_to_ccf'] and 'Spinal cord' not in allen_df['name'].values:
        structureId = np.max(allen_df['structureId'])+1
        spinal_cord_dict = {'acronym':'SPC',
                            'structureId':structureId,
                            'parentStructureId':8,
                            'depth':3,
                            'name':'Spinal cord',
                            'structureIdPath':'/997/8/{}'.format(structureId),
                            'volume':50000.,
                            'final_nucleus':True}
        allen_df = allen_df.append(spinal_cord_dict, ignore_index=True)
    #%
    stepsize_list_for_axon_endings = list()
    stepsize_list_for_somas = list()
    endpoint_num_outside_grey_matter = 0
    soma_num_outside_grey_matter = 0
    endpoint_num_outside_brain  = 0
    soma_num_outside_brain = 0
    midline_coordinate_x = 5700#AllenUtils.brainCenter().getX()
    
    resolution_coarse = 500 # microns
    normalizing_function_coarse = 'sum' #sum/max
    resolution_fine = 100# microns
    normalizing_function_fine = 'sum'#sum/max
    
    resolution_dendrite_coarse = 500# microns
    
    resolution_dendrite_fine = 100# microns
    
    sigma_step = 1 #sigma / resolution
    show_reconstruction = False
    random_displacement_n = 1 #number of iterations
    random_displacement_distance = 100 #microns
    
    
    dendrite_length_normalize_percentiles = [0,100]
    dendrite_length_cut_outliers = False
    dendrite_branch_point_numbers_normalize_percentiles = [0,100]
    dendrite_primary_branch_numbers_normalize_percentiles = [0,100]
    dendrite_branch_numbers_normalize_percentiles = [0,100]
    dendrite_end_point_numbers_normalize_percentiles = [0,100]
    
    axon_length_normalize_percentiles = [0,100]
    axon_length_cut_outliers = False
    axon_branch_point_numbers_normalize_percentiles = [0,100]
    axon_branch_numbers_normalize_percentiles = [0,100]
    axon_end_point_numbers_normalize_percentiles = [0,100]
    
    # ----------------------------- PARAMETERS -----------------------------♥
    displacements_x = np.random.uniform(-1*random_displacement_distance,random_displacement_distance,random_displacement_n)
    displacements_y = np.random.uniform(-1*random_displacement_distance,random_displacement_distance,random_displacement_n)
    displacements_z = np.random.uniform(-1*random_displacement_distance,random_displacement_distance,random_displacement_n)
    
    
    volumes_coarse = list()
    volumes_fine = list()
    volumes_dendrite_coarse = list()
    volumes_dendrite_fine = list()
    allen_axon_list = list()
    allen_dendrite_list = list()
    allen_soma_list = list()
    allen_soma_distribution_list = list()
    allen_axon_branch_point_list = list()
    allen_axon_end_point_list = list()
    allen_dendrite_branch_point_list = list()
    allen_dendrite_end_point_list = list()
    allen_axon_list_original = list()
    allen_dendrite_list_original = list()
    allen_soma_list_original = list()
    allen_soma_distribution_list_original = list()
    allen_axon_branch_point_list_original = list()
    allen_axon_end_point_list_original = list()
    allen_dendrite_branch_point_list_original = list()
    allen_dendrite_end_point_list_original = list()
    
    dendrite_bounding_ball_radius = []
    dendrite_bounding_ball_coordinates = []
    
    
    cell_names = list()#♥ <3 
    cell_juci_clusters = list()
    soma_coordinates = {'X':list(),
                        'Y':list(),
                        'Z':list()}
    soma_locations = list()
    target_areas = list()
    target_locations = list()
    if parameters_dict['project'] == 'medulla':
        soma_area_in_the_medulla = list() # sensory vs motor related
    subtle_features = {'Soma shape':list(),
                       'Axon origin point':list(),
                       'Bouton types':list(),
                       'Dendritic spines':list()}
    dendrite_lengths = list()
    axon_lengths=list()
    axon_branch_point_numbers = list()
    axon_branch_numbers = list()
    axon_end_point_numbers = list() # endpoint number with identified ontology
    axon_end_point_numbers_total = list() # endpoint number based on the json file
    
    dendrite_branch_point_numbers = list()
    dendrite_branch_numbers = list()
    dendrite_primary_branch_numbers = list()
    dendrite_end_point_numbers = list()
    
    xscale_coarse = np.arange(1000,10000,resolution_coarse)
    yscale_coarse = np.arange(500,8000,resolution_coarse)
    zscale_coarse = np.arange(6000,14000,resolution_coarse)
    
    xscale_fine = np.arange(1000,10000,resolution_fine)
    yscale_fine = np.arange(500,8000,resolution_fine)
    zscale_fine = np.arange(6000,14000,resolution_fine)
    
    xscale_dendrite_coarse = np.arange(1000,10000,resolution_dendrite_coarse)
    yscale_dendrite_coarse = np.arange(500,8000,resolution_dendrite_coarse)
    zscale_dendrite_coarse = np.arange(6000,14000,resolution_dendrite_coarse)
    
    xscale_dendrite_fine = np.arange(1000,10000,resolution_dendrite_fine)
    yscale_dendrite_fine = np.arange(500,8000,resolution_dendrite_fine)
    zscale_dendrite_fine = np.arange(6000,14000,resolution_dendrite_fine)
    #sigma = 1
    
    #%
    
    
    
    
    target_brain_part_names = os.listdir(parameters_dict['metadata_dir'])
    X_edges = [np.nan,np.nan]
    Y_edges = [np.nan,np.nan]
    Z_edges = [np.nan,np.nan]
    #%
    
    for target_brain_part_name in target_brain_part_names:
        if '_DO NOT ANALYSE' in target_brain_part_name:
            continue
        try:
            #print(target_brain_part_name)
            metadata_spreadsheet = pd.read_csv(os.path.join(parameters_dict['metadata_dir'],target_brain_part_name))
            metadata_spreadsheet = remove_nan_from_spreadsheet(metadata_spreadsheet)
            #%
            neuron_IDs = metadata_spreadsheet['Neuron ID'].unique()
        except: # there are no neurons in this spreadsheet, probably not apropriate
            continue
        #%
        for neuron_ID in neuron_IDs:
         #%
            #ccf_disagreements_now = 0
            cell_spreadsheet = metadata_spreadsheet.loc[metadata_spreadsheet['Neuron ID']==neuron_ID]
            target_areas.append(cell_spreadsheet['Target area'])
            target_locations.append(cell_spreadsheet['Target location'])
            #%
            #%
            soma_locations.append(cell_spreadsheet['Soma location'].unique()[0])
            
            
            
            if parameters_dict['project'] == 'medulla':
                soma_area_in_the_medulla.append(cell_spreadsheet['Soma area in the medulla'].unique()[0])
            subtle_features['Soma shape'].append(cell_spreadsheet['Soma shape'].unique()[0])
            subtle_features['Axon origin point'].append(cell_spreadsheet['Axon origin point'].unique()[0])
            bouton_types = list()
            bouton_types_raw = cell_spreadsheet['Bouton types'].unique()
            for bouton_type in bouton_types_raw:
                if '[' in bouton_type:
                    bouton_type = bouton_type.replace('[','[\'')
                    bouton_type = bouton_type.replace(',','\',\'')
                    bouton_type = bouton_type.replace(']','\']')
                    bouton_types.extend(eval(bouton_type))
                    
                else:
                    bouton_types.append(bouton_type)
                    #%
            bouton_types = np.unique(bouton_types ).tolist()
            bouton_types_final = list()
            for bouton in bouton_types:
                bouton_types_final.append(bouton.strip())
                
            subtle_features['Bouton types'].append(bouton_types_final)
            if cell_spreadsheet['Dendritic spines'].unique()[0].lower()=='no':
                subtle_features['Dendritic spines'].append(False)
            elif cell_spreadsheet['Dendritic spines'].unique()[0].lower()=='yes':
                subtle_features['Dendritic spines'].append(True)
            else:
                print('ERROR: Dendritic spines is not yes/no for {}'.format(neuron_ID) )
            #%
            #%
            cell_names.append(neuron_ID)
            cell_juci_clusters.append(target_brain_part_name[:-4])
            
            print(neuron_ID)
            if show_reconstruction:
                
                fig = plt.figure()
                ax_3d = fig.add_subplot(331, projection='3d')
                ax_yz = fig.add_subplot(334)
                ax_xz = fig.add_subplot(335)
                ax_xy = fig.add_subplot(336)
                ax_yz_fine = fig.add_subplot(337)
                ax_xz_fine = fig.add_subplot(338)
                ax_xy_fine = fig.add_subplot(339)
                
                fig_dendrite = plt.figure()
                ax_dendrite_3d = fig_dendrite.add_subplot(331, projection='3d')
                ax_dendrite_yz = fig_dendrite.add_subplot(334)
                ax_dendrite_xz = fig_dendrite.add_subplot(335)
                ax_dendrite_xy = fig_dendrite.add_subplot(336)
                ax_dendrite_yz_fine = fig_dendrite.add_subplot(337)
                ax_dendrite_xz_fine = fig_dendrite.add_subplot(338)
                ax_dendrite_xy_fine = fig_dendrite.add_subplot(339)
            sum_axon_length = 0 
            sum_dendrite_length = 0
            sum_axon_branch_point_num = 0
            sum_axon_end_point_num = 0
            sum_dendrite_branch_point_num = 0
            sum_dendrite_end_point_num = 0
            volume_coarse = np.zeros([len(xscale_coarse),len(yscale_coarse),len(zscale_coarse)])
            volume_fine = np.zeros([len(xscale_fine),len(yscale_fine),len(zscale_fine)])
            volume_dendrite_coarse = np.zeros([len(xscale_dendrite_coarse),len(yscale_dendrite_coarse),len(zscale_dendrite_coarse)])
            volume_dendrite_fine = np.zeros([len(xscale_dendrite_fine),len(yscale_dendrite_fine),len(zscale_dendrite_fine)])
            allen_axon_ipsi = np.zeros(len(allen_df))
            allen_axon_contra = np.zeros(len(allen_df))
            allen_dendrite = np.zeros(len(allen_df))
            allen_soma = np.zeros(len(allen_df))
            allen_soma_distribution = np.zeros([len(allen_df),random_displacement_n])
            allen_axon_branch_point_ipsi = np.zeros(len(allen_df))
            allen_axon_branch_point_contra = np.zeros(len(allen_df))
            allen_axon_end_point_ipsi = np.zeros(len(allen_df))
            allen_axon_end_point_contra = np.zeros(len(allen_df))
            allen_dendrite_branch_point_ipsi = np.zeros(len(allen_df))
            allen_dendrite_branch_point_contra = np.zeros(len(allen_df))
            allen_dendrite_end_point_ipsi = np.zeros(len(allen_df))
            allen_dendrite_end_point_contra = np.zeros(len(allen_df))
            
            xvals = X_edges
            yvals = Y_edges
            zvals = Z_edges
            jsonfile_with_path=os.path.join(parameters_dict['json_dir'],neuron_ID+'.json')
    # =============================================================================
    #         if not os.path.isfile(jsonfile_with_path):
    #             print('----------------------------------------------'+neuron_ID)
    #             missing_cells.append(neuron_ID)
    #         continue
    # =============================================================================
            axon = snt.Tree(jsonfile_with_path, "axon")
    # =============================================================================
    #         node_num_list.append(len(list(axon.getNodes().toArray())))
    #         continue
    # =============================================================================
            try:
                dendrite = snt.Tree(jsonfile_with_path, "dendrite")
            except:
                print('no dendrite found for  {}, using the dendrite of the last cell...'.format(neuron_ID))
    # =============================================================================
    #         break
    #     break
    # #%%
    # =============================================================================
            if parameters_dict['ccf_version'] == '3.0':
                if axon.getRoot().z>midline_coordinate_x:
                    mirror_cell = True
                else:
                    mirror_cell = False
            elif parameters_dict['ccf_version'] == '2.5':
                if axon.getRoot().x>midline_coordinate_x:
                    mirror_cell = True
                else:
                    mirror_cell = False
           ### mirror_cell = False
            
            if parameters_dict['ccf_version'] == '3.0':
                if mirror_cell:
                    soma_coordinates['Z'].append(2*midline_coordinate_x -axon.getRoot().z)
                else:
                    soma_coordinates['Z'].append(axon.getRoot().z)
                soma_coordinates['Y'].append(axon.getRoot().y)
                soma_coordinates['X'].append(axon.getRoot().x)
            elif parameters_dict['ccf_version'] == '2.5':
                if mirror_cell:
                    soma_coordinates['X'].append(2*midline_coordinate_x -axon.getRoot().x)
                else:
                    soma_coordinates['X'].append(axon.getRoot().x)
                soma_coordinates['Y'].append(axon.getRoot().y)
                soma_coordinates['Z'].append(axon.getRoot().z)
            try:
                #node_id = axon.getRoot().annotation.id()  
    # =============================================================================
    #             node_id = axon.getRoot().getAnnotation().id()#id
    #             node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
    # =============================================================================
                
                somanode =   axon.getRoot()
                somanode.x = soma_coordinates['X'][-1]
                somanode.y = soma_coordinates['Y'][-1]
                somanode.z = soma_coordinates['Z'][-1]
                if somanode.getAnnotation() == None:
                    node_id = 997
                    node_ids = [997]
                else:
                    node_id = somanode.getAnnotation().id()  
                    node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                final_nucleus = allen_df.loc[allen_df['structureId']==node_id,'final_nucleus'].values[0]
                if parameters_dict['allocate_soma_in_white_matter_to_closest_grey'] and not final_nucleus:#(8 not in node_ids or any(allen_df['parentStructureId']==node_id)): # 8 is the mother of all gray matters, and it should not be a parent # 8 not in node_ids:#
                    x = somanode.x
                    y = somanode.y
                    z = somanode.z
                    for stepsize in parameters_dict['correction_steps']:
                        for dx,dy,dz in zip([1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1],[0,0,1,-1,0,0,1,-1,-1,1,1,-1,-1,1],[0,0,0,0,1,-1,-1,-1,-1,-1,1,1,1,1]):
                            try:
                                node_id = get_allen_annotation(x+dx*stepsize,y+dy*stepsize,z+dz*stepsize,mirror_cell,parameters_dict['ccf_version'],allen_annotation_volume)
                                node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                            except:
                                node_id = 997
                                node_ids = [997]
                            final_nucleus = allen_df.loc[allen_df['structureId']==node_id,'final_nucleus'].values[0]
                            if final_nucleus:#8 in node_ids and not any(allen_df['parentStructureId']==node_id):
                                #print('found grey matter with a step of {}'.format(stepsize))
                                stepsize_list_for_somas.append(stepsize)
                                break
                        if final_nucleus:#8 in node_ids and not any(allen_df['parentStructureId']==node_id):
                            #print('node found')
                            break
                    if not final_nucleus:#8 not in node_ids or any(allen_df['parentStructureId']==node_id): #returning to original point if couldn't break out of white matter
                        node_id = somanode.getAnnotation().id() 
                        node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                        soma_num_outside_grey_matter += 1
                
                
                
                for node_id in node_ids:
                    allen_soma[(allen_df['structureId']==node_id).argmax()] += 1
            except:
                print('no annotation found for soma: {}'.format(neuron_ID))
                soma_num_outside_brain += 1
                #time.sleep(1000)
            #% displacements of soma
            for idx,(dx,dy,dz) in enumerate(zip(displacements_x,displacements_y,displacements_z)):
                try:
                    node_id = get_allen_annotation(soma_coordinates['X'][-1]+dx,soma_coordinates['Y'][-1]+dy,soma_coordinates['Z'][-1]+dz,mirror_cell,parameters_dict['ccf_version'],allen_annotation_volume)
                    node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                    for node_id in node_ids:
                        allen_soma_distribution[(allen_df['structureId']==node_id).argmax(),idx] += 1
                except:
                    pass
            allen_soma_distribution = np.mean(allen_soma_distribution,1)
            #%
            
            axonstats = snt.TreeStatistics(axon)
            branches = axonstats.getBranches().toArray() #branches
            branchlengths_local = list()
            branchlengths_nonlocal = list()
    
    
            axonendings =  axon.getGraph().getTips().toArray()
      
            
            if parameters_dict['use_probabilistic_axon_ending_map']:
                #%
                allen_axon_end_point_ipsi = np.zeros(len(allen_df))
                allen_axon_end_point_contra = np.zeros(len(allen_df))
                gaussian_kernel_dim = int(parameters_dict['probabilistic_axon_ending_map_sigma']*5/10+1)
                gaussian_kernel_offset = int(gaussian_kernel_dim/2)
                gaussian_kernel = np.zeros([gaussian_kernel_dim,gaussian_kernel_dim,gaussian_kernel_dim],float)
                gaussian_kernel[int(gaussian_kernel_dim/2),int(gaussian_kernel_dim/2),int(gaussian_kernel_dim/2)] = 1
                gaussian_kernel = scipy.ndimage.gaussian_filter(gaussian_kernel,parameters_dict['probabilistic_axon_ending_map_sigma']/10)
                probabilistic_axon_endings_filt = np.zeros_like(allen_annotation_volume,float)
                for axonending in axonendings :
                    if parameters_dict['ccf_version'] == '3.0': # it is an ugly hack but I cannot think it through.. this is a mess, thanks for switching up the axes..
                        z = axonending.z
                        x = axonending.x
                        axonending.z = x
                        axonending.x = z
                    if mirror_cell:
                        axonending.x  = 2*midline_coordinate_x -axonending.x
    # =============================================================================
    #                 if parameters_dict['ccf_version'] == '3.0':
    #                     if mirror_cell:
    #                         axonending.z  = 2*midline_coordinate_x -axonending.z
    #                 elif parameters_dict['ccf_version'] == '2.5':
    #                     if mirror_cell:
    #                         axonending.x  = 2*midline_coordinate_x -axonending.x
    # =============================================================================
                    allen_idxs= np.array((np.array([axonending.z,axonending.y,axonending.x])/10).round(),int)
                    #if all(allen_idxs>=gaussian_kernel_offset) and all(np.asarray(allen_annotation_volume.shape)>allen_idxs+gaussian_kernel_offset):
                    kernel_offset_x1 = np.min([gaussian_kernel_offset+1,allen_idxs[0]])
                    kernel_offset_x2 =np.min([gaussian_kernel_offset,allen_annotation_volume.shape[0]-allen_idxs[0]])
                    kernel_offset_y1 = np.min([gaussian_kernel_offset+1,allen_idxs[1]])
                    kernel_offset_y2 =np.min([gaussian_kernel_offset,allen_annotation_volume.shape[1]-allen_idxs[1]])
                    kernel_offset_z1 = np.min([gaussian_kernel_offset+1,allen_idxs[2]])
                    kernel_offset_z2 =np.min([gaussian_kernel_offset,allen_annotation_volume.shape[2]-allen_idxs[2]])
                    if any(np.asarray([kernel_offset_x1,kernel_offset_x2,kernel_offset_y1,kernel_offset_y2,kernel_offset_z1,kernel_offset_z2])<0):
                        print('axon ending outside allen map, skipping: {}'.format(allen_idxs))
                        #time.sleep(10000)
                        continue
                    probabilistic_axon_endings_filt[allen_idxs[0]-kernel_offset_x1:allen_idxs[0]+kernel_offset_x2,allen_idxs[1]-kernel_offset_y1:allen_idxs[1]+kernel_offset_y2,allen_idxs[2]-kernel_offset_z1:allen_idxs[2]+kernel_offset_z2]   += gaussian_kernel[gaussian_kernel_offset+1-kernel_offset_x1:gaussian_kernel_offset+kernel_offset_x2+1,gaussian_kernel_offset+1-kernel_offset_y1:gaussian_kernel_offset+kernel_offset_y2+1,gaussian_kernel_offset+1-kernel_offset_z1:gaussian_kernel_offset+kernel_offset_z2+1]
                    sum_axon_end_point_num += 1
    
                nonzero_idxes = np.where(probabilistic_axon_endings_filt>0)
                #%
                axon_ending_values = list()
                for x,y,z in zip(nonzero_idxes[0],nonzero_idxes[1],nonzero_idxes[2]):
                    axon_ending_value = probabilistic_axon_endings_filt[x,y,z]
                    node_id = allen_annotation_volume[x,y,z]
                    if node_id ==0:
                        continue
                    try:
                        node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                    except:
                        print('allen ID {} not found in allen dataframe'.format(node_id))
                        continue
    # =============================================================================
    #                 if parameters_dict['ccf_version'] == '3.0':
    #                     lateral_coord = z
    #                 elif parameters_dict['ccf_version'] == '2.5':
    #                     lateral_coord = x
    # =============================================================================
                    for node_id in node_ids:
                        if x>midline_coordinate_x/10:
                            allen_axon_end_point_contra[(allen_df['structureId']==node_id).argmax()] += axon_ending_value
                        else:
                            allen_axon_end_point_ipsi[(allen_df['structureId']==node_id).argmax()] += axon_ending_value
                            #%
            else:
                for axonending in axonendings :
                    try:
                        #%
                        if parameters_dict['ccf_version'] == '3.0':
                            if mirror_cell:
                                axonending.z  = 2*midline_coordinate_x -axonending.z
                            axonending_lateral = axonending.z
                        elif parameters_dict['ccf_version'] == '2.5':
                            if mirror_cell:
                                axonending.x  = 2*midline_coordinate_x -axonending.x                            
                            axonending_lateral = axonending.x
                           
                        #node_id = axonending.annotation.id()  
                        if axonending.getAnnotation() == None:
                            node_id = 997
                            node_ids = [997]
                        else:
                            node_id = axonending.getAnnotation().id()  
                            node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                        final_nucleus = allen_df.loc[allen_df['structureId']==node_id,'final_nucleus'].values[0]
                        if parameters_dict['allocate_axon_ending_in_white_matter_to_closest_grey'] and not final_nucleus:#(8 not in node_ids or any(allen_df['parentStructureId']==node_id)): # 8 is the mother of all gray matters, and it should not be a parent # 8 not in node_ids:#
                        #%
    # =============================================================================
    #                         mirroragain=False
    #                         while 8 not in node_ids:#this loops goes back along the axon to reach gray matter
    #                             #print('moving axon ending back')
    #                             mirroragain = True
    #                             axonending = axonending.getPreviousPoint()
    #                             if axonending.getAnnotation() == None:
    #                                 node_id = 997
    #                                 node_ids = [997]
    #                             else:
    #                                 node_id = axonending.getAnnotation().id()  
    #                                 node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
    #                         if mirroragain:
    #                             if parameters_dict['ccf_version'] == '3.0':
    #                                 if mirror_cell:
    #                                     axonending.z  = 2*midline_coordinate_x -axonending.z
    #                                 axonending_lateral = axonending.z
    #                             elif parameters_dict['ccf_version'] == '2.5':
    #                                 if mirror_cell:
    #                                     axonending.x  = 2*midline_coordinate_x -axonending.x                            
    #                                 axonending_lateral = axonending.x
    # =============================================================================
                            #%
                            x = axonending.x
                            y = axonending.y
                            z = axonending.z
                            for stepsize in parameters_dict['correction_steps']:
                                for dx,dy,dz in zip([1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1],[0,0,1,-1,0,0,1,-1,-1,1,1,-1,-1,1],[0,0,0,0,1,-1,-1,-1,-1,-1,1,1,1,1]):
                                    try:
                                        node_id = get_allen_annotation(x+dx*stepsize,y+dy*stepsize,z+dz*stepsize,mirror_cell,parameters_dict['ccf_version'],allen_annotation_volume)
                                        node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                                    except:
                                        node_id = 997
                                        node_ids = [997]
                                    final_nucleus = allen_df.loc[allen_df['structureId']==node_id,'final_nucleus'].values[0]
                                    if final_nucleus:#8 in node_ids and not any(allen_df['parentStructureId']==node_id):
                                        #print('found grey matter with a step of {}'.format(stepsize))
                                        stepsize_list_for_axon_endings.append(stepsize)
                                        break
                                if final_nucleus:#8 in node_ids and not any(allen_df['parentStructureId']==node_id):
                                    #print('node found')
                                    break
                                if parameters_dict['assign_spinal_cord_to_ccf'] and node_id == 997 and axonending.y>4700 and (axonending.x>13000 or axonending.y>7300):# ASSIGN to spinal cord
                                    node_id = allen_df.loc[allen_df['name']=='Spinal cord']['structureId'].values[0]
                                    node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                                    final_nucleus  = True
                                    break
                                
                            if not final_nucleus:#8 not in node_ids or any(allen_df['parentStructureId']==node_id): #returning to original point if couldn't break out of white matter
                                node_id = axonending.getAnnotation().id() 
                                node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                                endpoint_num_outside_grey_matter += 1
                        for node_id in node_ids:
                            if axonending_lateral>midline_coordinate_x:
                                allen_axon_end_point_contra[(allen_df['structureId']==node_id).argmax()] += 1
                            else:
                                allen_axon_end_point_ipsi[(allen_df['structureId']==node_id).argmax()] += 1
                        sum_axon_end_point_num += 1
                        #%
                    except:
                        #time.sleep(1000)
                        #print('axon ending outside brain: {}'.format([x,y,z]))
                        endpoint_num_outside_brain += 1
                        pass # no brain region for branch point
             #% 
            axonbranchpoints =  axon.getGraph().getBPs().toArray()    
            for axonbranchpoint in axonbranchpoints :
                try:
                    if parameters_dict['ccf_version'] == '3.0':
                        if mirror_cell:
                            axonbranchpoint.z  = 2*midline_coordinate_x -axonbranchpoint.z
                        axonbranchpoint_lateral = axonbranchpoint.z 
                    elif parameters_dict['ccf_version'] == '2.5':
                        if mirror_cell:
                            axonbranchpoint.x  = 2*midline_coordinate_x -axonbranchpoint.x
                        axonbranchpoint_lateral = axonbranchpoint.x 
                    #node_id = axonbranchpoint.annotation.id()  
                    node_id = axonbranchpoint.getAnnotation().id()  
                    node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                    for node_id in node_ids:
                        if axonbranchpoint_lateral>midline_coordinate_x:
                            allen_axon_branch_point_contra[(allen_df['structureId']==node_id).argmax()] += 1
                        else:
                            allen_axon_branch_point_ipsi[(allen_df['structureId']==node_id).argmax()] += 1
                    sum_axon_branch_point_num += 1
                except:
                    pass # no brain region for branch point
                
            for branch in branches: # axon loop
                branch_x = list()
                branch_y = list()
                branch_z = list()
                nodeids = list()
                node_idx = 0
                while True:
                    node_idx += 1
                    try:
                        
                        node = branch.getNode(node_idx)
                        node_prev = branch.getNode(node_idx-1)
        
                        node_x = node.x
                        node_y = node.y
                        node_z = node.z
                        
                        node_prev_x = node_prev.x
                        node_prev_y = node_prev.y
                        node_prev_z = node_prev.z
                        if parameters_dict['ccf_version'] == '3.0':
                            if mirror_cell:
                                node_z  = 2*midline_coordinate_x -node_z
                                node_prev_z  = 2*midline_coordinate_x - node_prev_z
                            node_lateral = node_z
                        elif parameters_dict['ccf_version'] == '2.5':
                            if mirror_cell:
                                node_x  = 2*midline_coordinate_x -node_x
                                node_prev_x  = 2*midline_coordinate_x - node_prev_x
                            node_lateral = node_x
                        
                        xvals.append(node_x)
                        yvals.append(node_y)
                        zvals.append(node_z)
                        if show_reconstruction:
                            branch_x.append(node_x)
                            branch_y.append(node_y)
                            branch_z.append(node_z)
                        
                        nodeidx_coarse = [np.argmax(node_x < xscale_coarse),np.argmax(node_y < yscale_coarse),np.argmax(node_z < zscale_coarse)]
                        nodeidx_coarse_prev = [np.argmax(node_prev_x < xscale_coarse),np.argmax(node_prev_y < yscale_coarse),np.argmax(node_prev_z < zscale_coarse)]
                        if nodeidx_coarse == nodeidx_coarse_prev:
                            volume_coarse[nodeidx_coarse[0]][nodeidx_coarse[1]][nodeidx_coarse[2]] += node.distanceTo(node_prev)
                            
                        nodeidx_fine = [np.argmax(node_x < xscale_fine),np.argmax(node_y < yscale_fine),np.argmax(node_z < zscale_fine)]
                        nodeidx_fine_prev = [np.argmax(node_prev_x < xscale_fine),np.argmax(node_prev_y < yscale_fine),np.argmax(node_prev_z < zscale_fine)]
                        if nodeidx_fine == nodeidx_fine_prev:
                            volume_fine[nodeidx_fine[0]][nodeidx_fine[1]][nodeidx_fine[2]] += node.distanceTo(node_prev)
                        
                #%    
                        try:
                            #node_id = node.annotation.id()  
                            node_id = node.getAnnotation().id()  
                            
                            
    # =============================================================================
    #                         node_id_measured = get_allen_annotation(node_x,node_y,node_z,mirror_cell,parameters_dict['ccf_version'],allen_annotation_volume)
    #                         if not node_id ==node_id_measured:
    #                             print('neuron {} x:{},y:{}z:{} -- node_id : {} => {}'.format(neuron_ID,int(node_x),int(node_y),int(node_z),node_id,node_id_measured))
    #                            # time.sleep(10000)
    #                             ccf_disagreements_now+=1
    # =============================================================================
                            
                            
                            node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                            #node_prev_id = node_prev.annotation.id()
                            node_prev_id = node_prev.getAnnotation().id()
                            node_prev_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_prev_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                            for node_id in node_ids:
                                if node_id in node_prev_ids:
                                    if node_lateral>midline_coordinate_x:
                                        allen_axon_contra[(allen_df['structureId']==node_id).argmax()] += node.distanceTo(node_prev)
                                    else:
                                        allen_axon_ipsi[(allen_df['structureId']==node_id).argmax()] += node.distanceTo(node_prev)
                            sum_axon_length += node.distanceTo(node_prev)
                        except:
                            pass
                        
                       #%
                    except: # this is a branch ending
                        break
                if show_reconstruction:
                    ax_3d.plot(branch_x,branch_y,branch_z)       
                    
            dendritestats = snt.TreeStatistics(dendrite)
            branches = dendritestats.getBranches().toArray()   
            dendriteendings =  dendrite.getGraph().getTips().toArray()
            dendrite_endpoints = []
            for dendriteending in dendriteendings :
                dendrite_endpoints.append([dendriteending.x,dendriteending.y,dendriteending.z])
                try:
                    if parameters_dict['ccf_version'] == '3.0':
                        if mirror_cell:
                            dendriteending.z  = 2*midline_coordinate_x -dendriteending.z
                        dendriteending_lateral = dendriteending.z
                    elif parameters_dict['ccf_version'] == '2.5':
                        if mirror_cell:
                            dendriteending.x  = 2*midline_coordinate_x -dendriteending.x
                        dendriteending_lateral = dendriteending.x
                    #node_id = dendriteending.annotation.id()  
                    node_id = dendriteending.getAnnotation().id() 
                    node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                    for node_id in node_ids:
                        if dendriteending_lateral>midline_coordinate_x:
                            allen_dendrite_end_point_contra[(allen_df['structureId']==node_id).argmax()] += 1
                        else:
                            allen_dendrite_end_point_ipsi[(allen_df['structureId']==node_id).argmax()] += 1
                    sum_dendrite_end_point_num += 1
                except:
                    pass # no brain region for branch point
            C, r2 = miniball.get_bounding_ball(np.asarray(dendrite_endpoints))
            dendrite_bounding_ball_radius.append(np.sqrt(r2))
            dendrite_bounding_ball_coordinates.append(C)
            dendritebranchpoints =  dendrite.getGraph().getBPs().toArray()    
            for dendritebranchpoint in dendritebranchpoints :
                try:
                    if parameters_dict['ccf_version'] == '3.0':
                        if mirror_cell:
                            dendritebranchpoint.z  = 2*midline_coordinate_x -dendritebranchpoint.z
                        dendritebranchpoint_lateral = dendritebranchpoint.z
                    elif parameters_dict['ccf_version'] == '2.5':
                        if mirror_cell:
                            dendritebranchpoint.x  = 2*midline_coordinate_x -dendritebranchpoint.x
                        dendritebranchpoint_lateral=  dendritebranchpoint.x
                    #node_id = dendritebranchpoint.annotation.id()  
                    node_id = dendritebranchpoint.getAnnotation().id()  
                    node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                    for node_id in node_ids:
                        if dendritebranchpoint_lateral>midline_coordinate_x:
                            allen_dendrite_branch_point_contra[(allen_df['structureId']==node_id).argmax()] += 1
                        else:
                            allen_dendrite_branch_point_ipsi[(allen_df['structureId']==node_id).argmax()] += 1
                    sum_dendrite_branch_point_num += 1
                except:
                    pass # no brain region for branch point
                    
            for branch in branches: # dendrite loop
                branch_x = list()
                branch_y = list()
                branch_z = list()
                nodeids = list()
                node_idx = 0
                while True:
                    node_idx += 1
                    try:
                        
                        node = branch.getNode(node_idx)
                        node_prev = branch.getNode(node_idx-1)
    
                        node_x = node.x
                        node_y = node.y
                        node_z = node.z
                        
                        node_prev_x = node_prev.x
                        node_prev_y = node_prev.y
                        node_prev_z = node_prev.z
                        
                        if parameters_dict['ccf_version'] == '3.0':
                            if mirror_cell:
                                node_z  = 2*midline_coordinate_x -node_z
                                node_prev_z  = 2*midline_coordinate_x - node_prev_z
                        elif parameters_dict['ccf_version'] == '2.5':
                            if mirror_cell:
                                node_x  = 2*midline_coordinate_x -node_x
                                node_prev_x  = 2*midline_coordinate_x - node_prev_x
                        
                        xvals.append(node_x)
                        yvals.append(node_y)
                        zvals.append(node_z)
                        if show_reconstruction:
                            branch_x.append(node_x)
                            branch_y.append(node_y)
                            branch_z.append(node_z)
                        
                        nodeidx_dendrite_coarse = [np.argmax(node_x < xscale_dendrite_coarse),np.argmax(node_y < yscale_dendrite_coarse),np.argmax(node_z < zscale_dendrite_coarse)]
                        nodeidx_dendrite_coarse_prev = [np.argmax(node_prev_x < xscale_dendrite_coarse),np.argmax(node_prev_y < yscale_dendrite_coarse),np.argmax(node_prev_z < zscale_dendrite_coarse)]
                        if nodeidx_dendrite_coarse == nodeidx_dendrite_coarse_prev:
                            volume_dendrite_coarse[nodeidx_dendrite_coarse[0]][nodeidx_dendrite_coarse[1]][nodeidx_dendrite_coarse[2]] += node.distanceTo(node_prev)
                            
                        nodeidx_dendrite_fine = [np.argmax(node_x < xscale_dendrite_fine),np.argmax(node_y < yscale_dendrite_fine),np.argmax(node_z < zscale_dendrite_fine)]
                        nodeidx_dendrite_fine_prev = [np.argmax(node_prev_x < xscale_dendrite_fine),np.argmax(node_prev_y < yscale_dendrite_fine),np.argmax(node_prev_z < zscale_dendrite_fine)]
                        if nodeidx_dendrite_fine == nodeidx_dendrite_fine_prev:
                            volume_dendrite_fine[nodeidx_dendrite_fine[0]][nodeidx_dendrite_fine[1]][nodeidx_dendrite_fine[2]] += node.distanceTo(node_prev)
                        try:    
                            #node_id = node.annotation.id()  
                            node_id = node.getAnnotation().id()
                            node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                            #node_prev_id = node_prev.annotation.id()
                            node_prev_id = node_prev.getAnnotation().id()
                            node_prev_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_prev_id,'structureIdPath'].values[0].strip('/').split('/'),int)
                            for node_id in node_ids:
                                if node_id in node_prev_ids:
                                    allen_dendrite[(allen_df['structureId']==node_id).argmax()] += node.distanceTo(node_prev)
        
                            sum_dendrite_length += node.distanceTo(node_prev)
                        except:
                            pass
                            #print('no annotation for dendrite segment')
                    except:
                        break
                if show_reconstruction:
                    ax_dendrite_3d.plot(branch_x,branch_y,branch_z)
                    
            
            
            allen_axon_list_original.append(np.concatenate([allen_axon_ipsi,allen_axon_contra]))
            allen_dendrite_list_original.append(allen_dendrite)
            allen_axon_branch_point_list_original.append(np.concatenate([allen_axon_branch_point_ipsi,allen_axon_branch_point_contra]))
            allen_axon_end_point_list_original.append(np.concatenate([allen_axon_end_point_ipsi,allen_axon_end_point_contra]))
            allen_dendrite_branch_point_list_original.append(np.concatenate([allen_dendrite_branch_point_ipsi,allen_dendrite_branch_point_contra]))
            allen_dendrite_end_point_list_original.append(np.concatenate([allen_dendrite_end_point_ipsi,allen_dendrite_end_point_contra]))
            
                     
            
            allen_axon_normalized = np.concatenate([allen_axon_ipsi,allen_axon_contra]) /sum_axon_length
            allen_axon_list.append(allen_axon_normalized)
            allen_dendrite_normalized = allen_dendrite /sum_dendrite_length
            allen_dendrite_list.append(allen_dendrite_normalized)
            allen_soma_list.append(allen_soma)#no need to normalize
            allen_soma_distribution_list.append(allen_soma_distribution)#no need to normalize
            
            allen_axon_branch_point_normalized = np.concatenate([allen_axon_branch_point_ipsi,allen_axon_branch_point_contra]) /sum_axon_branch_point_num
            allen_axon_branch_point_list.append(allen_axon_branch_point_normalized)
            allen_axon_end_point_normalized = np.concatenate([allen_axon_end_point_ipsi,allen_axon_end_point_contra]) /sum_axon_end_point_num
            allen_axon_end_point_list.append(allen_axon_end_point_normalized)
            
            allen_dendrite_branch_point_normalized = np.concatenate([allen_dendrite_branch_point_ipsi,allen_dendrite_branch_point_contra]) /sum_dendrite_branch_point_num
            allen_dendrite_branch_point_list.append(allen_dendrite_branch_point_normalized)
            allen_dendrite_end_point_normalized = np.concatenate([allen_dendrite_end_point_ipsi,allen_dendrite_end_point_contra]) /sum_dendrite_end_point_num
            allen_dendrite_end_point_list.append(allen_dendrite_end_point_normalized)
            
            
            
            volume_coarse_filt = ndimage.filters.gaussian_filter(volume_coarse,[sigma_step,sigma_step,sigma_step])
            volume_fine_filt = ndimage.filters.gaussian_filter(volume_fine,[sigma_step,sigma_step,sigma_step])
            if normalizing_function_coarse == 'sum':
                volume_coarse_filt  = volume_coarse_filt /sum(volume_coarse_filt.flatten())
            elif normalizing_function_coarse == 'max':
                volume_coarse_filt  = volume_coarse_filt /max(volume_coarse_filt.flatten())
            if normalizing_function_fine == 'sum':
                volume_fine_filt  = volume_fine_filt /sum(volume_fine_filt.flatten())
            elif normalizing_function_fine == 'max':
                volume_fine_filt  = volume_fine_filt /max(volume_fine_filt.flatten())
            
            #volumes.append(np.concatenate([volume_coarse_filt.flatten(),volume_fine_filt.flatten()]))
            volumes_coarse.append(volume_coarse_filt.flatten())
            volumes_fine.append(volume_fine_filt.flatten())
            
            volume_dendrite_coarse_filt = ndimage.filters.gaussian_filter(volume_dendrite_coarse,[sigma_step,sigma_step,sigma_step])
            volume_dendrite_fine_filt = ndimage.filters.gaussian_filter(volume_dendrite_fine,[sigma_step,sigma_step,sigma_step])
            if normalizing_function_coarse == 'sum':
                volume_dendrite_coarse_filt  = volume_dendrite_coarse_filt /sum(volume_dendrite_coarse_filt.flatten())
            elif normalizing_function_coarse == 'max':
                volume_dendrite_coarse_filt  = volume_dendrite_coarse_filt /max(volume_dendrite_coarse_filt.flatten())
            if normalizing_function_fine == 'sum':
                volume_dendrite_fine_filt  = volume_dendrite_fine_filt /sum(volume_dendrite_fine_filt.flatten())
            elif normalizing_function_fine == 'max':
                volume_dendrite_fine_filt  = volume_dendrite_fine_filt /max(volume_dendrite_fine_filt.flatten())
            
            #volumes.append(np.concatenate([volume_dendrite_coarse_filt.flatten(),volume_dendrite_fine_filt.flatten()]))
            volumes_dendrite_coarse.append(volume_dendrite_coarse_filt.flatten())
            volumes_dendrite_fine.append(volume_dendrite_fine_filt.flatten())
            
            dendrite_lengths.append(dendritestats.getCableLength())
            axon_lengths.append(axonstats.getCableLength())
            axon_branch_numbers.append(axonstats.getNBranches())
            axon_branch_point_numbers.append(len(axon.getGraph().getBPs().toArray()))
            axon_end_point_numbers.append(sum_axon_end_point_num)
            axon_end_point_numbers_total.append(len(axonendings))
            dendrite_branch_numbers.append(dendritestats.getNBranches())
            dendrite_primary_branch_numbers.append(len(dendritestats.getPrimaryBranches()))
            dendrite_branch_point_numbers.append(len(dendrite.getGraph().getBPs().toArray()))
            dendrite_end_point_numbers.append(len(dendriteendings))
            
            #ccf_disagreements_list.append(ccf_disagreements_now)
            #%
            if show_reconstruction:
                ax_3d.set_title(neuron_ID)
                ax_3d.set_xlabel('X')
                ax_3d.set_ylabel('Y')
                ax_3d.set_zlabel('Z')
                ax_yz.imshow(np.max(volume_coarse_filt,0)) #
                ax_yz.set_ylabel('Y')
                ax_yz.set_xlabel('Z')
                ax_xz.imshow(np.max(volume_coarse_filt,1))
                ax_xz.set_ylabel('X')
                ax_xz.set_xlabel('Z')
                ax_xy.imshow(np.max(volume_coarse_filt,2))
                ax_xy.set_ylabel('X')
                ax_xy.set_xlabel('Y')
                
                
                ax_yz_fine.imshow(np.max(volume_fine_filt,0)) #
                ax_yz_fine.set_ylabel('Y')
                ax_yz_fine.set_xlabel('Z')
                ax_xz_fine.imshow(np.max(volume_fine_filt,1))
                ax_xz_fine.set_ylabel('X')
                ax_xz_fine.set_xlabel('Z')
                ax_xy_fine.imshow(np.max(volume_fine_filt,2))
                ax_xy_fine.set_ylabel('X')
                ax_xy_fine.set_xlabel('Y')
                
                # dendrite plot
                
                
                ax_dendrite_3d.set_title(neuron_ID)
                ax_dendrite_3d.set_xlabel('X')
                ax_dendrite_3d.set_ylabel('Y')
                ax_dendrite_3d.set_zlabel('Z')
                ax_dendrite_yz.imshow(np.max(volume_dendrite_coarse_filt,0)) #
                ax_dendrite_yz.set_ylabel('Y')
                ax_dendrite_yz.set_xlabel('Z')
                ax_dendrite_xz.imshow(np.max(volume_dendrite_coarse_filt,1))
                ax_dendrite_xz.set_ylabel('X')
                ax_dendrite_xz.set_xlabel('Z')
                ax_dendrite_xy.imshow(np.max(volume_dendrite_coarse_filt,2))
                ax_dendrite_xy.set_ylabel('X')
                ax_dendrite_xy.set_xlabel('Y')
                
                
                ax_dendrite_yz_fine.imshow(np.max(volume_dendrite_fine_filt,0)) #
                ax_dendrite_yz_fine.set_ylabel('Y')
                ax_dendrite_yz_fine.set_xlabel('Z')
                ax_dendrite_xz_fine.imshow(np.max(volume_dendrite_fine_filt,1))
                ax_dendrite_xz_fine.set_ylabel('X')
                ax_dendrite_xz_fine.set_xlabel('Z')
                ax_dendrite_xy_fine.imshow(np.max(volume_dendrite_fine_filt,2))
                ax_dendrite_xy_fine.set_ylabel('X')
                ax_dendrite_xy_fine.set_xlabel('Y')
                
                plt.show()
            #%
            X_edges = [np.nanmin(xvals),np.nanmax(xvals)]
            Y_edges = [np.nanmin(yvals),np.nanmax(yvals)]
            Z_edges = [np.nanmin(zvals),np.nanmax(zvals)]
        #break
            if show_reconstruction:
                break
        #%
        if show_reconstruction:
            break
    #%
    allen_axon_matrix = np.asarray(allen_axon_list)
    allen_dendrite_matrix = np.asarray(allen_dendrite_list)
    allen_soma_matrix = np.asarray(allen_soma_list)
    allen_soma_distribution_matrix = np.asarray(allen_soma_distribution_list)
    allen_axon_branch_points_matrix = np.asarray(allen_axon_branch_point_list)
    allen_axon_end_points_matrix = np.asarray(allen_axon_end_point_list)
    allen_dendrite_branch_points_matrix = np.asarray(allen_dendrite_branch_point_list)
    allen_dendrite_end_points_matrix = np.asarray(allen_dendrite_end_point_list)
    volumes_coarse_matrix = np.asarray(volumes_coarse)   
    volumes_fine_matrix = np.asarray(volumes_fine)  
    volumes_dendrite_coarse_matrix = np.asarray(volumes_dendrite_coarse)   
    volumes_dendrite_fine_matrix = np.asarray(volumes_dendrite_fine)  
    soma_coordinates_matrix = np.asarray([soma_coordinates['X'],soma_coordinates['Y'],soma_coordinates['Z']]).T 
    soma_coordinates_matrix_normalized = (soma_coordinates_matrix-np.min(soma_coordinates_matrix,0))/(np.max(soma_coordinates_matrix,0)-np.min(soma_coordinates_matrix,0))
    
    
    
    
    
    dendrite_lengths_matrix = np.asarray(dendrite_lengths)   
    edge_vals = np.percentile(dendrite_lengths_matrix,dendrite_length_normalize_percentiles)
    dendrite_lengths_matrix  = (dendrite_lengths_matrix-edge_vals[0])/np.diff(edge_vals)
    if dendrite_length_cut_outliers:
        dendrite_lengths_matrix[dendrite_lengths_matrix>1] = 1
        dendrite_lengths_matrix[dendrite_lengths_matrix<0] = 0
    dendrite_branch_point_numbers_matrix = np.asarray(dendrite_branch_point_numbers)   
    edge_vals = np.percentile(dendrite_branch_point_numbers_matrix,dendrite_branch_point_numbers_normalize_percentiles)
    dendrite_branch_point_numbers_matrix   = (dendrite_branch_point_numbers_matrix -edge_vals[0])/np.diff(edge_vals)
    if dendrite_length_cut_outliers:
        dendrite_branch_point_numbers_matrix[dendrite_branch_point_numbers_matrix>1] = 1
        dendrite_branch_point_numbers_matrix[dendrite_branch_point_numbers_matrix<0] = 0    
    dendrite_branch_numbers_matrix = np.asarray(dendrite_branch_numbers)   
    edge_vals = np.percentile(dendrite_branch_numbers_matrix,dendrite_branch_numbers_normalize_percentiles)
    dendrite_branch_numbers_matrix   = (dendrite_branch_numbers_matrix -edge_vals[0])/np.diff(edge_vals)
    if dendrite_length_cut_outliers:
        dendrite_branch_numbers_matrix[dendrite_branch_numbers_matrix>1] = 1
        dendrite_branch_numbers_matrix[dendrite_branch_numbers_matrix<0] = 0           
        
    dendrite_primary_branch_numbers_matrix = np.asarray(dendrite_primary_branch_numbers)   
    edge_vals = np.percentile(dendrite_primary_branch_numbers_matrix,dendrite_primary_branch_numbers_normalize_percentiles)
    dendrite_primary_branch_numbers_matrix   = (dendrite_primary_branch_numbers_matrix -edge_vals[0])/np.diff(edge_vals)        
        
    dendrite_end_point_numbers_matrix = np.asarray(dendrite_end_point_numbers)   
    edge_vals = np.percentile(dendrite_end_point_numbers_matrix,dendrite_end_point_numbers_normalize_percentiles)
    dendrite_end_point_numbers_matrix   = (dendrite_end_point_numbers_matrix -edge_vals[0])/np.diff(edge_vals)
    if dendrite_length_cut_outliers:
        dendrite_end_point_numbers_matrix[dendrite_end_point_numbers_matrix>1] = 1
        dendrite_end_point_numbers_matrix[dendrite_end_point_numbers_matrix<0] = 0
    
    
    
    axon_lengths_matrix = np.asarray(axon_lengths)   
    edge_vals = np.percentile(axon_lengths_matrix,axon_length_normalize_percentiles)
    axon_lengths_matrix  = (axon_lengths_matrix-edge_vals[0])/np.diff(edge_vals)
    if axon_length_cut_outliers:
        axon_lengths_matrix[axon_lengths_matrix>1] = 1
        axon_lengths_matrix[axon_lengths_matrix<0] = 0
    axon_branch_point_numbers_matrix = np.asarray(axon_branch_point_numbers)   
    edge_vals = np.percentile(axon_branch_point_numbers_matrix,axon_branch_point_numbers_normalize_percentiles)
    axon_branch_point_numbers_matrix   = (axon_branch_point_numbers_matrix -edge_vals[0])/np.diff(edge_vals)
    if axon_length_cut_outliers:
        axon_branch_point_numbers_matrix[axon_branch_point_numbers_matrix>1] = 1
        axon_branch_point_numbers_matrix[axon_branch_point_numbers_matrix<0] = 0    
    axon_branch_numbers_matrix = np.asarray(axon_branch_numbers)   
    edge_vals = np.percentile(axon_branch_numbers_matrix,axon_branch_numbers_normalize_percentiles)
    axon_branch_numbers_matrix   = (axon_branch_numbers_matrix -edge_vals[0])/np.diff(edge_vals)
    if axon_length_cut_outliers:
        axon_branch_numbers_matrix[axon_branch_numbers_matrix>1] = 1
        axon_branch_numbers_matrix[axon_branch_numbers_matrix<0] = 0    
    axon_end_point_numbers_matrix = np.asarray(axon_end_point_numbers)   
    edge_vals = np.percentile(axon_end_point_numbers_matrix,axon_end_point_numbers_normalize_percentiles)
    axon_end_point_numbers_matrix   = (axon_end_point_numbers_matrix -edge_vals[0])/np.diff(edge_vals)
    if axon_length_cut_outliers:
        axon_end_point_numbers_matrix[axon_end_point_numbers_matrix>1] = 1
        axon_end_point_numbers_matrix[axon_end_point_numbers_matrix<0] = 0
    
    
    soma_locations_header = np.unique(soma_locations)
    
    uniquesomashapes = np.unique(subtle_features['Soma shape'])
    uniqueboutontypes = np.unique(np.concatenate(subtle_features['Bouton types']))
    uniqueaxonoriginpoints = np.unique(subtle_features['Axon origin point'])
    subtle_features_list = list()
    soma_locations_list = list()
    for neuron_idx,neuron_ID in enumerate(cell_names):
        subtle_features_header = list()
        subtle_features_neuron = list()
        
        
        
        for somashape in uniquesomashapes:
            subtle_features_header.append('Soma shape - {}'.format(somashape))
            if subtle_features['Soma shape'][neuron_idx] == somashape:
                subtle_features_neuron.append(True)
            else:
                subtle_features_neuron.append(False)
        for axonorigin in uniqueaxonoriginpoints:
            subtle_features_header.append('Axon origin point - {}'.format(axonorigin))
            if subtle_features['Axon origin point'][neuron_idx] == axonorigin:
                subtle_features_neuron.append(True)
            else:
                subtle_features_neuron.append(False)
                
        for boutontype in uniqueboutontypes:
            subtle_features_header.append('Bouton type - {}'.format(boutontype))
            if boutontype in subtle_features['Bouton types'][neuron_idx]:
                subtle_features_neuron.append(True)
            else:
                subtle_features_neuron.append(False)        
        
        subtle_features_header.append('Dendritic spine')
        try:
            subtle_features_neuron.append(subtle_features['Dendritic spines'][neuron_idx])
        except:
            print('error with dendritic spine data in tabel for neuron {}'.format(neuron_ID))
            subtle_features_neuron.append('-')
        subtle_features_list.append(subtle_features_neuron)
        
        
        soma_locations_neuron = list()
        for somalocation in soma_locations_header:
            if soma_locations[neuron_idx] == somalocation:
                soma_locations_neuron.append(True)
            else:
                soma_locations_neuron.append(False)
        
        soma_locations_list.append(soma_locations_neuron)
        
    subtle_features_matrix = np.asarray(subtle_features_list)    
    soma_locations_matrix = np.asarray(soma_locations_list)  
    
    #% manually annotated target locations and targe
    target_areas_all = np.asarray([])
    target_locations_all = np.asarray([])
    for target_area, target_location in zip(target_areas,target_locations):
        target_locations_all = np.concatenate([target_locations_all,target_location])
        target_areas_all = np.concatenate([target_areas_all,target_area])
    try:
        target_areas_header = np.unique(target_areas_all )
    except:
        print('target areas not filled out properly')
        target_areas_header = [target_areas_all[0]]
    target_locations_header = np.unique(target_locations_all )
    #% construct header for manually selected target locations
    target_location_potential_parents =  ['Cerebellar cortex','Cerebellar nuclei','Medulla','Pons','Midbrain','Thalamus','Hypothalamus','fiber tracts', 'ventricular systems']
    target_locations_header_with_parent = list()
    todel = list()
    parents_list = list()
    for i,target_locations_head in enumerate(target_locations_header):
        try:
            node_prev_ids = np.asarray(allen_df.loc[allen_df['name']==target_locations_head,'structureIdPath'].values[0].strip('/').split('/'),int)
            for node_prev_id in node_prev_ids:
                parentname = allen_df.loc[allen_df['structureId']==int(node_prev_id),'name'].values[0]
                if parentname in target_location_potential_parents:
                    break
            target_locations_header_with_parent.append('{} - {}'.format(parentname,target_locations_head))
            parents_list.append(parentname)
        except:
            todel.append(i)
            print('no parents found for {}, please check, deleting this nucleus..'.format(target_locations_head))
    
    neededhead = np.asarray(np.ones(len(target_locations_header)),bool)
    neededhead[todel]=0
    target_locations_header=target_locations_header[neededhead]
    target_locations_header_with_parent = np.asarray(target_locations_header_with_parent)
    order = np.argsort(target_locations_header_with_parent)
    target_locations_header = target_locations_header[order]
    target_locations_header_with_parent = target_locations_header_with_parent[order]
    parents_list = np.asarray(parents_list)[order]
    #% reorganizing headers to follow the order of the parents in the user defined list
    order = list()
    for parent_now in target_location_potential_parents:
        order.append(np.where(parents_list == parent_now)[0])
    order = np.concatenate(order)
    target_locations_header = target_locations_header[order]
    target_locations_header_with_parent = target_locations_header_with_parent[order]
    #%
    target_locations_list = list()
    target_areas_list = list()
    for target_areas_now, target_locations_now in zip(target_areas,target_locations):
        target_areas_array = np.zeros(len(target_areas_header))
        for target_area in target_areas_now:
            target_areas_array[target_areas_header==target_area]=1
        target_areas_list.append(target_areas_array)
        target_locations_array = np.zeros(len(target_locations_header))
        for target_location in target_locations_now:
            target_locations_array[target_locations_header==target_location]=1
        target_locations_list.append(target_locations_array)
    #%
    original_data = dict()
    if not parameters_dict['project'] == 'medulla':
        soma_area_in_the_medulla = list()
        
    original_data={'volumes_coarse':volumes_coarse.copy(),
                  'volumes_fine':volumes_fine.copy(),
                  'volumes_dendrite_coarse':volumes_dendrite_coarse.copy(),
                  'volumes_dendrite_fine':volumes_dendrite_fine.copy(),
                  'allen_axon_list':allen_axon_list.copy(),
                  'allen_dendrite_list':allen_dendrite_list.copy(),
                  'allen_soma_list':allen_soma_list.copy(),
                  'allen_soma_distribution_list':allen_soma_distribution_list.copy(),
                  'allen_axon_branch_point_list':allen_axon_branch_point_list.copy(),
                  'allen_axon_end_point_list':allen_axon_end_point_list.copy(),
                  'allen_dendrite_branch_point_list' :allen_dendrite_branch_point_list.copy(),
                  'allen_dendrite_end_point_list' :allen_dendrite_end_point_list.copy(),
                  'allen_axon_list_original' :allen_axon_list_original.copy(),
                  'allen_dendrite_list_original' :allen_dendrite_list_original.copy(),
                  'allen_soma_list_original' :allen_soma_list_original.copy(),
                  'allen_soma_distribution_list_original' :allen_soma_distribution_list_original.copy(),
                  'allen_axon_branch_point_list_original' :allen_axon_branch_point_list_original.copy(),
                  'allen_axon_end_point_list_original' :allen_axon_end_point_list_original.copy(),
                  'allen_dendrite_branch_point_list_original' :allen_dendrite_branch_point_list_original.copy(),
                  'allen_dendrite_end_point_list_original':allen_dendrite_end_point_list_original.copy(),
                  'cell_names' :cell_names.copy(),
                  'cell_juci_clusters' :cell_juci_clusters.copy(),
                  'soma_coordinates' :soma_coordinates.copy(),
                  'soma_locations' :soma_locations.copy(),
                  'soma_area_in_the_medulla' :soma_area_in_the_medulla.copy(),
                  'subtle_features' :subtle_features.copy(),
                  'dendrite_lengths' :dendrite_lengths.copy(),
                  'axon_lengths' :axon_lengths.copy(),
                  'axon_branch_point_numbers' :axon_branch_point_numbers.copy(),
                  'axon_end_point_numbers' :axon_end_point_numbers.copy(),
                  'axon_end_point_numbers_total' :axon_end_point_numbers_total.copy(),
                  'dendrite_branch_point_numbers' :dendrite_branch_point_numbers.copy(),
                  'dendrite_branch_numbers' :dendrite_branch_numbers.copy(),
                  'dendrite_primary_branch_numbers' :dendrite_primary_branch_numbers.copy(),
                  'dendrite_end_point_numbers' :dendrite_end_point_numbers.copy(),
                  'dendrite_bounding_sphere_radius':dendrite_bounding_ball_radius,
                  'dendrite_bounding_sphere_coordinates':dendrite_bounding_ball_coordinates,
                  'allen_axon_matrix' :allen_axon_matrix.copy(),
                  'allen_dendrite_matrix' :allen_dendrite_matrix.copy(),
                  'allen_soma_matrix' :allen_soma_matrix.copy(),
                  'allen_soma_distribution_matrix' :allen_soma_distribution_matrix.copy(),
                  'allen_axon_branch_points_matrix' :allen_axon_branch_points_matrix.copy(),
                  'allen_axon_end_points_matrix' :allen_axon_end_points_matrix.copy(),
                  'allen_dendrite_branch_points_matrix' :allen_dendrite_branch_points_matrix.copy(),
                  'allen_dendrite_end_points_matrix' :allen_dendrite_end_points_matrix.copy(),
                  'volumes_coarse_matrix' :volumes_coarse_matrix.copy(),
                  'volumes_fine_matrix' :volumes_fine_matrix.copy(),
                  'volumes_dendrite_coarse_matrix' :volumes_dendrite_coarse_matrix.copy(),
                  'volumes_dendrite_fine_matrix' :volumes_dendrite_fine_matrix.copy(),
                  'soma_coordinates_matrix' :soma_coordinates_matrix.copy(),
                  'soma_coordinates_matrix_normalized' :soma_coordinates_matrix_normalized.copy(),
                  'dendrite_lengths_matrix' :dendrite_lengths_matrix.copy(),
                  'dendrite_branch_point_numbers_matrix' :dendrite_branch_point_numbers_matrix.copy(),
                  'dendrite_branch_numbers_matrix' :dendrite_branch_numbers_matrix.copy(),
                  'dendrite_primary_branch_numbers_matrix' :dendrite_primary_branch_numbers_matrix.copy(),
                  'dendrite_end_point_numbers_matrix' :dendrite_end_point_numbers_matrix.copy(),
                  'axon_lengths_matrix' :axon_lengths_matrix.copy(),
                  'axon_branch_point_numbers_matrix' :axon_branch_point_numbers_matrix.copy(),
                  'axon_branch_numbers_matrix' :axon_branch_numbers_matrix.copy(),
                  'axon_end_point_numbers_matrix' :axon_end_point_numbers_matrix.copy(),
                  'subtle_features_list' :subtle_features_list.copy(),
                  'subtle_features_header':subtle_features_header,
                  'soma_locations_list' :soma_locations_list.copy(),
                  'subtle_features_matrix' :subtle_features_matrix.copy(),
                  'soma_locations_matrix' :soma_locations_matrix.copy(),
                  'target_locations_matrix': np.asarray(target_locations_list),
                  'target_areas_matrix': np.asarray(target_areas_list),
                  'target_areas_header':target_areas_header,
                  'target_locations_header':target_locations_header,
                  'target_locations_header_with_parent':target_locations_header_with_parent,
                  'parameters_dict':parameters_dict,
                  'allen_df':allen_df}
    
    
    
    
    if parameters_dict['use_probabilistic_axon_ending_map']:
        np.save(os.path.join(parameters_dict['save_dir'],'{}_exported_{}_probabilistic_axon_endings.npy'.format(parameters_dict['project'],str(datetime.date.today()))),original_data)
    else:
        np.save(os.path.join(parameters_dict['save_dir'],'{}_exported_{}.npy'.format(parameters_dict['project'],str(datetime.date.today()))),original_data)

def generate_big_matrix(data,allen_df,parameters_dict):
    weight_dict = parameters_dict['weights_dict']
    data['allen_axon_matrix_now'] = data['allen_axon_matrix'].copy()*data['allen_needed_nuclei']
    data['allen_axon_branch_points_matrix_now'] = data['allen_axon_branch_points_matrix'].copy()*data['allen_needed_nuclei']
    
    data['allen_axon_end_points_matrix_now'] = data['allen_axon_end_points_matrix'].copy()
    if parameters_dict['merge_ipsi_contra_projections']:
        data['allen_axon_end_points_matrix_now'][:,:len(allen_df)]=data['allen_axon_end_points_matrix_now'][:,:len(allen_df)] + data['allen_axon_end_points_matrix_now'][:,len(allen_df):]
        data['allen_axon_end_points_matrix_now'][:,len(allen_df):] = 0
    if parameters_dict['normalize_endpoints_to_basic_cell_groups']:
        basic_cell_groups_idx = np.where((allen_df['structureId']==8).values)[0][0]
        values_ipsi = data['allen_axon_end_points_matrix_now'][:,basic_cell_groups_idx]
        values_contra = data['allen_axon_end_points_matrix_now'][:,basic_cell_groups_idx+len(allen_df)]
        values_contra = data['allen_axon_end_points_matrix_now'] = data['allen_axon_end_points_matrix_now']*(values_contra+values_ipsi)[:,np.newaxis]
    data['allen_axon_end_points_matrix_now'] = data['allen_axon_end_points_matrix_now']*data['allen_needed_nuclei']
    data['allen_dendrite_matrix_now'] = data['allen_dendrite_matrix'].copy()*data['allen_needed_nuclei'][:len(allen_df)]
    data['allen_dendrite_branch_points_matrix_now'] = data['allen_dendrite_branch_points_matrix'].copy()*data['allen_needed_nuclei']
    data['allen_dendrite_end_points_matrix_now'] = data['allen_dendrite_end_points_matrix'].copy()*data['allen_needed_nuclei']
    data['allen_soma_matrix_now'] = data['allen_soma_matrix'].copy()*data['allen_needed_nuclei'][:len(allen_df)]
    data['allen_soma_distribution_matrix_now'] = data['allen_soma_distribution_matrix'].copy()*data['allen_needed_nuclei'][:len(allen_df)]

    for key_now in ['allen_axon_matrix_now',
                    'allen_axon_branch_points_matrix_now',
                    'allen_axon_end_points_matrix_now',
                    'allen_dendrite_matrix_now',
                    'allen_dendrite_branch_points_matrix_now',
                    'allen_dendrite_end_points_matrix_now',
                    'allen_soma_matrix_now',
                    'allen_soma_distribution_matrix_now']:
        if np.any(np.isnan(data[key_now])):
            print('nan in {} - zeroing'.format(key_now))
            data[key_now][np.isnan(data[key_now])] = 0

    if parameters_dict['normalize_endpoints_to_volume']:
        volumes = np.concatenate([allen_df['volume'].values]*2)
        data['allen_axon_end_points_matrix_now'] = data['allen_axon_end_points_matrix_now']*data['axon_end_point_numbers'][:,None]
        data['allen_axon_end_points_matrix_now'] = data['allen_axon_end_points_matrix_now'] /volumes
    elif parameters_dict['use_actual_end_point_numbers']:
        data['allen_axon_end_points_matrix_now'] = data['allen_axon_end_points_matrix_now']*data['axon_end_point_numbers'][:,None]
        
        
    # calculate dendrite metrics
    axon_internodes = np.asarray(data['axon_lengths'])/(np.asarray(data['axon_end_point_numbers'])+np.asarray(data['axon_branch_point_numbers']))
    dendrite_internodes = np.asarray(data['dendrite_lengths'])/(np.asarray(data['dendrite_end_point_numbers'])+np.asarray(data['dendrite_branch_point_numbers']))

    data['big_matrix'] = np.concatenate([data['subtle_features_matrix']*weight_dict['subtle_features_weight'],
                                         data['soma_coordinates_matrix_normalized']*weight_dict['soma_coordinates_weight'],
                                         data['target_locations_matrix']*weight_dict['target_locations_weight'],
                                         data['target_areas_matrix']*weight_dict['target_areas_weight'],
                                         data['volumes_coarse_matrix']*weight_dict['axon_projection_coarse_weight'],
                                         data['volumes_fine_matrix']*weight_dict['axon_projection_fine_weight'],
                                         data['volumes_dendrite_coarse_matrix']*weight_dict['dendrite_projection_coarse_weight'],
                                         data['volumes_dendrite_fine_matrix']*weight_dict['dendrite_projection_fine_weight'],
                                         data['axon_lengths_matrix'][:,np.newaxis]*weight_dict['axon_length_weight'],
                                         data['axon_branch_point_numbers_matrix'][:,np.newaxis]*weight_dict['axon_branch_point_number_weight'],
                                         data['axon_branch_numbers_matrix'][:,np.newaxis]*weight_dict['axon_branch_number_weight'],
                                         data['axon_end_point_numbers_matrix'][:,np.newaxis]*weight_dict['axon_end_point_number_weight'],
                                         data['dendrite_lengths_matrix'][:,np.newaxis]*weight_dict['dendrite_length_weight'],
                                         axon_internodes[:,np.newaxis]*weight_dict['axon_internode_distance_weight'],
                                         dendrite_internodes[:,np.newaxis]*weight_dict['dendirte_internode_distance_weight'],
                                         data['dendrite_bounding_sphere_radius'][:,np.newaxis]*weight_dict['dendrite_bounding_sphere_weight'],
                                         data['dendrite_branch_point_numbers_matrix'][:,np.newaxis]*weight_dict['dendrite_branch_point_number_weight'],
                                         data['dendrite_branch_numbers_matrix'][:,np.newaxis]*weight_dict['dendrite_branch_number_weight'],
                                         data['dendrite_end_point_numbers_matrix'][:,np.newaxis]*weight_dict['dendrite_end_point_number_weight'],
                                         data['dendrite_primary_branch_numbers_matrix'][:,np.newaxis]*weight_dict['dendrite_primary_branch_number_weight'],
                                         data['soma_locations_matrix']*weight_dict['soma_locations_weight'],
                                         data['allen_axon_matrix_now']*weight_dict['allen_axon_weight'],
                                         data['allen_axon_branch_points_matrix_now']*weight_dict['allen_axon_branch_points_weight'],
                                         data['allen_axon_end_points_matrix_now']*weight_dict['allen_axon_end_points_weight'],
                                         data['allen_dendrite_matrix_now']*weight_dict['allen_dendrite_weight'],                             
                                         data['allen_dendrite_branch_points_matrix_now']*weight_dict['allen_dendrite_branch_points_weight'],
                                         data['allen_dendrite_end_points_matrix_now']*weight_dict['allen_dendrite_end_points_weight'],
                                         data['allen_soma_matrix_now']*weight_dict['allen_soma_weight'],
                                         data['allen_soma_distribution_matrix_now']*weight_dict['allen_soma_distribution_weight']],1)
    data['big_matrix'] = data['big_matrix'][:,np.sum(data['big_matrix'],0)!=0]
    if parameters_dict['zscore_big_matrix']:
        
        data['big_matrix']=scipy.stats.zscore(data['big_matrix'],0)
    return data