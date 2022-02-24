#%% 1 - Import packages & set parameters
import os
os.chdir(r'C:\Users\judith.baka\OneDrive - Allen Institute\Desktop\Scripts\Mouselight_scripts')
import NeuronMorphologyAnalysis as nma
import numpy as np
import matplotlib.pyplot as plt
%matplotlib qt
parameters_dict = {'ccf_version':'3.0',#'3.0' or '2.5'
                   'project': 'medulla',#'thalamus'#'cortex'#'medulla' #  
                   'load_volumes_with_allen_df' : True,                 
                   'allocate_soma_in_white_matter_to_closest_grey' : True,
                   'allocate_axon_ending_in_white_matter_to_closest_grey': True,
                   'use_probabilistic_axon_ending_map' : False, # this smears the axon endings
                   'probabilistic_axon_ending_map_sigma' : 40,
                   'correction_steps' : list(range(5,95,5))#[5,10,15,20,30,40,50,60,70,80]}
                   }

#minimum depth is the CCF depth you minimally have to be to assign an axon ending/soma location to a brain region (e.g. Thalamus is not specific enough)
minimum_depth_dict = {'Hindbrain':6, 
                      'Interbrain':6,
                      'Midbrain':5,
                      'Cerebellar cortex':5,
                      'Cerebellar nuclei':4,
                      'Cerebral cortex': 7,
                      'Cerebral nuclei':6}
parameters_dict['minimum_depth_dict'] = minimum_depth_dict
parameters_dict = nma.utils.parameters.set_parameters(parameters_dict) # - 
del minimum_depth_dict
#============================================================================================================================
#============================================================================================================================
#%% 2 - Update metadata (optional)
# Set parameters on top and run.
nma.utils.online_notebook.update_metadata(parameters_dict['online_notebook_ID'], parameters_dict['metadata_dir'])  
#============================================================================================================================
#============================================================================================================================
#%% 3 - Run the actual analysis (optional) #kiszedi a json fileokbol az adatokat
allen_df =  nma.ccf.load_allen_dataframe(parameters_dict)
allen_df = allen_df[allen_df['volume']>0]
allen_df = nma.ccf.define_final_nuclei_for_endpoint_search(allen_df,parameters_dict)
nma.analysis.analyze_json_files(allen_df,parameters_dict)
#============================================================================================================================
#============================================================================================================================
#%% 4 - Load data #mindig meg kell csinalni! megjelenik egy ablak amiben kivalasztod a talan legujabb npy filet
original_data, parameters_dict, allen_df = nma.io.load_raw_data_matrices(parameters_dict)
#============================================================================================================================
#============================================================================================================================
#%% 5 - select cell group to analyze
# =============================================================================
# choose from this list:
# 'all'
# 'light red'
# 'dark red'
# 'dark green'
# 'light green'
# 'light blue'
# 'dark blue'
# 'gray'
# 'motor'
# 'soma location'   - also set 'soma_locations_needed'
# 'axon projection' - also set 'axon_projections_needed'
# 'Cerebellum':
# 'Thalamus_Hypothalamus':
# 'Midbrain':
# 'Medulla_Pons':
# =============================================================================
subselect_parameters = {'needed_group' : 'all',
                        'soma_locations_needed' :None, #['Intermediate reticular nucleus', 'Medullary reticular nucleus'] or None
                        'axon_projections_needed':None #['VII','XII','AMB','V'] or None
                        }
data = nma.analysis.subselect_cells(allen_df,
                                    original_data,
                                    subselect_parameters['needed_group'],
                                    soma_locations_needed=subselect_parameters['soma_locations_needed'], 
                                    axon_projections_needed=subselect_parameters['axon_projections_needed']) 
parameters_dict['cell_subselection'] = subselect_parameters
del subselect_parameters
#============================================================================================================================
#============================================================================================================================
#%% 6- Set weights - compose big properties matrix
parameters_dict['normalize_endpoints_to_volume']= False ############### EZZEL CSAK JATSZADOZUNK!

max_depth_dict = {'Hindbrain':6, #%  defining maximum depth to be considered in analysis for EVERYTHING
                    'Interbrain':6,
                    'Midbrain':5,
                    'Cerebellar cortex':5,
                    'Cerebellar nuclei':4,
                    'Cerebral cortex': 7,
                    'Cerebral nuclei':6}
weight_dict = {'allen_axon_weight':0,
              'allen_axon_branch_points_weight':0,
              'allen_axon_end_points_weight':2,
              'allen_dendrite_branch_points_weight':0,
              'allen_dendrite_end_points_weight':0,
              'allen_dendrite_weight':0,
              'allen_soma_weight':0,
              'allen_soma_distribution_weight':0,
              'soma_coordinates_weight':0,
              'soma_locations_weight':0,
              'axon_projection_coarse_weight':0,
              'axon_projection_fine_weight':0,
              'axon_length_weight':0,
              'axon_branch_point_number_weight':0,
              'axon_end_point_number_weight':1,
              'dendrite_branch_point_number_weight':0,
              'dendrite_end_point_number_weight':0,
              'dendrite_projection_coarse_weight':0,
              'dendrite_projection_fine_weight':0,
              'dendrite_length_weight':0,
              'subtle_features_weight':0,#1.5#1#0.5#1,
              'target_areas_weight':0,
              'target_locations_weight':0}
parameters_dict['maximum_depth_dict'] = max_depth_dict
parameters_dict['weights_dict'] = weight_dict

data['allen_needed_nuclei'] = nma.ccf.generate_list_of_nuclei_by_depth(allen_df,parameters_dict)
data = nma.analysis.generate_big_matrix(data,allen_df,parameters_dict)

del max_depth_dict,weight_dict
#============================================================================================================================
#============================================================================================================================ 
#%% 7 - Set plot parameters
plot_parameters = {'juci_cluster_names' : {'Cerebellum':'Cerebellum',
                                           'Medulla_Pons':'Medulla/Pons',
                                           'Midbrain':'Midbrain', 
                                           'Thalamus_Hypothalamus':'Thalamus/Hypothalamus'},
                   'juci_cluster_colors' : {'Cerebellum':'red', 
                                            'Medulla_Pons':'dodgerblue',
                                            'Midbrain':'limegreen',  #yellowgreen
                                            'Thalamus_Hypothalamus':'gray'}} #whitesmoke

#============================================================================================================================
#============================================================================================================================ 
#%% 8 - PCA (optional)
from sklearn.decomposition import PCA 
    
#volumes_matrix_zscore = scipy.stats.zscore(volumes_matrix[:,np.sum(volumes_matrix,0)>0])
data['big_matrix_pca'] = PCA(n_components=np.min([50,len(data['cell_names'])]))    
data['big_matrix_pcs'] = data['big_matrix_pca'].fit_transform(data['big_matrix'])##volumes_matrix
nma.utils.plot.plot_pca(data, parameters_dict, plot_parameters)
#============================================================================================================================
#============================================================================================================================ 
#%% 9 - tSNE & uMAP (optional)
from sklearn.manifold import TSNE  
import umap

# ----------------------------- PARAMETERS -----------------------------
parameters_dict['pca_parameters'] = {'perplexity': 3,#4
                                     'learning_rate' : 20}#14#10#5#14
parameters_dict['umap_parameters'] = {'learning_rate': .3,#4
                                     'n_neigbors' : 14,#14#10#5#14
                                     'min_dist' : .01,
                                     'n_epochs':500}
# ----------------------------- PARAMETERS -----------------------------

tsne = TSNE(n_components=2, verbose=1, perplexity=parameters_dict['pca_parameters']['perplexity'], n_iter=1000,learning_rate=parameters_dict['pca_parameters']['learning_rate'])
data['big_matrix_tsne'] = tsne.fit_transform(data['big_matrix'])#pcs

reducer = umap.UMAP()
reducer.verbose=True
reducer.learning_rate = parameters_dict['umap_parameters']['learning_rate']
reducer.n_neigbors = parameters_dict['umap_parameters']['n_neigbors']
reducer.min_dist = parameters_dict['umap_parameters']['min_dist']
reducer.n_epochs=parameters_dict['umap_parameters']['n_epochs']
data['big_matrix_umap'] = reducer.fit_transform(data['big_matrix'])
#%
nma.utils.plot.plot_tsne_umap(data, parameters_dict, plot_parameters)
#============================================================================================================================
#============================================================================================================================ 
#%% 10 - clustering (optional)
parameters_dict['clustering_method'] = 'weighted'#'ward')#'weighted')# 'average')#'complete')#pcs
data = nma.utils.plot.plot_clustering(data,parameters_dict,plot_parameters)
#============================================================================================================================
#============================================================================================================================ 
#%% 11 - projection pattern
#data['allen_axon_length_matrix'] = data['allen_axon_matrix']*data['axon_lengths'][:,None]/1000
plot_parameters['projection_pattern']={'separate_ipsi_contra_projections':False,
                                        'normalize_to_volume':False,
                                        'add_subtle_features':False,
                                        'base_nuclei' : ['Cerebellar cortex','Cerebellar nuclei','Medulla','Pons','Midbrain','Thalamus','Hypothalamus','fiber tracts', 'ventricular systems'],
                                        'nuclei_ordering':'abc',#'rastermap','clustering'
                                        'projection_matrix_to_use': 'allen_axon_end_points_matrix', #'allen_axon_length_matrix',#*axon_end_point_numbers[:,None] #allen_axon_matrix*axon_lengths[:,None]/1000 #allen_axon_matrix_now*axon_lengths[:,None]/1000 #allen_axon_end_point_list_original#allen_axon_end_points_matrix#allen_axon_end_points_matrix*axon_end_point_numbers[:,None] #np.asarray(allen_axon_branch_point_list_original)#allen_axon_end_points_matrix'
                                        'use_log_scale':False,
                                        'vmin':0.000001,
                                        'gridcolor' : 'black',
                                        'gridwidth' : .5,
                                        }
nma.utils.plot.plot_projection_patterns(data, allen_df, plot_parameters)
#============================================================================================================================
#============================================================================================================================ 
#%% 12 - ground truth projection pattern (curated by Juci)
nma.utils.plot.plot_ground_truth_output_nuclei(data,allen_df)
#============================================================================================================================
#============================================================================================================================ 
#%% 13 - initialize viewer
snt_viewer_3d, snt_viewer_3d_cells_added = nma.utils.snt.init_viewer()
#%% 14 - add cells to the viewer


#%% add a list of cells to the viewer
# =============================================================================
# choose from this list:
# 'all'
# 'light red'
# 'dark red'
# 'dark green'
# 'light green'
# 'light blue'
# 'dark blue'
# 'gray'
# 'yellow'  gold
# 'motor'
# 'soma location'   - also set 'soma_locations_needed'
# 'axon projection' - also set 'axon_projections_needed'
# 'Cerebellum': red
# 'Thalamus_Hypothalamus': gray
# 'Midbrain': limegreen
# 'Medulla_Pons': blue
# =============================================================================

subselect_parameters = {'needed_group' : 'yellow',
                        'soma_locations_needed' :None, #['Intermediate reticular nucleus', 'Medullary reticular nucleus'] or None
                        'axon_projections_needed':None #['VII','XII','AMB','V'] or None
                        }
needed_cells = nma.analysis.generate_cell_list(allen_df,
                                               original_data,
                                               subselect_parameters['needed_group'],
                                               soma_locations_needed=subselect_parameters['soma_locations_needed'], 
                                               axon_projections_needed=subselect_parameters['axon_projections_needed']) 
axon_color = 'gray'
dendrite_color = 'gray'



cells_added = {}
for cell_to_add in needed_cells:
    nma.utils.snt.add_cell_to_viewer(snt_viewer_3d,
                                     jsonfolder = parameters_dict['json_dir'],
                                     cell_name = cell_to_add,
                                     cells_added = {},
                                     axon_color= axon_color,
                                     dendrite_color = dendrite_color,
                                     show_axon = True,
                                     show_dendrite = True,
                                     all_cells_on_left = True,)