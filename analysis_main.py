#%% 1 - Import packages & set parameters
import os
os.chdir(r'C:\Users\judith.baka\OneDrive - Allen Institute\Desktop\Scripts\Mouselight_scripts')
#os.chdir(r'C:\Users\Judith\OneDrive - Allen Institute\Desktop\Scripts\Mouselight_scripts')
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
                   'correction_steps' : list(range(5,95,5)),#[5,10,15,20,30,40,50,60,70,80]}
                   'assign_spinal_cord_to_ccf':True
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
parameters_dict = nma.utils.parameters.set_parameters(parameters_dict) # - correct paths

#============================================================================================================================
#============================================================================================================================
#%% 5 - select cell group to analyze - mindig meg kell csinalni
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
#%% 6- Set weights - compose big properties matrix - mindig me kell csinalni
#---------PROJECTIONS---------------------------
parameters_dict['normalize_endpoints_to_volume']= False ############### EZZEL CSAK JATSZADOZUNK!
parameters_dict['use_actual_end_point_numbers'] = False
parameters_dict['merge_ipsi_contra_projections'] = True
parameters_dict['normalize_endpoints_to_basic_cell_groups'] = True
parameters_dict['minimum_depth']=3
max_depth_dict = {'Hindbrain':6, #%  defining maximum depth to be considered in analysis for EVERYTHING
                    'Interbrain':6,
                    'Midbrain':5,
                    'Cerebellar cortex':5,
                    'Cerebellar nuclei':4,
                    'Cerebral cortex': 7,
                    'Cerebral nuclei':6,
                    'Spinal cord':4}
max_depth_dict = {'Hindbrain':3, #%  defining maximum depth to be considered in analysis for EVERYTHING
                    'Interbrain':3,
                    'Midbrain':3,
                    'Cerebellar cortex':3,
                    'Cerebellar nuclei':3,
                    'Cerebral cortex': 3,
                    'Cerebral nuclei':3,
                    'Spinal cord':4}
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
              'axon_branch_number_weight':0,
              'axon_end_point_number_weight':0,#1,
              'axon_internode_distance_weight':0,
              'dendrite_branch_point_number_weight':0,
              'dendrite_branch_number_weight':0,
              'dendrite_primary_branch_number_weight':0,
              'dendrite_bounding_sphere_weight':0,
              'dendrite_end_point_number_weight':0,
              'dendrite_projection_coarse_weight':0,
              'dendrite_projection_fine_weight':0,
              'dendrite_length_weight':0,
              'dendirte_internode_distance_weight':0,
              'subtle_features_weight':0,#1.5#1#0.5#1,
              'target_areas_weight':0,
              'target_locations_weight':0}
parameters_dict['maximum_depth_dict'] = max_depth_dict
parameters_dict['weights_dict'] = weight_dict
parameters_dict['zscore_big_matrix'] = False  # NOT for projection analysis

data['allen_needed_nuclei'] = nma.ccf.generate_list_of_nuclei_by_depth(allen_df,parameters_dict)
data = nma.analysis.generate_big_matrix(data,allen_df,parameters_dict)

del max_depth_dict,weight_dict
#============================================================================================================================
#============================================================================================================================ 

#%% 6- Set weights - compose big properties matrix - mindig me kell csinalni
#---------Dendrites---------------------------
parameters_dict['normalize_endpoints_to_volume']= False ############### EZZEL CSAK JATSZADOZUNK!
parameters_dict['use_actual_end_point_numbers'] = False
parameters_dict['merge_ipsi_contra_projections'] = True
parameters_dict['normalize_endpoints_to_basic_cell_groups'] = True
parameters_dict['minimum_depth']=3
max_depth_dict = {'Hindbrain':6, #%  defining maximum depth to be considered in analysis for EVERYTHING
                    'Interbrain':6,
                    'Midbrain':5,
                    'Cerebellar cortex':5,
                    'Cerebellar nuclei':4,
                    'Cerebral cortex': 7,
                    'Cerebral nuclei':6,
                    'Spinal cord':4}
weight_dict = {'allen_axon_weight':0,
              'allen_axon_branch_points_weight':0,
              'allen_axon_end_points_weight':0,
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
              'axon_branch_number_weight':0,
              'dendrite_primary_branch_number_weight':0,
              'axon_end_point_number_weight':0,#1,
              'axon_internode_distance_weight':0,
              'dendrite_branch_point_number_weight':1,
              'dendrite_branch_number_weight':1,
              'dendrite_bounding_sphere_weight':1,
              'dendrite_end_point_number_weight':1,
              'dendrite_projection_coarse_weight':0,
              'dendrite_projection_fine_weight':0,
              'dendrite_length_weight':1,
              'dendirte_internode_distance_weight':1,
              'subtle_features_weight':0,#1.5#1#0.5#1,
              'target_areas_weight':0,
              'target_locations_weight':0}
parameters_dict['maximum_depth_dict'] = max_depth_dict
parameters_dict['weights_dict'] = weight_dict
parameters_dict['zscore_big_matrix'] = True  # NOT for projection analysis

data['allen_needed_nuclei'] = nma.ccf.generate_list_of_nuclei_by_depth(allen_df,parameters_dict)
data = nma.analysis.generate_big_matrix(data,allen_df,parameters_dict)

del max_depth_dict,weight_dict
#============================================================================================================================
#============================================================================================================================ 
#%% 7 - Set plot parameters
plot_parameters = {'juci_cluster_names' : {'Cerebellum':'Cerebellum',
                                           'Medulla_Pons':'Medulla/Pons',
                                           'Midbrain':'Midbrain', 
                                           'Thalamus_Hypothalamus':'Thalamus/Hypothalamus',
                                           'Spinal_cord':'Spinal cord'},
                   'juci_cluster_colors' : {'Cerebellum':'red', 
                                            'Medulla_Pons':'blue',
                                            'Midbrain':'limegreen',  #yellowgreen
                                            'Thalamus_Hypothalamus':'gray',
                                            'Spinal_cord': 'gold'}} #whitesmoke

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
parameters_dict['pca_parameters'] = {'perplexity': 4,#3
                                     'learning_rate' : 10}#14#10#5#14
parameters_dict['umap_parameters'] = {'learning_rate': .4,#4
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
#%% 10 - hierarchical clustering (optional)
parameters_dict['clustering_method'] = 'weighted'#'ward')#'weighted')# 'average')#'complete')#pcs
data = nma.utils.plot.plot_clustering(data,parameters_dict,plot_parameters)
#============================================================================================================================
#============================================================================================================================ 
# =============================================================================
# #%% 10.5 - k-means clustering
# plot_parameters['kmeans']= {'cluster_num':5}
# from sklearn.cluster import KMeans
# inertias = []
# for k in range(1,10):
#     plot_parameters['kmeans']['cluster_num'] = k
#     kmeans = KMeans(n_clusters=plot_parameters['kmeans']['cluster_num'], init='k-means++', max_iter=300, n_init=10, random_state=0,verbose = True)
#     kmeans.fit(data['big_matrix'])
#     inertias.append(kmeans.inertia_)
# plt.figure()
# #%%
# plt.plot(range(1,10),inertias)
# 
# =============================================================================
#%% 10.5  -2 GMM clustering
from sklearn import mixture
from rastermap import Rastermap

plot_parameters['gmm']= {'cluster_num':2}
data_gmm = data['big_matrix_pcs'][:,:5]#data['big_matrix'][:,np.sum(data['big_matrix'],0)!=0]
gmm = mixture.GaussianMixture(n_components=plot_parameters['gmm']['cluster_num'], covariance_type='full').fit(data_gmm)
labels = gmm.predict(data_gmm)
probs = gmm.predict_proba(data_gmm)
model = Rastermap(n_components=1, n_X=40, nPC=20, init='pca')
model.fit(probs)
cellorder = model.isort
#%
plt.figure()
plt.imshow(probs[cellorder[::-1],:])
data['cell_order_gmm_clust'] = cellorder
for i,clustname in enumerate(data['cell_juci_clusters'][cellorder]):
    plt.plot(-1,len(data['cell_names'])-i-1,'o',color = plot_parameters['juci_cluster_colors'][clustname])
#%% 11 - projection pattern
#data['allen_axon_length_matrix'] = data['allen_axon_matrix']*data['axon_lengths'][:,None]/1000
plot_parameters['projection_pattern']={'separate_ipsi_contra_projections':False,
                                        'normalize_to_volume':False,
                                        'add_subtle_features':False,
                                        'base_nuclei' : ['Spinal cord','Cerebellar cortex','Cerebellar nuclei','Medulla','Pons','Midbrain','Thalamus','Hypothalamus','fiber tracts', 'ventricular systems'],
                                        'nuclei_ordering':'gmm_clustering',#'clustering',#'abc',#'rastermap','clustering'
                                        'projection_matrix_to_use': 'allen_axon_end_points_matrix', #'allen_axon_length_matrix',#*axon_end_point_numbers[:,None] #allen_axon_matrix*axon_lengths[:,None]/1000 #allen_axon_matrix_now*axon_lengths[:,None]/1000 #allen_axon_end_point_list_original#allen_axon_end_points_matrix#allen_axon_end_points_matrix*axon_end_point_numbers[:,None] #np.asarray(allen_axon_branch_point_list_original)#allen_axon_end_points_matrix'
                                        'use_log_scale':False,
                                        'vmin':0.000001,
                                        'gridcolor' :'black',#'white',
                                        'gridwidth' :.5,
                                        'style':'dark_background',#'dark_background', #'default'
                                        'colormap':plt.cm.RdPu_r,
                                        'colormap_set_under':'dimgray',#'white'
                                        'label_size' : 6
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
#%% add a colored cluster to the viewer
cluster_dict = nma.analysis.generate_sensory_motor_groups_with_projections(original_data)
for cluster_id in cluster_dict['cluster_names'].keys():
    data['cell_names'][cluster_dict['cluster_indices'][cluster_id]]
    subselect_parameters = {'needed_group' :  list(data['cell_names'][cluster_dict['cluster_indices'][cluster_id]]), #['AA1329'],#'Thalamus_Hypothalamus',
                            'soma_locations_needed' :None, #['Intermediate reticular nucleus', 'Medullary reticular nucleus'] or None
                            'axon_projections_needed':None #['VII','XII','AMB','V'] or None
                            }
    needed_cells = nma.analysis.generate_cell_list(allen_df,
                                                   original_data,
                                                   subselect_parameters['needed_group'],
                                                   soma_locations_needed=subselect_parameters['soma_locations_needed'], 
                                                   axon_projections_needed=subselect_parameters['axon_projections_needed']) 
    axon_color = 'limegreen'
    dendrite_color = cluster_dict['cluster_colors'][cluster_id]



    cells_added = {}
    for cell_to_add in needed_cells:
        nma.utils.snt.add_cell_to_viewer(snt_viewer_3d,
                                         jsonfolder = parameters_dict['json_dir'],
                                         cell_name = cell_to_add,
                                         cells_added = {},
                                         axon_color= axon_color,
                                         dendrite_color = dendrite_color,
                                         show_axon = False,
                                         show_dendrite = True,
                                         show_soma_location = False,
                                         soma_marker_size = 10,
                                         soma_marker_color = 'red',
                                         all_cells_on_left = False)
    

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
# 'yellow'  'gold'
# 'motor'
# 'soma location'   - also set 'soma_locations_needed'
# 'axon projection' - also set 'axon_projections_needed'
# 'Cerebellum': 'red'
# 'Thalamus_Hypothalamus': 'darkgray'
## 'Midbrain': 'limegreen'
# 'Medulla_Pons': 'blue'
# ['cellname1','cellname2'] just list the cells you want to add
# =============================================================================

subselect_parameters = {'needed_group' : 'Spinal_cord', #['AA1329'],#'Thalamus_Hypothalamus',
                        'soma_locations_needed' :None, #['Intermediate reticular nucleus', 'Medullary reticular nucleus'] or None
                        'axon_projections_needed':None #['VII','XII','AMB','V'] or None
                        }
needed_cells = nma.analysis.generate_cell_list(allen_df,
                                               original_data,
                                               subselect_parameters['needed_group'],
                                               soma_locations_needed=subselect_parameters['soma_locations_needed'], 
                                               axon_projections_needed=subselect_parameters['axon_projections_needed']) 
axon_color = 'limegreen'
dendrite_color = 'green'



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
                                     show_soma_location = False,
                                     soma_marker_size = 10,
                                     soma_marker_color = 'red',
                                     all_cells_on_left = False)

#%% 15 - add brain regions to the viewer - this is extremely slow for some reason...
brain_regions_to_add = [['Whole Brain',                   [86, 101, 115],    50],
                        ['Cerebellar cortex',             [220, 118, 51],    50],
                        ['Pyramus (VIII), Purkinje layer',[231, 76, 606],    50],
                        ['Medulla',                       [235, 222, 240],   50]]


for brain_region_i, brain_region in enumerate(brain_regions_to_add):
    try:
        structureId = allen_df.loc[allen_df['name'] == brain_region[0],'structureId'].values[0]
        allen_idx = (allen_df['structureId']==structureId).argmax()
        depth = allen_df['depth'][allen_idx]
        allenmesh = nma.utils.snt.AllenCompartment(structureId).getMesh()
        #color_0_1 = (depth - np.min(depth_range))/(np.max(depth_range) - np.min(depth_range))
        #color_now = np.asarray(256*np.asarray(cmap(color_0_1)))
        color_now=brain_region[1]
        allenmesh.setColor(nma.utils.snt.sntrgbcolor(color_now[0],color_now[1],color_now[2]),brain_region[2])
        snt_viewer_3d.addMesh(allenmesh)
    except:
        print('{} mesh was not found'.format(brain_region[0]))
#%% 16 - show soma location in CCF images
neuron_ID = 'AA1329'
try:
    allen_annotation_volume[0,0,0]
except:
    allen_annotation_volume = nma.utils.ccf.load_ccf_volume(parameters_dict)                           
jsonfile_with_path=os.path.join(parameters_dict['json_dir'],neuron_ID+'.json')
soma= nma.utils.snt.Tree(jsonfile_with_path, "soma")                              

xcoord =int(soma.getRoot().x/10)
zcoord = int(soma.getRoot().z/10)
ycoord = int(soma.getRoot().y/10)

fig = plt.figure(figsize = [15,15])
ax_sagittal = fig.add_subplot(211)
ax_coronal = fig.add_subplot(212)
img = allen_annotation_volume[xcoord,:,:].copy()
for i,val in enumerate(np.unique(img.flatten())):
    img[img==val]=i


im = ax_coronal.imshow(img,alpha = 1,cmap = 'Blues_r')#,cmap = 'jet'
#im.set_clim([0,i*2])
ax_coronal.plot(zcoord,ycoord,'ro')
img = allen_annotation_volume[:,:,zcoord].copy()
for i,val in enumerate(np.unique(img.flatten())):
    img[img==val]=i

im = ax_sagittal.imshow(img.T,alpha = 1,cmap = 'Blues_r')#,cmap = 'jet'
#im.set_clim([0,i*2])
ax_sagittal.plot(xcoord,ycoord,'ro')
ax_sagittal.set_title(neuron_ID)

#%% 17 axon/dendrite statistics for various groups
cluster_dict = nma.analysis.generate_main_projection_groups(original_data)
nma.utils.plot.plot_axon_dendrit_statistics(original_data, cluster_dict)


cluster_dict = nma.analysis.generate_soma_location_groups(original_data)
nma.utils.plot.plot_axon_dendrit_statistics(original_data, cluster_dict)


cluster_dict = nma.analysis.generate_projection_groups(original_data)
nma.utils.plot.plot_axon_dendrit_statistics(original_data, cluster_dict)

#%%
cluster_dict = nma.analysis.generate_sensory_motor_groups_with_projections(original_data)
nma.utils.plot.plot_axon_dendrit_statistics(original_data, cluster_dict)

#%% 18 - show dendrites in viewer

#% initialize viewer
snt_viewer_3d, snt_viewer_3d_cells_added = nma.utils.snt.init_viewer()
#%add cells to the viewer

for cell, soma_area in zip(data['cell_names'],data['soma_area_in_the_medulla']):
    subselect_parameters = {'needed_group' : [cell], #['AA1329'],#'Thalamus_Hypothalamus',
                            'soma_locations_needed' :None, #['Intermediate reticular nucleus', 'Medullary reticular nucleus'] or None
                            'axon_projections_needed':None #['VII','XII','AMB','V'] or None
                            }
    needed_cells = nma.analysis.generate_cell_list(allen_df,
                                                   original_data,
                                                   subselect_parameters['needed_group'],
                                                   soma_locations_needed=subselect_parameters['soma_locations_needed'], 
                                                   axon_projections_needed=subselect_parameters['axon_projections_needed']) 
    axon_color = 'limegreen'
    if 'sensory' in soma_area.lower():
        dendrite_color = 'green'
    elif 'motor' in soma_area.lower():
        dendrite_color = 'red'
    else:
        dendrite_color = 'yellow'

    cells_added = {}
    for cell_to_add in needed_cells:
        nma.utils.snt.add_cell_to_viewer(snt_viewer_3d,
                                         jsonfolder = parameters_dict['json_dir'],
                                         cell_name = cell_to_add,
                                         cells_added = {},
                                         axon_color= axon_color,
                                         dendrite_color = dendrite_color,
                                         show_axon = False,
                                         show_dendrite = True,
                                         show_soma_location = False,
                                         soma_marker_size = 10,
                                         soma_marker_color = 'red',
                                         all_cells_on_left = False)
        
        
        
#%% WTF????
cell_name = 'AA1329'
jsonfile_with_path = os.path.join(parameters_dict['json_dir'],cell_name+'.json')
axon = nma.snt.Tree(jsonfile_with_path, "axon")
dendrite = Tree(jsonfile_with_path, "dendrite")
dendritestats = nma.snt.TreeStatistics(dendrite)
branches = dendritestats.getBranches().toArray()   