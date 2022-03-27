#%%  Axonal group

#%% figure 4.A Cerebellum projecting group



#%% 13 - initialize viewer
snt_viewer_3d, snt_viewer_3d_cells_added = nma.utils.snt.init_viewer()

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

subselect_parameters = {'needed_group' : 'Cerebellum', #['AA1329'],#'Thalamus_Hypothalamus',
                        'soma_locations_needed' :None, #['Intermediate reticular nucleus', 'Medullary reticular nucleus'] or None
                        'axon_projections_needed':None #['VII','XII','AMB','V'] or None
                        }
needed_cells = nma.analysis.generate_cell_list(allen_df,
                                               original_data,
                                               subselect_parameters['needed_group'],
                                               soma_locations_needed=subselect_parameters['soma_locations_needed'], 
                                               axon_projections_needed=subselect_parameters['axon_projections_needed']) 
axon_color = 'red'
dendrite_color = 'darkred'



cells_added = {}
for cell_to_add in needed_cells:
    nma.utils.snt.add_cell_to_viewer(snt_viewer_3d,
                                     jsonfolder = parameters_dict['json_dir'],
                                     cell_name = cell_to_add,
                                     cells_added = {},
                                     axon_color= axon_color,
                                     dendrite_color = dendrite_color,
                                     show_axon = True,
                                     show_dendrite = False,
                                     show_soma_location = True,
                                     soma_marker_size = 5,
                                     soma_marker_color = 'red',
                                     all_cells_on_left = False)

#%% 15 - add brain regions to the viewer - this is extremely slow for some reason...
brain_regions_to_add = [['Whole Brain',                   [86, 101, 115],    98],
                        ['Medulla',                       [207, 215, 255],    85],
                        ['Pons',                          [225, 216, 253],    85],
                        ['Cerebellum',                    [255, 234, 249],   85],
                        ['Midbrain',                      [204, 255, 153], 85],
                        ['Thalamus',                      [229, 214, 157], 90],
                        ['Hypothalamus',                  [233, 215, 105], 90]]


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



#=============================================================================


#%% figure 4.B Medulla Pons projectiong group



#%% 13 - initialize viewer
snt_viewer_3d, snt_viewer_3d_cells_added = nma.utils.snt.init_viewer()


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

subselect_parameters = {'needed_group' : 'Medulla_Pons', #['AA1329'],#'Thalamus_Hypothalamus',
                        'soma_locations_needed' :None, #['Intermediate reticular nucleus', 'Medullary reticular nucleus'] or None
                        'axon_projections_needed':None #['VII','XII','AMB','V'] or None
                        }
needed_cells = nma.analysis.generate_cell_list(allen_df,
                                               original_data,
                                               subselect_parameters['needed_group'],
                                               soma_locations_needed=subselect_parameters['soma_locations_needed'], 
                                               axon_projections_needed=subselect_parameters['axon_projections_needed']) 
axon_color = 'blue'
dendrite_color = 'darkblue'



cells_added = {}
for cell_to_add in needed_cells:
    nma.utils.snt.add_cell_to_viewer(snt_viewer_3d,
                                     jsonfolder = parameters_dict['json_dir'],
                                     cell_name = cell_to_add,
                                     cells_added = {},
                                     axon_color= axon_color,
                                     dendrite_color = dendrite_color,
                                     show_axon = True,
                                     show_dendrite = False,
                                     show_soma_location = True,
                                     soma_marker_size = 5,
                                     soma_marker_color = 'blue',
                                     all_cells_on_left = False)

#%% 15 - add brain regions to the viewer - this is extremely slow for some reason...
brain_regions_to_add = [['Whole Brain',                   [86, 101, 115],    98],
                        ['Medulla',                       [207, 215, 255],    85],
                        ['Pons',                          [225, 216, 253],    85],
                        ['Cerebellum',                    [255, 234, 249],   85],
                        ['Midbrain',                      [204, 255, 153], 85],
                        ['Thalamus',                      [229, 214, 157], 90],
                        ['Hypothalamus',                  [233, 215, 105], 90]]


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



#=============================================================================


#%% figure 4.C Midbrain projectiong group



#%% 13 - initialize viewer
snt_viewer_3d, snt_viewer_3d_cells_added = nma.utils.snt.init_viewer()


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

subselect_parameters = {'needed_group' : 'Midbrain', #['AA1329'],#'Thalamus_Hypothalamus',
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
                                     show_dendrite = False,
                                     show_soma_location = True,
                                     soma_marker_size = 5,
                                     soma_marker_color = 'limegreen',
                                     all_cells_on_left = False)

#%% 15 - add brain regions to the viewer - this is extremely slow for some reason...
brain_regions_to_add = [['Whole Brain',                   [86, 101, 115],    98],
                        ['Medulla',                       [207, 215, 255],    85],
                        ['Pons',                          [225, 216, 253],    85],
                        ['Cerebellum',                    [255, 234, 249],   85],
                        ['Midbrain',                      [204, 255, 153], 85],
                        ['Thalamus',                      [229, 214, 157], 90],
                        ['Hypothalamus',                  [233, 215, 105], 90]]


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



#=============================================================================


#%% figure 4.D Thalamus Hypothalamus projectiong group



#%% 13 - initialize viewer
snt_viewer_3d, snt_viewer_3d_cells_added = nma.utils.snt.init_viewer()


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

subselect_parameters = {'needed_group' : 'Thalamus_Hypothalamus', #['AA1329'],#'Thalamus_Hypothalamus',
                        'soma_locations_needed' :None, #['Intermediate reticular nucleus', 'Medullary reticular nucleus'] or None
                        'axon_projections_needed':None #['VII','XII','AMB','V'] or None
                        }
needed_cells = nma.analysis.generate_cell_list(allen_df,
                                               original_data,
                                               subselect_parameters['needed_group'],
                                               soma_locations_needed=subselect_parameters['soma_locations_needed'], 
                                               axon_projections_needed=subselect_parameters['axon_projections_needed']) 
axon_color = 'grey'
dendrite_color = 'black'



cells_added = {}
for cell_to_add in needed_cells:
    nma.utils.snt.add_cell_to_viewer(snt_viewer_3d,
                                     jsonfolder = parameters_dict['json_dir'],
                                     cell_name = cell_to_add,
                                     cells_added = {},
                                     axon_color= axon_color,
                                     dendrite_color = dendrite_color,
                                     show_axon = True,
                                     show_dendrite = False,
                                     show_soma_location = True,
                                     soma_marker_size = 5,
                                     soma_marker_color = 'grey',
                                     all_cells_on_left = False)

#%% 15 - add brain regions to the viewer - this is extremely slow for some reason...
brain_regions_to_add = [['Whole Brain',                   [86, 101, 115],    98],
                        ['Medulla',                       [207, 215, 255],    85],
                        ['Pons',                          [225, 216, 253],    85],
                        ['Cerebellum',                    [255, 234, 249],   85],
                        ['Midbrain',                      [204, 255, 153], 85],
                        ['Thalamus',                      [229, 214, 157], 90],
                        ['Hypothalamus',                  [233, 215, 105], 90]]


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



#=============================================================================


#%% figure 4.E Spinal cord projecting group




#%% 13 - initialize viewer
snt_viewer_3d, snt_viewer_3d_cells_added = nma.utils.snt.init_viewer()


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

subselect_parameters = {'needed_group' : 'yellow', #['AA1329'],#'Thalamus_Hypothalamus',
                        'soma_locations_needed' :None, #['Intermediate reticular nucleus', 'Medullary reticular nucleus'] or None
                        'axon_projections_needed':None #['VII','XII','AMB','V'] or None
                        }
needed_cells = nma.analysis.generate_cell_list(allen_df,
                                               original_data,
                                               subselect_parameters['needed_group'],
                                               soma_locations_needed=subselect_parameters['soma_locations_needed'], 
                                               axon_projections_needed=subselect_parameters['axon_projections_needed']) 
axon_color = 'gold'
dendrite_color = 'orange'



cells_added = {}
for cell_to_add in needed_cells:
    nma.utils.snt.add_cell_to_viewer(snt_viewer_3d,
                                     jsonfolder = parameters_dict['json_dir'],
                                     cell_name = cell_to_add,
                                     cells_added = {},
                                     axon_color= axon_color,
                                     dendrite_color = dendrite_color,
                                     show_axon = True,
                                     show_dendrite = False,
                                     show_soma_location = True,
                                     soma_marker_size = 5,
                                     soma_marker_color = 'gold',
                                     all_cells_on_left = False)

#%% 15 - add brain regions to the viewer - this is extremely slow for some reason...
brain_regions_to_add = [['Whole Brain',                   [86, 101, 115],    98],
                        ['Medulla',                       [207, 215, 255],    85],
                        ['Pons',                          [225, 216, 253],    85],
                        ['Cerebellum',                    [255, 234, 249],   85],
                        ['Midbrain',                      [204, 255, 153], 85],
                        ['Thalamus',                      [229, 214, 157], 90],
                        ['Hypothalamus',                  [233, 215, 105], 90]]


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



#=============================================================================


#%% figure 4.F  All group together

#%% 13 - initialize viewer
snt_viewer_3d, snt_viewer_3d_cells_added = nma.utils.snt.init_viewer()


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

subselect_parameters = {'needed_group' : 'Cerebellum', #['AA1329'],#'Thalamus_Hypothalamus',
                        'soma_locations_needed' :None, #['Intermediate reticular nucleus', 'Medullary reticular nucleus'] or None
                        'axon_projections_needed':None #['VII','XII','AMB','V'] or None
                        }
needed_cells = nma.analysis.generate_cell_list(allen_df,
                                               original_data,
                                               subselect_parameters['needed_group'],
                                               soma_locations_needed=subselect_parameters['soma_locations_needed'], 
                                               axon_projections_needed=subselect_parameters['axon_projections_needed']) 
axon_color = 'red'
dendrite_color = 'red'



cells_added = {}
for cell_to_add in needed_cells:
    nma.utils.snt.add_cell_to_viewer(snt_viewer_3d,
                                     jsonfolder = parameters_dict['json_dir'],
                                     cell_name = cell_to_add,
                                     cells_added = {},
                                     axon_color= axon_color,
                                     dendrite_color = dendrite_color,
                                     show_axon = True,
                                     show_dendrite = False,
                                     show_soma_location = True,
                                     soma_marker_size = 5,
                                     soma_marker_color = 'red',
                                     all_cells_on_left = False)

#%% 15 - add brain regions to the viewer - this is extremely slow for some reason...
brain_regions_to_add = [['Whole Brain',                   [86, 101, 115],    98],
                        ['Medulla',                       [207, 215, 255],    85],
                        ['Pons',                          [225, 216, 253],    85],
                        ['Cerebellum',                    [255, 234, 249],   85],
                        ['Midbrain',                      [204, 255, 153], 85],
                        ['Thalamus',                      [229, 214, 157], 90],
                        ['Hypothalamus',                  [233, 215, 105], 90]]


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



#=============================================================================

#============================================================================================================================ 
#%% 11 - projection pattern
#data['allen_axon_length_matrix'] = data['allen_axon_matrix']*data['axon_lengths'][:,None]/1000
plot_parameters['projection_pattern']={'separate_ipsi_contra_projections':False,
                                        'normalize_to_volume':False,
                                        'add_subtle_features':False,
                                        'base_nuclei' : ['Cerebellar cortex','Cerebellar nuclei','Medulla','Pons','Midbrain','Thalamus','Hypothalamus','fiber tracts', 'ventricular systems'],
                                        'nuclei_ordering':'clustering',#'abc',#'rastermap','clustering'
                                        'projection_matrix_to_use': 'allen_axon_end_points_matrix', #'allen_axon_length_matrix',#*axon_end_point_numbers[:,None] #allen_axon_matrix*axon_lengths[:,None]/1000 #allen_axon_matrix_now*axon_lengths[:,None]/1000 #allen_axon_end_point_list_original#allen_axon_end_points_matrix#allen_axon_end_points_matrix*axon_end_point_numbers[:,None] #np.asarray(allen_axon_branch_point_list_original)#allen_axon_end_points_matrix'
                                        'use_log_scale':False,
                                        'vmin':0.000001,
                                        'gridcolor' : 'white',#'black',
                                        'gridwidth' : .5,
                                        'style':'default', #'dark_background', #'default'
                                        'colormap': plt.cm.viridis_r,#magma_r,
                                        'colormap_set_under':'white'
                                        }
nma.utils.plot.plot_projection_patterns(data, allen_df, plot_parameters)
#============================================================================================================================
#============================================================================================================================ 
