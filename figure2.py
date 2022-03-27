#%%  Example cells based on axonal projection
# Dendrit thickness 4, axon 3
#%% figure 2.A Cerebellum projectiong cell



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

subselect_parameters = {'needed_group' : ['AA1344'], #['AA1329'],#'Thalamus_Hypothalamus',
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
                                     show_dendrite = True,
                                     show_soma_location = False,
                                     soma_marker_size = 10,
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


#%% figure 2.B Midbrain projectiong cell with nuclei



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

subselect_parameters = {'needed_group' : ['AA1343'], #['AA1329'],#'Thalamus_Hypothalamus',
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
                                     soma_marker_color = 'green',
                                     all_cells_on_left = False)

#%% 15 - add brain regions to the viewer - this is extremely slow for some reason...
brain_regions_to_add = [['Whole Brain',                   [86, 101, 115],    50],
                        ['Medulla',                       [220, 118, 51],    50],
                        ['Pons',                          [231, 76, 606],    50],
                        ['Cerebellum'                     [235, 222, 240],   50],
                        ['Midbrain'                       [], 50],
                        ['Thalamus'                       [], 50],
                        ['Hypothalamus'                   [], 50]]


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

#%% figure 2.C Medulla Pons projectiong cell with nuclei



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

subselect_parameters = {'needed_group' : ['AA1514'], #['AA1329'],#'Thalamus_Hypothalamus',
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
                                     show_dendrite = True,
                                     show_soma_location = False,
                                     soma_marker_size = 10,
                                     soma_marker_color = 'blue',
                                     all_cells_on_left = False)

#%% 15 - add brain regions to the viewer - this is extremely slow for some reason...
brain_regions_to_add = [['Whole Brain',                   [86, 101, 115],    50],
                        ['Medulla',                       [220, 118, 51],    50],
                        ['Pons',                          [231, 76, 606],    50],
                        ['Cerebellum'                     [235, 222, 240],   50],
                        ['Midbrain'                       [], 50],
                        ['Thalamus'                       [], 50],
                        ['Hypothalamus'                   [], 50]]


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


#%% figure 2.D Thalamus projectiong cell with nuclei



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

subselect_parameters = {'needed_group' : ['AA1511'], #['AA1329'],#'Thalamus_Hypothalamus',
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
                                     show_dendrite = True,
                                     show_soma_location = False,
                                     soma_marker_size = 10,
                                     soma_marker_color = 'gray',
                                     all_cells_on_left = False)

#%% 15 - add brain regions to the viewer - this is extremely slow for some reason...
brain_regions_to_add = [['Whole Brain',                   [86, 101, 115],    50],
                        ['Medulla',                       [220, 118, 51],    50],
                        ['Pons',                          [231, 76, 606],    50],
                        ['Cerebellum'                     [235, 222, 240],   50],
                        ['Midbrain'                       [], 50],
                        ['Thalamus'                       [], 50],
                        ['Hypothalamus'                   [], 50]]


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


#%% figure 2.E Spinal cord projectiong cell with nuclei



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

subselect_parameters = {'needed_group' : ['AA1517'], #['AA1329'],#'Thalamus_Hypothalamus',
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
                                     show_dendrite = True,
                                     show_soma_location = False,
                                     soma_marker_size = 10,
                                     soma_marker_color = 'orange',
                                     all_cells_on_left = False)

#%% 15 - add brain regions to the viewer - this is extremely slow for some reason...
brain_regions_to_add = [['Whole Brain',                   [86, 101, 115],    50],
                        ['Medulla',                       [220, 118, 51],    50],
                        ['Pons',                          [231, 76, 606],    50],
                        ['Cerebellum'                     [235, 222, 240],   50],
                        ['Midbrain'                       [], 50],
                        ['Thalamus'                       [], 50],
                        ['Hypothalamus'                   [], 50]]


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


#%% figure 2.F Spinal cord projectiong cell with nuclei




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

subselect_parameters = {'needed_group' : 'all', #['AA1329'],#'Thalamus_Hypothalamus',
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
                                     show_axon = False,
                                     show_dendrite = True,
                                     show_soma_location = False,
                                     soma_marker_size = 10,
                                     soma_marker_color = 'red',
                                     all_cells_on_left = False)

#%% 15 - add brain regions to the viewer - this is extremely slow for some reason...
brain_regions_to_add = [['Whole Brain',                   [86, 101, 115],    50],
                        ['Medulla',                       [220, 118, 51],    50],
                        ['Pons',                          [231, 76, 606],    50],
                        ['Cerebellum'                     [235, 222, 240],   50],
                        ['Midbrain'                       [], 50],
                        ['Thalamus'                       [], 50],
                        ['Hypothalamus'                   [], 50]]


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





