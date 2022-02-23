import numpy as np
import pandas as pd
import nrrd
import h5py
from scyjava import jimport
AllenCompartment = jimport('sc.fiji.snt.annotation.AllenCompartment')#, AllenUtils, VFBUtils, ZBAtlasUtils)
AllenUtils = jimport('sc.fiji.snt.annotation.AllenUtils')#, , VFBUtils, ZBAtlasUtils)

def generate_list_of_nuclei_by_depth(allen_df,parameters_dict):
    max_depth_dict = parameters_dict['maximum_depth_dict']
    base_brainregion_names = list()
    base_brainregion_ids = list()
    base_brainregion_max_depths = list()
    parentids = list()
    for brainregion_now in max_depth_dict.keys():
        base_brainregion_ids.append(allen_df.loc[allen_df['name']==brainregion_now,'structureId'].values[0])
        parentids.extend(np.asarray(allen_df.loc[allen_df['name']==brainregion_now,'structureIdPath'].values[0].strip('/').split('/'),int))
        base_brainregion_names.append(brainregion_now)
        base_brainregion_max_depths.append(max_depth_dict[brainregion_now])
    parentids = np.unique(parentids)
    allen_df_needed_nucleus = list()
    for allen_row in allen_df.iterrows():
        allen_row = allen_row[1]
        id_path = np.asarray(allen_row['structureIdPath'].strip('/').split('/'),int)
        if 8 not in id_path:
            allen_df_needed_nucleus.append(False)
            continue
        if id_path[-1] in parentids:
            allen_df_needed_nucleus.append(True)
            continue
        for base_brainregion_id,base_brainregion_name,base_brainregion_max_depth in zip(base_brainregion_ids,base_brainregion_names,base_brainregion_max_depths):
            brainregion_found = False
            if base_brainregion_id in id_path:
                if allen_row['depth'] <= base_brainregion_max_depth:
                    allen_df_needed_nucleus.append(True)
                else:
                    allen_df_needed_nucleus.append(False)
                brainregion_found = True
                break
        if brainregion_found:
            continue
        allen_df_needed_nucleus.append(False)

    allen_needed_nuclei = np.concatenate([allen_df_needed_nucleus,allen_df_needed_nucleus])
    return allen_needed_nuclei

def calculate_and_save_ccf_volumes(allen_df,
                               allen_annotation_volume,
                               savefile = r'C:\Users\judith.baka\OneDrive - Allen Institute\Desktop\MouseLight\Allen_p56_mouse_annotation\v2.5\allen_df_with_volumes.csv'):
    """
    Calculate volumes of each brain region in the CCF by summing up voxels. Then adds them to the allen dataframe and saves it as a .csv file
    Parameters
    ----------
    allen_df : pandas dataframe
    allen_annotation_volume : ndarray float
        3d volume of the mouse brain showing brain regions
    savefile : str, optional
        where output should be saved.

    Returns
    -------
    None.

    """
    
    brain_region_volumes = np.zeros(len(allen_df))
    for z,plane in enumerate(allen_annotation_volume):
        print(z)
        plane = plane.flatten()
        unique_ids = np.unique(plane)
        for region_i,allen_df_row in enumerate(allen_df.iterrows()):
    #        print(region_i)
            node_id = allen_df_row[1]['structureId']
            if node_id not in unique_ids:
                continue
            volume_now = sum(plane == node_id)
            #if volume_now>0:
            node_ids = np.asarray(allen_df.loc[allen_df['structureId']==node_id,'structureIdPath'].values[0].strip('/').split('/'),int)
            for node_id in node_ids:
                brain_region_volumes[np.argmax(allen_df['structureId']==node_id)] = brain_region_volumes[np.argmax(allen_df['structureId']==node_id)] + volume_now
    
    allen_df_with_volumes = allen_df.assign(volume = brain_region_volumes/1000)
    allen_df_with_volumes.to_csv(savefile)
    
def load_allen_dataframe(parameters_dict):
    """
    Loads our constructs the Allen dataframe that contains all the ontologies

    Parameters
    ----------
    parameters_dict : dict
        dictionary containing all the parameters

    Returns
    -------
    allen_df :dataframe
        The allen dataframe

    """
    if parameters_dict['load_volumes_with_allen_df']:
        if parameters_dict['ccf_version'] == '3.0':
            allen_df = pd.read_csv(parameters_dict['path_allen_df_with_volumes'])
        elif parameters_dict['ccf_version'] == '2.5':
            allen_df = pd.read_csv(parameters_dict['path_allen_df_with_volumes'])
        else:
            print('UNKNOWN CCF version {}'.format(parameters_dict['ccf_version']))
    else:
        #%
        header = ['acronym','structureId','parentStructureId','depth','name','structureIdPath']
        allen_dict = dict()
        for head in header: allen_dict[head] = list()
        areas = AllenUtils.getOntologies()
        for area in areas:#.iterator():
            structureidpath = ''
            ancestorlist = list()
            try:
                parentid = int(area.getAncestor(1).id())
                ancestors = area.getAncestors()
                for ancestor in ancestors: ancestorlist.append(ancestor.id()); structureidpath += '/{}'.format(ancestor.id())
            except:
                print('{} has no parents'.format(area.name()))
                parentid = None
            structureidpath += '/{}'.format(int(area.id()))
            allen_dict['acronym'].append(str(area.acronym()))
            allen_dict['structureId'].append(int(area.id()))
            allen_dict['parentStructureId'].append(parentid)
            allen_dict['depth'].append(len(ancestorlist))
            allen_dict['name'].append(str(area.name()))
            allen_dict['structureIdPath'].append(structureidpath)
            #break
        allen_df = pd.DataFrame.from_dict(allen_dict)
    return allen_df

def load_ccf_volume(parameters_dict):
    if float(parameters_dict['ccf_version']) == 3.:
        allen_annotation_volume, allen_volume_header_v3 = nrrd.read(r'C:\Users\judith.baka\OneDrive - Allen Institute\Desktop\MouseLight\Allen_p56_mouse_annotation\annotation_10.nrrd')    #3.0 v annotation_10.nrrd
    elif float(parameters_dict['ccf_version']) == 2.5:
        #%
        f = h5py.File(r'C:\Users\judith.baka\OneDrive - Allen Institute\Desktop\MouseLight\Allen_p56_mouse_annotation\v2.5\OntologyAtlas.h5', 'r')
        allen_annotation_volume = f['OntologyAtlas']
        allen_annotation_volume = allen_annotation_volume[:]
        f.close()
        #%
    else:
        print('UNKNOWN CCF version {}'.format(parameters_dict['ccf_version']))  
    return allen_annotation_volume

#%  defining final nuclei for  axon endpoint search
def define_final_nuclei_for_endpoint_search(allen_df,parameters_dict):
    minimum_depth_dict = parameters_dict['minimum_depth_dict']
    base_brainregion_names = list()
    base_brainregion_ids = list()
    base_brainregion_min_depths = list()
    for brainregion_now in minimum_depth_dict.keys():
        base_brainregion_ids.append(allen_df.loc[allen_df['name']==brainregion_now,'structureId'].values[0])
        base_brainregion_names.append(brainregion_now)
        base_brainregion_min_depths.append(minimum_depth_dict[brainregion_now])
    allen_df_final_nucleus = list()
    for allen_row in allen_df.iterrows():
        allen_row = allen_row[1]
        id_path = np.asarray(allen_row['structureIdPath'].strip('/').split('/'),int)
        if 8 not in id_path:
            allen_df_final_nucleus.append(False)
            continue
        if not any(allen_df['parentStructureId'].values==allen_row['structureId']):
            allen_df_final_nucleus.append(True)
            continue
        for base_brainregion_id,base_brainregion_name,base_brainregion_min_depth in zip(base_brainregion_ids,base_brainregion_names,base_brainregion_min_depths):
            brainregion_found = False
            if base_brainregion_id in id_path:
                if allen_row['depth'] >= base_brainregion_min_depth:
                    allen_df_final_nucleus.append(True)
                else:
                    allen_df_final_nucleus.append(False)
                brainregion_found = True
                break
        if brainregion_found:
            continue
        allen_df_final_nucleus.append(False)
        
    allen_df['final_nucleus']= allen_df_final_nucleus
    allen_df = allen_df.reset_index() 
    return allen_df