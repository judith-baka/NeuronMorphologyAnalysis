import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from rastermap import Rastermap
from sklearn.cluster import KMeans

def plot_params(paramsdict,ax_parameters,ncolumns=6):
    plt.style.use('default')
    ncolumns = 7
    nrows = np.ceil(len(paramsdict)/ncolumns)
    params_table=list()
    params_row = list()
    sumweight = 0
    for key in paramsdict.keys():
        if 'weight' in key:
           sumweight +=  paramsdict[key]
    for key in paramsdict.keys():
        if len(params_row)==nrows:
            params_table.append(params_row)
            params_row = list()
        params_row.append(key+': - {:.0f}%'.format(paramsdict[key]/sumweight*100))
    if len(params_row)>0:
        while len(params_row)<nrows:
            params_row.append('')
    params_table.append(params_row)        
            
    table = ax_parameters.table(cellText=params_table, rowLoc = 'center',colLoc = 'center',edges = 'open',loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.1, 1.1)
    ax_parameters.axis('off')
    
def plot_pca(data,parameters_dict,plot_parameters,cluster_dict = None):
    plt.style.use('default') #atirni dark_background-ra ha fekete hatteret akarok
    pc0_idx = 0
    pc1_idx = 1
    pc2_idx = 2
    fig = plt.figure(figsize = [10,15])
    ax_pca = fig.add_subplot(311, projection='3d')
    ax_pca_var = fig.add_subplot(312)
    ax_parameters = fig.add_subplot(313)
    if cluster_dict == None:
        for juci_cluster in plot_parameters['juci_cluster_colors'].keys():
            color = plot_parameters['juci_cluster_colors'][juci_cluster]
            name = plot_parameters['juci_cluster_names'][juci_cluster]
            idx = (np.asarray(data['cell_juci_clusters']) == juci_cluster) & (np.array(data['soma_area_in_the_medulla'],str) == 'Motor related')
            ax_pca.plot(data['big_matrix_pcs'][idx,pc0_idx],data['big_matrix_pcs'][idx,pc1_idx],data['big_matrix_pcs'][idx,pc2_idx], 'o',color = color,label=name, alpha = 1) #+' projecting, motor related'
            idx = (np.asarray(data['cell_juci_clusters']) == juci_cluster) & (np.array(data['soma_area_in_the_medulla'],str) == 'Sensory related')
            ax_pca.plot(data['big_matrix_pcs'][idx,pc0_idx],data['big_matrix_pcs'][idx,pc1_idx],data['big_matrix_pcs'][idx,pc2_idx], 'o', color = color, alpha = .25)#label=juci_cluster+' projecting, sensory related',
    else:
        for cluster_id in cluster_dict['cluster_names'].keys():
            color = cluster_dict['cluster_colors'][cluster_id]
            name = cluster_dict['cluster_names'][cluster_id]
            idx = cluster_dict['cluster_indices'][cluster_id]
            ax_pca.plot(data['big_matrix_pcs'][idx,pc0_idx],data['big_matrix_pcs'][idx,pc1_idx],data['big_matrix_pcs'][idx,pc2_idx], 'o',color = color,label=name, alpha = 1) #+' projecting, motor related'
    ax_pca.legend()
    ax_pca.set_xlabel('PC{}'.format(pc0_idx+1))
    ax_pca.set_ylabel('PC{}'.format(pc1_idx+1))
    ax_pca.set_zlabel('PC{}'.format(pc2_idx+1))
    ax_pca_var.plot(np.cumsum(data['big_matrix_pca'].explained_variance_ratio_),'ko-')
    ax_pca_var.set_ylabel('Cumulative explained variance')
    ax_pca_var.set_xlabel('PCs')
    ax_pca_var.set_ylim([0,1])
    plot_params(parameters_dict['weights_dict'],ax_parameters)
    
def plot_tsne_umap(data,parameters_dict,plot_parameters,cluster_dict = None):
    plt.style.use('default')
    fig_tsne_umap = plt.figure(figsize = [10,15])
    ax_tsne = fig_tsne_umap.add_subplot(311)
    
    if cluster_dict == None:
        for juci_cluster in plot_parameters['juci_cluster_colors'].keys():
            color = plot_parameters['juci_cluster_colors'][juci_cluster]
            name = plot_parameters['juci_cluster_names'][juci_cluster]
            idx = (np.asarray(data['cell_juci_clusters']) == juci_cluster) & (np.array(data['soma_area_in_the_medulla'],str) == 'Motor related')
            ax_tsne.plot(data['big_matrix_tsne'][idx,0],data['big_matrix_tsne'][idx,1],'o',color = color,label=name,ms=12, alpha = 1)#+' projecting, motor related'
            idx = (np.asarray(data['cell_juci_clusters']) == juci_cluster) & (np.array(data['soma_area_in_the_medulla'],str) == 'Sensory related')
            ax_tsne.plot(data['big_matrix_tsne'][idx,0],data['big_matrix_tsne'][idx,1],'o',color = color,ms=12, alpha = .25) #label=juci_cluster+' projecting, sensory related',
    else:
        for cluster_id in cluster_dict['cluster_names'].keys():
            color = cluster_dict['cluster_colors'][cluster_id]
            name = cluster_dict['cluster_names'][cluster_id]
            idx = cluster_dict['cluster_indices'][cluster_id]
            ax_tsne.plot(data['big_matrix_tsne'][idx,0],data['big_matrix_tsne'][idx,1],'o',color = color,label=name,ms=12, alpha = 1)#+' projecting, motor related'
            
    for i,txt in enumerate(data['cell_names']):
        ax_tsne.annotate(txt,(data['big_matrix_tsne'][i,0],data['big_matrix_tsne'][i,1]))
    ax_tsne.set_title('tSNE - perplexity: {}, learning rate: {}'.format(parameters_dict['pca_parameters']['perplexity'],parameters_dict['pca_parameters']['learning_rate']))    
    plt.legend()
    
    ax_umap = fig_tsne_umap.add_subplot(312)
    if cluster_dict == None:
        for juci_cluster in plot_parameters['juci_cluster_colors'].keys():
            color = plot_parameters['juci_cluster_colors'][juci_cluster]
            name = plot_parameters['juci_cluster_names'][juci_cluster]
            idx = (np.asarray(data['cell_juci_clusters']) == juci_cluster) & (np.array(data['soma_area_in_the_medulla'],str) == 'Motor related')
            ax_umap.plot(data['big_matrix_umap'][idx,0],data['big_matrix_umap'][idx,1],'o',color = color,label=name,ms=12, alpha = 1)#+' projecting, motor related'
            idx = (np.asarray(data['cell_juci_clusters']) == juci_cluster) & (np.array(data['soma_area_in_the_medulla'],str) == 'Sensory related')
            ax_umap.plot(data['big_matrix_umap'][idx,0],data['big_matrix_umap'][idx,1],'o',color = color,ms=12, alpha = .25) #label=juci_cluster+' projecting, sensory related',
    else:
        for cluster_id in cluster_dict['cluster_names'].keys():
            color = cluster_dict['cluster_colors'][cluster_id]
            name = cluster_dict['cluster_names'][cluster_id]
            idx = cluster_dict['cluster_indices'][cluster_id]
            ax_umap.plot(data['big_matrix_umap'][idx,0],data['big_matrix_umap'][idx,1],'o',color = color,label=name,ms=12, alpha = 1)#+' projecting, motor related'
        
    for i,txt in enumerate(data['cell_names']):
        ax_umap.annotate(txt,(data['big_matrix_umap'][i,0],data['big_matrix_umap'][i,1]))
    ax_umap.set_title('uMAP - neighbors: {}, min_dist: {}, learning rate: {}'.format(parameters_dict['umap_parameters']['n_neigbors'],
                                                                                     parameters_dict['umap_parameters']['min_dist'],
                                                                                     parameters_dict['umap_parameters']['learning_rate']))        

    plt.legend()
    ax_parameters = fig_tsne_umap.add_subplot(313)
    plot_params(parameters_dict['weights_dict'],ax_parameters)

def plot_clustering(data,parameters_dict,plot_parameters):
    

    plt.style.use('dark_background')
    hierarchy.set_link_color_palette(['lightcoral','dodgerblue','yellowgreen',])
    linked = linkage(data['big_matrix'],parameters_dict['clustering_method'])

    fig_cluster = plt.figure(figsize=(15, 15))
    ax_clust = fig_cluster.add_subplot(111)
    R = dendrogram(linked,
                   orientation='top',
                   labels=data['cell_names'],
                   color_threshold=4)

    #%
    cell_order_clust = list()
    for i,cell_ID in enumerate(R['ivl']):
        cell_idx = np.asarray(data['cell_names']) == cell_ID
        cell_order_clust.append(np.where(cell_idx)[0][0])
        try:
            color = plot_parameters['juci_cluster_colors'][np.asarray(data['cell_juci_clusters'])[np.argmax(cell_idx)]]
        except:
            print(cell_idx)
        motor_vs_sensory = data['soma_area_in_the_medulla'][np.argmax(cell_idx)]
        if motor_vs_sensory == 'Sensory related':
            marker = 'o'
        elif motor_vs_sensory == 'Motor related':
            marker = '^'
        else:
            marker = 'v'
        ax_clust.plot(5+i*10,0,marker,color = color,markersize = 14)
    ax_clust.tick_params(axis='x', which='major', labelsize=8) 
    cell_order_clust = np.asarray(cell_order_clust)
    plt.grid(False)
    data['cell_order_clust'] = cell_order_clust
    return data

# =============================================================================
# def plot_kmeans(data,parameters_dict,plot_parameters):
#     #%%
#     kmeans = KMeans(n_clusters=plot_parameters['cluster_num'], init='k-means++', max_iter=300, n_init=10, random_state=0)
#     kmeans.fit(data['big_matrix'])
# =============================================================================
#%%
def plot_projection_patterns(data,allen_df,plot_parameters):
    cmap = plot_parameters['projection_pattern']['colormap']
    cmap.set_under(color=plot_parameters['projection_pattern']['colormap_set_under'])
    plt.style.use(plot_parameters['projection_pattern']['style'])

    projection_matrix_to_use = data[plot_parameters['projection_pattern']['projection_matrix_to_use']]
    if plot_parameters['projection_pattern']['normalize_to_volume']:
        volumes = np.concatenate([allen_df['volume'].values]*2)
        projection_matrix_to_use = projection_matrix_to_use /volumes
    if plot_parameters['projection_pattern']['use_log_scale']:
        projection_matrix_to_use = np.log(projection_matrix_to_use*10000000) # TODO what's this?
        projection_matrix_to_use[projection_matrix_to_use <0] = 0


    cellname_nucleus = list()
    soma_locations_acronyms = list()
    try:
        for cellname,nucleus in zip(data['cell_names'],data['soma_locations']):
            acronym = allen_df.loc[allen_df['name']==nucleus,'acronym'].values[0]
            cellname_nucleus.append('{}   -   {}'.format(acronym,cellname))
            soma_locations_acronyms.append(acronym)
    except:
        cellname_nucleus = data['cell_names']
        soma_locations_acronyms = data['cell_names']
    if plot_parameters['projection_pattern']['nuclei_ordering'] == 'abc':
        cellorder = np.argsort(cellname_nucleus)
    elif plot_parameters['projection_pattern']['nuclei_ordering'] == 'rastermap':
        model = Rastermap(n_components=1, n_X=40, nPC=20, init='pca')
        model.fit(projection_matrix_to_use)
        cellorder = model.isort
    elif plot_parameters['projection_pattern']['nuclei_ordering'] == 'clustering':
        cellorder = data['cell_order_clust']#cellorder_manual
    elif plot_parameters['projection_pattern']['nuclei_ordering'] == 'gmm_clustering':
        cellorder = data['cell_order_gmm_clust']#cellorder_manual
    else:
        cellorder =np.arange(len(cellname_nucleus))
        

    cellname_nucleus = np.asarray(cellname_nucleus)
    base_nuclei_idx = list()
    base_nuclei_structureId = list()
    header = list()
    projection_pattern= list()
    for base_nucleus in plot_parameters['projection_pattern']['base_nuclei']:
        base_nucleus_idx_ipsi = (allen_df['name']==base_nucleus).argmax()
        base_nucleus_idx_contra = base_nucleus_idx_ipsi+len(allen_df)
        base_nucleus_structureId = allen_df['structureId'][base_nucleus_idx_ipsi]
        base_nuclei_idx.append(base_nucleus_idx_ipsi)
        base_nuclei_structureId.append(base_nucleus_structureId)
        if plot_parameters['projection_pattern']['separate_ipsi_contra_projections']:
            header.append('{}-Y'.format(base_nucleus))
            header.append('{}-C'.format(base_nucleus))
            projection_pattern.append(projection_matrix_to_use[:,base_nucleus_idx_ipsi])#+projection_matrix_to_use[:,base_nucleus_idx_contra]
            projection_pattern.append(projection_matrix_to_use[:,base_nucleus_idx_contra])#+projection_matrix_to_use[:,base_nucleus_idx_contra]
        else:
            header.append('{}'.format(base_nucleus))
            projection_pattern.append(projection_matrix_to_use[:,base_nucleus_idx_ipsi]+projection_matrix_to_use[:,base_nucleus_idx_contra])
    if plot_parameters['projection_pattern']['add_subtle_features']:
        projection_pattern_now = np.concatenate([np.asarray(projection_pattern),data['subtle_features_matrix'].T],0)
        header_now = np.concatenate([header,data['subtle_features_header']])
    else:
        projection_pattern_now = projection_pattern
        header_now = header

    fig = plt.figure()
    ax=fig.add_subplot(111)
    im = ax.pcolormesh(np.asarray(projection_pattern_now).T[cellorder,:],
                       cmap = cmap,
                       edgecolor= plot_parameters['projection_pattern']['gridcolor'],
                       linewidth=plot_parameters['projection_pattern']['gridwidth'],
                       vmin=plot_parameters['projection_pattern']['vmin']) #ax.imshow(np.asarray(projection_pattern).T[cellorder,:],aspect='auto',cmap = 'RdPu')
    ax.set_yticks(np.arange(0,len(cellname_nucleus))+.5)
    ax.set_yticklabels(cellname_nucleus[cellorder])
    ax.set_xticks(np.arange(0,len(header_now))+.5)
    ax.set_xticklabels(header_now)
    fig.colorbar(im, ax=ax)
    ax.tick_params(axis='x',labelsize = plot_parameters['projection_pattern']['label_size'])
    ax.tick_params(axis='y',labelsize = plot_parameters['projection_pattern']['label_size'])
    if plot_parameters['projection_pattern']['add_subtle_features']:
        plt.xticks(rotation = 70)
        ax.tick_params(axis='x',labelsize = plot_parameters['projection_pattern']['label_size'])


    child_nuclei_structureId = list()
    children_structureId = list()
    projection_pattern = list()
    header = list()
    for base_structureid in base_nuclei_structureId:
        children_idxs = np.where(allen_df['parentStructureId'] == base_structureid)[0]
        if len(children_idxs) == 0: # in case there are no children, we keep the parents
            children_idxs = np.where(allen_df['structureId'] == base_structureid)[0]
            print('end of tree reached at {}'.format(allen_df['name'][children_idxs[0]]))
        for children_idx in children_idxs:
            children_idx_contra = children_idx+len(allen_df)
            children_nucleus = allen_df['acronym'][children_idx]
            children_structureId = allen_df['structureId'][children_idx]
            child_nuclei_structureId.append(children_structureId)
            if plot_parameters['projection_pattern']['separate_ipsi_contra_projections']:
                header.append('{}-Y'.format(children_nucleus))
                header.append('{}-C'.format(children_nucleus))
                projection_pattern.append(projection_matrix_to_use[:,children_idx])
                projection_pattern.append(projection_matrix_to_use[:,children_idx_contra])
            else:
                header.append('{}'.format(children_nucleus))
                projection_pattern.append(projection_matrix_to_use[:,children_idx]+projection_matrix_to_use[:,children_idx_contra])
    tokeep = np.sum(np.asarray(projection_pattern),1)>0    
    projection_pattern_orig = np.asarray(projection_pattern).T
    projection_pattern = projection_pattern_orig[:,tokeep] 
    header=np.asarray(header)[tokeep]  
    #% plot
    if plot_parameters['projection_pattern']['add_subtle_features']:
        projection_pattern_now = np.concatenate([np.asarray(projection_pattern),data['subtle_features_matrix']],1)
        header_now = np.concatenate([header,data['subtle_features_header']])
    else:
        projection_pattern_now = projection_pattern
        header_now = header

    fig = plt.figure()
    ax=fig.add_subplot(111)
# =============================================================================
#     im = ax.pcolormesh(projection_pattern_now[cellorder,:],cmap = cmap,edgecolor=gridcolor,
#                        linewidth=gridwidth,
#                        vmin=vmin)
# =============================================================================
    im = ax.pcolormesh(projection_pattern_now[cellorder,:],
                       cmap = cmap,
                       edgecolor= plot_parameters['projection_pattern']['gridcolor'],
                       linewidth=plot_parameters['projection_pattern']['gridwidth'],
                       vmin=plot_parameters['projection_pattern']['vmin'])
    ax.set_yticks(np.arange(0,len(cellname_nucleus))+.5)
    ax.set_yticklabels(cellname_nucleus[cellorder])
    ax.set_xticks(np.arange(0,len(header_now))+.5)
    ax.set_xticklabels(header_now)
    fig.colorbar(im, ax=ax)
    ax.tick_params(axis='x',labelsize = plot_parameters['projection_pattern']['label_size'])
    ax.tick_params(axis='y',labelsize = plot_parameters['projection_pattern']['label_size'])
    if plot_parameters['projection_pattern']['add_subtle_features']:
        plt.xticks(rotation = 70)
        ax.tick_params(axis='x',labelsize = plot_parameters['projection_pattern']['label_size'])



    #% grandchildren
    grandchild_nuclei_structureId = list()
    projection_pattern = list()
    grandchildren_name_acronym = list()
    header = list()
    for children_structureid in child_nuclei_structureId:
        grandchildren_idxs = np.where(allen_df['parentStructureId'] == children_structureid)[0]
        if len(grandchildren_idxs) == 0: # in case there are no children, we keep the parents
            grandchildren_idxs = np.where(allen_df['structureId'] == children_structureid)[0]
            print('end of tree reached at {}'.format(allen_df['name'][grandchildren_idxs[0]]))
        if children_structureid in [856,864]: #these two structures (thalamus) are further divided
            #print(grandchildren_idxs)
            grandchildren_idxs_real = list()
            for grandchildren_idx in grandchildren_idxs:
                grandchildren_idxs_real.append(np.where(allen_df['parentStructureId'] == allen_df['structureId'][grandchildren_idx])[0])
            grandchildren_idxs = np.concatenate(grandchildren_idxs_real)
            #print(grandchildren_idxs)
        for grandchildren_idx in grandchildren_idxs:
            grandchildren_idx_contra = grandchildren_idx+len(allen_df)
            #grandchildren_nucleus = allen_df['name'][grandchildren_idx]
            grandchildren_acronym = allen_df['acronym'][grandchildren_idx]
            grandchildren_name_acronym = grandchildren_acronym#'{} - {}'.format(grandchildren_acronym,grandchildren_nucleus)
            grandchildren_structureId = allen_df['structureId'][grandchildren_idx]
            grandchild_nuclei_structureId.append(grandchildren_structureId)
            if plot_parameters['projection_pattern']['separate_ipsi_contra_projections']:
                header.append('{}-Y'.format(grandchildren_name_acronym))
                header.append('{}-C'.format(grandchildren_name_acronym))
                projection_pattern.append(projection_matrix_to_use[:,grandchildren_idx])
                projection_pattern.append(projection_matrix_to_use[:,grandchildren_idx_contra])
            else:
                header.append('{}'.format(grandchildren_name_acronym))
                projection_pattern.append(projection_matrix_to_use[:,grandchildren_idx]+projection_matrix_to_use[:,grandchildren_idx_contra])
    #%
    tokeep = np.sum(np.asarray(projection_pattern),1)>0    
    projection_pattern_orig = np.asarray(projection_pattern).T
    projection_pattern = projection_pattern_orig[:,tokeep] 
    header=np.asarray(header)[tokeep]   

    if plot_parameters['projection_pattern']['add_subtle_features']:
        projection_pattern_now = np.concatenate([np.asarray(projection_pattern),data['subtle_features_matrix']],1)
        header_now = np.concatenate([header,data['subtle_features_header']])
    else:
        projection_pattern_now = projection_pattern
        header_now = header

    fig = plt.figure()
    ax=fig.add_subplot(111)
    im = ax.pcolormesh(projection_pattern_now[cellorder,:],
                       cmap = cmap,
                       edgecolor= plot_parameters['projection_pattern']['gridcolor'],
                       linewidth=plot_parameters['projection_pattern']['gridwidth'],
                       vmin=plot_parameters['projection_pattern']['vmin'])
    ax.set_yticks(np.arange(0,len(cellname_nucleus))+.5)
    ax.set_yticklabels(cellname_nucleus[cellorder])
    plt.xticks(np.arange(0,len(header_now))+.5,header_now,rotation='vertical')
    #ax.set_xticklabels(header)
    fig.colorbar(im, ax=ax)   
    ax.tick_params(axis='x',labelsize = plot_parameters['projection_pattern']['label_size'])     
    ax.tick_params(axis='y',labelsize = plot_parameters['projection_pattern']['label_size'])
    plt.xticks(rotation = 70)
    
    
def plot_ground_truth_output_nuclei(data,allen_df):
    #%% plot ground truth output nuclei
    from rastermap import Rastermap
    order_nuclei_alphabetically = True #source nuclei
    order_with_rastermap = False
    use_acronyms = True
    use_full_nucleus_names = False
    use_full_nucleus_names_with_parents = False
    projection_matrix = data['target_locations_matrix']#/np.sum(target_locations_matrix,1)[:,None]

    if use_acronyms :
        header = data['target_locations_header']#_with_parent
        acronyms = list()
        for head in data['target_locations_header']:
            acronyms.append(allen_df.loc[allen_df['name'] == head,'acronym'].values[0])
        header = np.asarray(acronyms)
    elif use_full_nucleus_names:
        header = data['target_locations_header']
    else:
        header = data['target_locations_header_with_parent']
        header_new = list()
        for head in header:
            header_new.append(head.replace('-','\n'))    
        header = np.asarray(header_new)
        
    cellname_nucleus = list()
    soma_locations_acronyms = list()
    try:
        for cellname,nucleus in zip(data['cell_names'],data['soma_locations']):
            acronym = allen_df.loc[allen_df['name']==nucleus,'acronym'].values[0]
            #cellname_nucleus.append('{}        {}   '.format(cellname,acronym.rjust(10, ' ')))
            cellname_nucleus.append('{}   -   {}'.format(acronym,cellname))
            soma_locations_acronyms.append(acronym)
    except:
        cellname_nucleus = data['cell_names']
        soma_locations_acronyms = data['cell_names']
    if order_nuclei_alphabetically:
        cellorder = np.argsort(cellname_nucleus)
    else:
        cellorder =np.arange(len(cellname_nucleus))
    if order_with_rastermap:
        model = Rastermap(n_components=1, n_X=40, nPC=20, init='pca')
        model.fit(projection_matrix)
        cellorder = model.isort
    cellorder_manual = cellorder
    cellname_nucleus = np.asarray(cellname_nucleus)

    #% plot
    plt.style.use('dark_background')
    fig = plt.figure(figsize = [20,20])
    ax=fig.add_subplot(111)

    #im = ax.imshow(projection_matrix[cellorder,:],aspect='auto',cmap = 'RdPu')

    tokeep = np.sum(np.asarray(projection_matrix),0)>0    #removing target nuclei that don't receive input
    projection_matrix = projection_matrix[:,tokeep] 
    header=np.asarray(header)[tokeep]   

    im = ax.pcolormesh(projection_matrix[cellorder,:],cmap = 'tab20b',edgecolor='k', linewidth=1) #
    ax.set_yticks(np.arange(0,len(cellname_nucleus))+.5)
    ax.set_yticklabels(cellname_nucleus[cellorder])
    ax.set_xticks(np.arange(0,len(header))+.5)
    ax.set_xticklabels(header)
    #fig.colorbar(im, ax=ax)
    ax.tick_params(axis='x',labelsize = 14)
    ax.tick_params(axis='y',labelsize = 14)
    if not use_acronyms:
        plt.xticks(rotation = 'vertical')
    else:
        plt.xticks(rotation = 70)


#%%
def plot_axon_dendrit_statistics(original_data,cluster_dict):
    
   #%% 
    fig = plt.figure(figsize = [20,15])
    plt.style.use('default')

    ax_juci_axon_len = fig.add_subplot(3,4,1) #(2,2,1)
    ax_juci_axon_len.set_xlabel(r'Axon length ($\mu$m)')
    ax_juci_axon_len.set_ylabel('Proportion of cells')
    ax_juci_axon_end_point = fig.add_subplot(3,4,2)#(2,2,2)
    ax_juci_axon_end_point.set_xlabel('Axon end point number')
    ax_juci_dendrite_len = fig.add_subplot(3,4,3)#(2,2,3)
    ax_juci_dendrite_len.set_xlabel(r'Dendrite length ($\mu$m)')
    #ax_juci_dendrite_len.set_ylabel('Proportion of cells')
    ax_juci_dendrite_end_point = fig.add_subplot(3,4,4)#(2,2,4)
    ax_juci_dendrite_end_point.set_xlabel('Dendrite end point number')


    ax_juci_axon_len_vs_end_point = fig.add_subplot(3,4,5)
    ax_juci_axon_len_vs_end_point.set_xlabel('Axon end point number')
    ax_juci_axon_len_vs_end_point.set_ylabel('Axon length')
    
    ax_juci_axon_len_end_ratio = fig.add_subplot(3,4,6)
    ax_juci_axon_len_end_ratio.set_xlabel('Mean axon internode distance')
    
    

    ax_juci_dendrite_len_vs_end_point = fig.add_subplot(3,4,7)
    ax_juci_dendrite_len_vs_end_point.set_xlabel('Dendrite end point number')
    ax_juci_dendrite_len_vs_end_point.set_ylabel('Dendrite length')

    ax_juci_dendrite_len_end_ratio = fig.add_subplot(3,4,8)
    ax_juci_dendrite_len_end_ratio.set_xlabel('Mean dendrite internode distance')

    
    ax_dendrite_end_point_vs_internode = fig.add_subplot(3,4,9)
    ax_dendrite_end_point_vs_internode.set_xlabel('Dendrite bounding sphere radius (microns)')
    ax_dendrite_end_point_vs_internode.set_ylabel('Dendrite branch number')
    
    ax_dendrite_bounding_sphere = fig.add_subplot(3,4,10)
    ax_dendrite_bounding_sphere.set_xlabel('Dendrite bounding sphere radius (microns)')
    
    ax_dendrite_internode_number = fig.add_subplot(3,4,11)
    ax_dendrite_internode_number.set_xlabel('Dendrite internode number')    

    axon_internodes = np.asarray(original_data['axon_lengths'])/(np.asarray(original_data['axon_end_point_numbers'])+np.asarray(original_data['axon_branch_point_numbers']))
    dendrite_internodes = np.asarray(original_data['dendrite_lengths'])/(np.asarray(original_data['dendrite_end_point_numbers'])+np.asarray(original_data['dendrite_branch_point_numbers']))
    
    histrange_axon_end_point = np.arange(np.min(original_data['axon_end_point_numbers']),np.max(original_data['axon_end_point_numbers']),1)
    histrange_axon_len = np.arange(np.min(original_data['axon_lengths']),np.max(original_data['axon_lengths']),10)
    histrange_axon_len_end_ratio = np.arange(np.min(axon_internodes),np.max(axon_internodes),10)
    histrange_dendrite_len = np.arange(np.min(original_data['dendrite_lengths']),np.max(original_data['dendrite_lengths']),10)
    histrange_dendrite_endpoint = np.arange(np.min(original_data['dendrite_end_point_numbers']),np.max(original_data['dendrite_end_point_numbers']),1)
    histrange_dendrite_len_end_ratio = np.arange(np.min(dendrite_internodes),np.max(dendrite_internodes),10)
    histrange_dendrite_bounding_sphere_radius = np.arange(np.min(original_data['dendrite_bounding_sphere_radius']),np.max(original_data['dendrite_bounding_sphere_radius']),10)
    histrange_dendrite_internode_number = np.arange(np.min(np.asarray(original_data['dendrite_end_point_numbers'])+np.asarray(original_data['dendrite_branch_point_numbers'])),np.max(np.asarray(original_data['dendrite_end_point_numbers'])+np.asarray(original_data['dendrite_branch_point_numbers'])),10)
    for juci_cluster in cluster_dict['cluster_indices'].keys():
        idx = cluster_dict['cluster_indices'][juci_cluster]
        
        hist_y,hist_x = np.histogram(np.asarray(original_data['axon_end_point_numbers'])[idx],histrange_axon_end_point)
        hist_x = np.mean([hist_x[:-1],hist_x[1:]],0)
        ax_juci_axon_end_point.plot(hist_x,np.cumsum(hist_y)/sum(hist_y),
                                    color = cluster_dict['cluster_colors'][juci_cluster],
                                    label = cluster_dict['cluster_names'][juci_cluster],
                                    linewidth = 4)
        
        
        hist_y,hist_x = np.histogram(np.asarray(original_data['axon_lengths'])[idx],histrange_axon_len)
        hist_x = np.mean([hist_x[:-1],hist_x[1:]],0)
        ax_juci_axon_len.plot(hist_x,np.cumsum(hist_y)/sum(hist_y),
                              color = cluster_dict['cluster_colors'][juci_cluster],
                              label = cluster_dict['cluster_names'][juci_cluster],
                              linewidth = 4)
        
        hist_y,hist_x = np.histogram(axon_internodes[idx],histrange_axon_len_end_ratio)
        hist_x = np.mean([hist_x[:-1],hist_x[1:]],0)
        ax_juci_axon_len_end_ratio.plot(hist_x,np.cumsum(hist_y)/sum(hist_y),
                                        color = cluster_dict['cluster_colors'][juci_cluster],
                                        label = cluster_dict['cluster_names'][juci_cluster],
                                        linewidth = 4)
        
        
        hist_y,hist_x = np.histogram(np.asarray(original_data['dendrite_lengths'])[idx],histrange_dendrite_len)
        hist_x = np.mean([hist_x[:-1],hist_x[1:]],0)
        ax_juci_dendrite_len.plot(hist_x,np.cumsum(hist_y)/sum(hist_y),
                                  color = cluster_dict['cluster_colors'][juci_cluster],
                                  label = cluster_dict['cluster_names'][juci_cluster],
                                  linewidth = 4)
        
        hist_y,hist_x = np.histogram(np.asarray(original_data['dendrite_end_point_numbers'])[idx],histrange_dendrite_endpoint)
        hist_x = np.mean([hist_x[:-1],hist_x[1:]],0)
        ax_juci_dendrite_end_point.plot(hist_x,np.cumsum(hist_y)/sum(hist_y),
                                        color = cluster_dict['cluster_colors'][juci_cluster],
                                        label = cluster_dict['cluster_names'][juci_cluster],
                                        linewidth = 4)
        
        hist_y,hist_x = np.histogram(dendrite_internodes[idx],histrange_dendrite_len_end_ratio)
        hist_x = np.mean([hist_x[:-1],hist_x[1:]],0)
        ax_juci_dendrite_len_end_ratio.plot(hist_x,np.cumsum(hist_y)/sum(hist_y),
                                            color = cluster_dict['cluster_colors'][juci_cluster],
                                            label = cluster_dict['cluster_names'][juci_cluster],
                                            linewidth = 4)
        
        ax_juci_axon_len_vs_end_point.plot(np.asarray(original_data['axon_end_point_numbers'])[idx],np.asarray(original_data['axon_lengths'])[idx],'o',color = cluster_dict['cluster_colors'][juci_cluster],label = cluster_dict['cluster_names'][juci_cluster])
        
        ax_juci_dendrite_len_vs_end_point.plot(np.asarray(original_data['dendrite_end_point_numbers'])[idx],np.asarray(original_data['dendrite_lengths'])[idx],'o',color = cluster_dict['cluster_colors'][juci_cluster],label = cluster_dict['cluster_names'][juci_cluster])
        
        ax_dendrite_end_point_vs_internode.plot(np.asarray(original_data['dendrite_bounding_sphere_radius'])[idx],
                                                np.asarray(original_data['dendrite_branch_numbers'])[idx],
                                                'o',color = cluster_dict['cluster_colors'][juci_cluster],
                                                label = cluster_dict['cluster_names'][juci_cluster],
                                                alpha = .8)
        
        
        
        hist_y,hist_x = np.histogram(np.asarray(original_data['dendrite_bounding_sphere_radius'])[idx],histrange_dendrite_bounding_sphere_radius)
        hist_x = np.mean([hist_x[:-1],hist_x[1:]],0)
        ax_dendrite_bounding_sphere.plot(hist_x,np.cumsum(hist_y)/sum(hist_y),
                                            color = cluster_dict['cluster_colors'][juci_cluster],
                                            label = cluster_dict['cluster_names'][juci_cluster],
                                            linewidth = 4)
        
        
        hist_y,hist_x = np.histogram((np.asarray(original_data['dendrite_end_point_numbers'])+np.asarray(original_data['dendrite_branch_point_numbers']))[idx],histrange_dendrite_internode_number)
        hist_x = np.mean([hist_x[:-1],hist_x[1:]],0)
        ax_dendrite_internode_number.plot(hist_x,np.cumsum(hist_y)/sum(hist_y),
                                            color = cluster_dict['cluster_colors'][juci_cluster],
                                            label = cluster_dict['cluster_names'][juci_cluster],
                                            linewidth = 4)
        
    #ax_juci_dendrite_end_point.legend()
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
