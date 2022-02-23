from scyjava import jimport
import os
import numpy as np

Tree = jimport('sc.fiji.snt.Tree')
TreeStatistics = jimport('sc.fiji.snt.analysis.TreeStatistics')





AllenUtils = jimport('sc.fiji.snt.annotation.AllenUtils')
AllenCompartment = jimport('sc.fiji.snt.annotation.AllenCompartment')
sntrgbcolor = jimport('org.scijava.util.ColorRGB')

def init_viewer():
    Viewer3D = jimport('sc.fiji.snt.viewer.Viewer3D')
    sntcolor = jimport('org.scijava.util.ColorRGB')
    viewer_3d = Viewer3D(True)#
    viewer_3d.setEnableDarkMode(True)
    brain = viewer_3d.loadRefBrain('mouse')
    #viewer_3d.setViewMode('sagittal')
    viewer_3d.setAnimationEnabled(True)
    viewer_3d.setDefaultThickness(1.5)
    viewer_3d.show()
    cells_added = {}  
    return viewer_3d, cells_added


def add_cell_to_viewer(viewer_3d,
                       jsonfolder,
                       cell_name = 'AA0922',
                       cells_added = {},
                       axon_color= 'red',
                       dendrite_color = 'blue',
                       show_axon = True,
                       show_dendrite = True,
                       all_cells_on_left = True,):
    #%
    if cell_name in cells_added.keys():
        viewer_3d.remove(cells_added[cell_name][0])
        viewer_3d.remove(cells_added[cell_name][1])
        del cells_added[cell_name]
    else:
        #Tree = autoclass('sc.fiji.snt.Tree')  # ezzel nyitod meg a rekonstrukciot
        Tree = jimport('sc.fiji.snt.Tree')  # ezzel nyitod meg a rekonstrukciot
        
        jsonfile_with_path = os.path.join(jsonfolder,cell_name+'.json')
        if show_axon:
            axon = Tree(jsonfile_with_path, "axon")
            axon.setColor(axon_color)
            axon.setLabel(axon.getLabel()+' axon')
            axon.setRadii(8)
        else:
            axon=None
        dendrite = Tree(jsonfile_with_path, "dendrite")
        dendrite.setLabel(dendrite.getLabel()+' dendrite')
        dendrite.setColor(dendrite_color)
        
        if not AllenUtils.isLeftHemisphere(dendrite) and all_cells_on_left:
                AllenUtils.assignToLeftHemisphere(axon)
                AllenUtils.assignToLeftHemisphere(dendrite)
        if show_axon:
            viewer_3d.add(axon)
        if show_dendrite:
            viewer_3d.add(dendrite)
        cells_added[cell_name] = [axon,dendrite]
    return cells_added#{cell_name:[axon,dendrite]}
# =============================================================================
# def plot_points_mean_std(ax_axon_length,axon_length_now,cell_names_now,color = 'red'):
#     ax_axon_length.clear()
#     ax_axon_length.plot(np.ones(len(axon_length_now)),axon_length_now,'o',color = color, markersize = 12)    
#     if len(axon_length_now)>1:
#         ax_axon_length.errorbar(1,np.mean(axon_length_now),np.std(axon_length_now),ecolor = 'black',elinewidth=3)
#         ax_axon_length.bar(1,np.mean(axon_length_now),color = 'w',edgecolor = 'k',linewidth = 3)
#     for name,y in zip(cell_names_now,axon_length_now):
#         ax_axon_length.text(1.2, y, name, fontsize=8)
# def add_cell_to_descriptive_plot(cells_added):
#     cell_idxs = list()
#     for cell_id in cells_added.keys():
#         cell_idxs.append(np.argmax(np.asarray(cell_names)==cell_id))
#     axon_length_now = np.asarray(axon_lengths)[cell_idxs]
#     plot_points_mean_std(ax_axon_length,axon_length_now,cells_added.keys(),color = 'red')
#     ax_axon_length.set_title('Axon length')
#     ax_axon_length.set_ylabel('microns')
#     
#     axon_branch_point_numbers_now = np.asarray(axon_branch_point_numbers)[cell_idxs]
#     plot_points_mean_std(ax_axon_branch_points,axon_branch_point_numbers_now,cells_added.keys(),color = 'red')
#     ax_axon_branch_points.set_title('Axon branch point number')
#     ax_axon_branch_points.set_ylabel('#')
#     
#     axon_end_point_numbers_now = np.asarray(axon_end_point_numbers)[cell_idxs]
#     plot_points_mean_std(ax_axon_end_points,axon_end_point_numbers_now,cells_added.keys(),color = 'red')
#     ax_axon_end_points.set_title('Axon end point number')
#     ax_axon_end_points.set_ylabel('#')
#         
#     dendrite_length_now = np.asarray(dendrite_lengths)[cell_idxs]
#     plot_points_mean_std(ax_dendrite_length,dendrite_length_now,cells_added.keys(),color = 'blue')
#     ax_dendrite_length.set_title('Dendrite length')
#     ax_dendrite_length.set_ylabel('microns')
#     
#     dendrite_branch_point_numbers_now = np.asarray(dendrite_branch_point_numbers)[cell_idxs]
#     plot_points_mean_std(ax_dendrite_branch_points,dendrite_branch_point_numbers_now,cells_added.keys(),color = 'blue')
#     ax_dendrite_branch_points.set_title('Dendrite branch point number')
#     ax_dendrite_branch_points.set_ylabel('#')
#     
#     dendrite_end_point_numbers_now = np.asarray(dendrite_end_point_numbers)[cell_idxs]
#     plot_points_mean_std(ax_dendrite_end_points,dendrite_end_point_numbers_now,cells_added.keys(),color = 'blue')
#     ax_dendrite_end_points.set_title('Dendrite end point number')
#     ax_dendrite_end_points.set_ylabel('#')
#     
#     
#     fig_cell_description.canvas.draw()
# =============================================================================
    
def add_brain_region_to_viewer(structureId,color,transparency,viewer_3d):
    allenmesh = AllenCompartment(structureId).getMesh()
    allenmesh.setColor(sntrgbcolor(color[0],color[1],color[2]),transparency) # color, transparency
    viewer_3d.addMesh(allenmesh)
    