# -*- coding: utf-8 -*-
"""
@author: Alexandre JANIN
@aim:    Graphical routines
"""

# External dependencies:
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
import h5py
from pathlib import Path

# Internal dependencies
from .generics import im
from .geotransform import displacement_sdt_to_xyz, normal2fault


# ----------------- FUNCTIONS -----------------


def plotFault2D(patches, patchIDs=None, group=None,\
                v=None, normal=True, refID=0,\
                view_vec=None, cmap=plt.cm.viridis,\
                vmin=None, vmax=None,\
                patches2hatch=None,hatchColor='white',\
                center=True,center_size=30,center_color='white',\
                vec=None, vec_scale=1e2, vec_width=4e-3, vec_color='k',\
                vec_decim=1, edgecolor="k", linewidth=0.2,\
                figsize=(10, 4), cbar_title=None,\
                aspect='equal',invert_xaxis=False, invert_yaxis=False,\
                title='Fault patches (projected view)',\
                xlabel='Projected X', ylabel='Projected Y',
                show=True,savefig=False,fname='myFigure.png',path='./',\
                dpi=300,verbose=True):
    """
    Fast 2D plot of fault patches projected onto a plane
    orthogonal to a given 3D vector.

    Args:
        patches (PatchCollection object):
                Fault patches
        patchIDs (None or numpy.ndarray, dtype=int, optional):
                List of patch IDs that will be display.
                The others will be ignored.
                (Default, None)
        group (None, int, list, optional):
                Display only the patches of the asked group:
                If None then, display all the group (default)
                If int displays only the patches of this group
                If list displays the pathces of all the groups listed
        v (array-like or None, optional):
                Values to color patches (e.g. Displacement)
                (Defaults: None)
        normal (bool, optional):
                Condition to plot the collection of patches
                in a plane normal to the normal to the patch number 'refID'.
                Normal computed with the function fault2normal.
                (Defaults: True)
        refID (int, optional):
                ID of the reference patch when normal is set to True.
                (Defaults: refID = 0)
        view_vec (array-like, shape (3,) or None, optional):
                3D vector defining projection direction.
                If not None, take the priority over the argument 'normal'
                (Default: None)
        cmap (matplotlib.colors.ListedColormap, optional):
                Matplotlib colormap when v is not None
                (Default: matplotlib.pyplot.cm.viridis)
        vmin (scalar, optional):
                minimum value for the color map. When None, set automatically
                (Default=None)
        vmax (scalar, optional):
                maximum value for the color map. When None, set automatically
                (Default=None)
        vec (py3Ddef.ElementStrDispl or None, optional):
                Vectorial field to be display as a quiver.
                The shape as to be the same as the geometry of the input PatchCollection
                (Default: None)
        vec_scale (scalar, optional):
                Scale of the quiver when vec is not None. (Default: 1e2)
        vec_width (scalar, optional):
                Width of the quiver when vec is not None. (Default: 2e-3)
        vec_color (str, optional):
                Color of the quiver when vec is not None. (Default: 'k')
        vec_decim (int, optional):
                Decimation factor to display the vectorial field.
                Display 1 every 'vec_decim' vectors.
                (Default: vec_decim=1, i.e. no decimation)
        edgecolor (str, optional):
                Patch edge color (Default: 'k')
        linewidth (float, optional):
                Patch edge linewidth. (Default, 0.2)
        figsize (tuple, optional):
                Shape of the figure. (Default: (7, 5))
        cbar_title (str, optional):
                Title of the colorbar
        aspect (str, optional):
                Aspect of the axis.
                See options of matplotlib.pyplot axis.
                (Default: 'equal')
        invert_xaxis (bool, optional):
                Option controlling the reversal of the x axis.
                (Default: False)
        invert_yaxis (bool, optional):
                Option controlling the reversal of the y axis.
                (Default: False)
        title (str, optional):
                Title of the figure
                (Default, 'Fault patches (projected view)')
        xlabel (str, optional):
                Label of the x axis.
                (Default: 'Projected X')
        ylabel (str, optional):
                Label of the y axis.
                (Default: 'Projected Y')
        verbose (bool, optional):
                Option controlling the verbose output of the function
    """
    pName = 'plotFault2D'
    im('Projected 2D figure of the patch collection', pName, verbose)

    fig = plt.figure(figsize=figsize)
    ax  = fig.add_subplot(111)

    if normal and view_vec is None:
        strikeRef  = patches.strike[refID]
        dipRef     = patches.dip[refID]
        view_vec = normal2fault(strikeRef, dipRef)
        if xlabel == 'Projected X':
            xlabel = 'Along strike'
        if ylabel == 'Projected Y':
            ylabel = 'Depth'
            invert_yaxis = True
    else:
        normal = False
        view_vec = np.asarray(view_vec, dtype=float)
    
    im('Projection vector:', pName, verbose)
    im('    - reference patch: '+str(refID), pName, verbose)
    im('    - normal to the reference patch: '+str(normal), pName, verbose)
    im('    - view vector: '+str(view_vec), pName, verbose)

    im('Patches displayed:', pName, verbose)
    if patchIDs is None:
        patchIDs = patches.ids.copy()
        im('    - Patch IDs: all', pName, verbose)
    else:
        if isinstance(patchIDs, list) or isinstance(patchIDs, tuple) or isinstance(patchIDs, np.ndarray):
            im('    - Patch IDs: '+str(patchIDs), pName, verbose)
        else:
            raise TypeError('PatchIDs needs to be a list or a tuple or a numpy.ndarray.')
    if group is None:
        group = np.unique(patches.group)
        im('    - Patch group: all', pName, verbose)
    else:
        if np.isscalar(group): # transform to list
            group = [group]
        if isinstance(group, list) or isinstance(group, tuple) or isinstance(group, np.ndarray):
            im('    - group: '+str(group), pName, verbose)
        else:
            raise TypeError('group needs to be a list or a tuple or a numpy.ndarray.')

    view_vec /= np.linalg.norm(view_vec)

    # Build orthonormal basis (e1, e2) for projection plane
    tmp = np.array([0.0, 0.0, 1.0])
    if np.allclose(np.abs(np.dot(tmp, view_vec)), 1.0):
        tmp = np.array([0.0, 1.0, 0.0])

    e1 = np.cross(view_vec, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(view_vec, e1)

    polygons = []
    polygons_hatch = []

    if vec is not None:

        if vec_decim < 1 or not isinstance(vec_decim, int):
            raise TypeError('vec_decim must be an integer >=1')
        mask_vec = np.zeros(vec.nop, dtype=bool)
        mask_vec[::vec_decim] = True

        X, Y, U, V = [], [], [], []

    j = 0
    v2plot = []
    v2plot_hatch = []
    any_hatches = False

    xc2display = []
    yc2display = []

    im('Iterate on patches', pName, verbose)
    for el in patches.patches:

        if patches.ids[j] in patchIDs and patches.group[j] in group:

            strike = np.deg2rad(el.strike)
            dip = np.deg2rad(el.dip)

            # Unit vectors
            s = np.array([np.sin(strike), np.cos(strike), 0.0])  # strike
            d = np.array([
             np.cos(strike) * np.cos(dip),
            -np.sin(strike) * np.cos(dip),
            -np.sin(dip)
            ])  # dip (downward)

            # Patch corners in 3D
            p0 = np.array([el.x0, el.y0, el.z0])
            p1 = p0 + el.L * s
            p2 = p1 + el.W * d
            p3 = p0 + el.W * d

            patch3d = [p0, p1, p2, p3]

            # Project onto (e1, e2)
            patch2d = [
                [np.dot(p, e1), np.dot(p, e2)]
                for p in patch3d
            ]

            if patches2hatch is None:
                polygons.append(patch2d)
                v2plot.append(v[j])
            else:
                polygons.append(patch2d)
                v2plot.append(v[j])
                if patches2hatch[j]:
                    any_hatches = True
                    polygons_hatch.append(patch2d)
                    v2plot_hatch.append(v[j])
            
            if center:
                p = np.array(el.center)
                xc2display.append(np.dot(p, e1))
                yc2display.append(np.dot(p, e2))

            if vec is not None and mask_vec[j]:

                element_j = patches.patches[j]
                usj = vec.us[j]
                udj = vec.ud[j]
                utj = vec.un[j]

                # Element center
                p = np.array(el.center)

                # Convert displacement to xyz
                u_xyz = displacement_sdt_to_xyz(usj, udj, utj, element_j.strike, element_j.dip)

                # Project position
                X.append(np.dot(p, e1))
                Y.append(np.dot(p, e2))

                # Project vector
                U.append(np.dot(u_xyz, e1))
                V.append(np.dot(u_xyz, e2))

        j += 1

    # display the reprojected patches
    polycoll = PolyCollection(
        polygons,
        array=v2plot,
        cmap=cmap,
        edgecolors=edgecolor,
        linewidths=linewidth
    )
    polycoll.set_clim([vmin, vmax])

    im('Plot the patch collection', pName, verbose)
    ax.add_collection(polycoll)

    # if any hatched polygon, displayed on the top
    if any_hatches:
        polycoll_hatch = PolyCollection(
            polygons_hatch,
            array=v2plot_hatch,
            cmap=cmap,
            edgecolors=hatchColor,
            linewidths=linewidth,
            hatch='///'
        )
        vmin, vmax = polycoll.get_clim() # set the same vmin, vmax
        polycoll_hatch.set_clim([vmin, vmax])
        ax.add_collection(polycoll_hatch)

    # plot quiver
    if vec is not None:
        im('Plot the displacement field', pName, verbose)
        ax.quiver(X, Y, U, V,
                scale=vec_scale,
                width=vec_width,
                color=vec_color,
                zorder=1)
    
    # display the centers
    if center:
        ax.scatter(xc2display, yc2display, s=center_size, facecolors='None', edgecolors=center_color)

    # Autoscale
    all_xy = np.vstack([np.array(p) for poly in polygons for p in poly])
    ax.set_xlim(all_xy[:, 0].min(), all_xy[:, 0].max())
    ax.set_ylim(all_xy[:, 1].min(), all_xy[:, 1].max())

    ax.set_aspect(aspect)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if invert_xaxis:
        ax.invert_xaxis()
    if invert_yaxis:
        ax.invert_yaxis()

    if v is not None:
        plt.colorbar(polycoll, ax=ax, label=cbar_title, shrink=0.5)

    fig.tight_layout()

    im('Figure: Done.', pName, verbose)

    if savefig:
        if path[-1] != '/':
            path += '/'
        fig.savefig(path+fname,dpi=dpi)

    if show:
        plt.show()
    else:
        plt.close()



def plotFault3D(patches, centers=True, refpts= True, displayGroup=True,
                cmap=plt.cm.viridis,figsize=(9,7),
                camera=(30,-45,0), aspect='equal'):
    """
    Interactive display of a set of dislocation in 3D.
    
    Args:
        patches (py3Ddef.geometry.PatchCollection):
                Patches that will be display in 3D
        centers (bool, optional):
                Option controlling the display of the centers of each patch.
                (Default, True)
        refpts (bool, optional):
                Option controlling the display of the reference points of each patch.
                (Default, True)
        displayGroup (bool, optional):
                Option controlling the display of each patch according to their group ID
                Will color patches according to the value of their group.
                (Default: True)
        cmap (matplotlib.colors.ListedColormap, optional):
                Matplotlib colormap used for the group when required.
                (Default matplotlib.pyplot.cm.viridis)
        figsize (tuple, optional):
                Shape of the figure. (Default: (9, 7))
        camera (list of length 3, optional):
                Parameter for the default camara orientation in 3D
                camera[0]: elevation
                camera[1]: azimuth
                camera[2]: roll
                (Default: (30, -45, 0))
        aspect (str, optional):
                Aspect of the axis.
                See options of matplotlib.pyplot axis.
                (Default: 'equal')
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')

    if displayGroup:
        if np.unique(patches.group).shape[0] > 1:
            norm = Normalize(vmin=patches.group.min(), vmax=patches.group.max())
        else:
            mygroup = np.unique(patches.group)
            norm = Normalize(vmin=mygroup-0.1, vmax=mygroup+0.1)

    i = 0
    for patch in patches.patches:
        x0, y0, z0 = patch.x0, patch.y0, patch.z0

        # Compute four corners of the patch
        corner1, corner2, corner3, corner4 = patch.corners

        # Create polygon
        poly = Poly3DCollection([[corner1, corner2, corner3, corner4]])
        if not displayGroup:
            color = 'lightblue'
        else:
            color = cmap(norm(patches.group[i]))
        poly.set_facecolor(color)
        poly.set_alpha(0.7)
        poly.set_edgecolor('k')
        ax.add_collection3d(poly)

        if refpts:
            ax.scatter(x0, y0, z0, color='red')
            
        if centers:
            xc, yc, zc = patch.center
            ax.scatter(xc, yc, zc, color='blue')
        
        i += 1
    
    # Set labels
    ax.set_xlabel('X, East (km)')
    ax.set_ylabel('Y, North (km)')
    ax.set_zlabel('Z, Depth (km)')
    ax.set_title('3D Fault Visualization')

    ax.set_aspect(aspect)

    # colorbar
    if displayGroup:
        mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        mappable.set_array(patches.group)
        fig.colorbar(mappable, ax=ax, label="Group", shrink=0.3)

    cam_elevation = camera[0]
    cam_azimuth   = camera[1]
    cam_roll      = camera[2]
    ax.view_init(elev=cam_elevation, azim=cam_azimuth, roll=cam_roll, vertical_axis='z')
    
    fig.tight_layout()
    plt.show()



def patches2paraview(fname, patches, fields=None, fieldnames=None, path='./'):
    """
    Write fault patches to XDMF + HDF5 for ParaView visualization.

    Adds automatically three vector fields at patch centers:
        - dir_strike
        - dir_dip
        - dir_normal

    Structure:
        Grid 1 : FaultPatches (quad mesh)
            - geometry
            - user scalar/tensor fields (cell data)

        Grid 2 : PatchCenters (point cloud)
            - centers of patches
            - dir_strike, dir_dip, dir_normal vectors (node data)

    In Paraview:
        - Open the XDMF file with the "XDMF reader" option (not the "Xdmf3 Reader S" or "Xdmf3 Reader T")
        - Separate "FaultPatches" and "PatchCenters" with "Extract Block" (avoid issues with "Glyph" representation)

        
    Args:
        patches (PatchCollection object): Dislocation patches
        fields (list of arrays, optional):
            Each array must be shape:
                (N,)    -> scalar field
                (N,3)   -> vector field
                (N,6)   -> symmetric tensor
                (N,9)   -> full tensor
                (N,k)   -> generic attribute
        fieldnames (list of str, optional):
            Names associated to 'fields' if fields is a list.
            Must have same length as fields.
        fname : str
            Base filename (without extension)
        path : str
            Output directory
    """
    if path[-1] != '/':
        path += '/'

    outdir = Path(path)
    outdir.mkdir(parents=True, exist_ok=True)

    h5file = outdir / f"{fname}.h5"
    xdmf_file = outdir / f"{fname}.xdmf"

    N = patches.nop

    # ============================================================
    # 1) BUILD PATCH GEOMETRY
    # ============================================================

    points = []
    cells = []

    centers = np.zeros((N, 3))
    strike_dirs = np.zeros((N, 3))
    dip_dirs = np.zeros((N, 3))
    normal_dirs = np.zeros((N, 3))

    for i, el in enumerate(patches.patches):

        strike = np.deg2rad(el.strike)
        dip = np.deg2rad(el.dip)

        # local orthonormal frame
        s = np.array([np.sin(strike), np.cos(strike), 0.0])
        d = np.array([
            np.cos(strike) * np.cos(dip),
           -np.sin(strike) * np.cos(dip),
           -np.sin(dip)
        ])
        n = np.cross(s, d)

        s /= np.linalg.norm(s)
        d /= np.linalg.norm(d)
        n /= np.linalg.norm(n)

        strike_dirs[i] = s
        dip_dirs[i] = d
        normal_dirs[i] = n

        centers[i] = [el.xc, el.yc, el.zc]

        # quad corners
        p0 = np.array([el.x0, el.y0, el.z0])
        p1 = p0 + el.L * s
        p2 = p1 + el.W * d
        p3 = p0 + el.W * d

        base = len(points)
        points.extend([p0, p1, p2, p3])
        cells.append([base, base + 1, base + 2, base + 3])

    points = np.asarray(points)
    cells = np.asarray(cells, dtype=np.int32)

    # ---- Normalize fields input ----
    field_dict = {}

    # Empty input field
    if len(fields) == 1 and fields[0] is None:
        fields = [np.ones(patches.shape, dtype=np.int32)]
        fieldnames = ['empty']

    # constrain the list type
    fields = list(fields)
    fieldnames = list(fieldnames)

    if fieldnames is not None:
        if len(fieldnames) != len(fields):
            raise ValueError("fieldnames must have same length as fields")
        for name, arr in zip(fieldnames, fields):
            field_dict[name] = arr
    else:
        for i, arr in enumerate(fields):
            field_dict[f"field_{i}"] = arr

    # --- Helper
    def infer_attribute_type(arr):
        arr = np.asarray(arr)
        if arr.ndim == 1:
            return "Scalar", (arr.shape[0],)
        elif arr.ndim == 2:
            k = arr.shape[1]
            if k == 3:
                return "Vector", arr.shape
            elif k == 6:
                return "Tensor6", arr.shape
            elif k == 9:
                return "Tensor", arr.shape
            else:
                return "None", arr.shape
        else:
            raise ValueError(f"Unsupported field shape {arr.shape}")

    # ============================================================
    # 2) WRITE HDF5
    # ============================================================

    with h5py.File(h5file, "w") as f:

        # mesh
        f.create_dataset("Geometry", data=points)
        f.create_dataset("Connectivity", data=cells)

        # centers
        f.create_dataset("Centers", data=centers)
        f.create_dataset("Polyvertex", data=np.arange(N, dtype=np.int32))

        # directions
        f.create_dataset("dir_strike", data=strike_dirs)
        f.create_dataset("dir_dip", data=dip_dirs)
        f.create_dataset("dir_normal", data=normal_dirs)

        # user fields
        for name, arr in field_dict.items():
            arr = np.asarray(arr)
            if arr.shape[0] != N:
                raise ValueError(f"Field '{name}' size mismatch")
            f.create_dataset(name, data=arr)

    # ============================================================
    # 3) WRITE XDMF
    # ============================================================

    with open(xdmf_file, "w") as f:

        f.write(f"""<?xml version="1.0" ?>
<Xdmf Version="3.0">
  <Domain>

    <!-- ======================= -->
    <!-- Fault patches surface  -->
    <!-- ======================= -->
    <Grid Name="FaultPatches" GridType="Uniform">

      <Topology TopologyType="Quadrilateral" NumberOfElements="{N}">
        <DataItem Format="HDF" Dimensions="{N} 4">
          {h5file.name}:/Connectivity
        </DataItem>
      </Topology>

      <Geometry GeometryType="XYZ">
        <DataItem Format="HDF" Dimensions="{points.shape[0]} 3">
          {h5file.name}:/Geometry
        </DataItem>
      </Geometry>
""")

        for name, arr in field_dict.items():
            attr_type, dims = infer_attribute_type(arr)
            dims_str = " ".join(map(str, dims))

            f.write(f"""
      <Attribute Name="{name}" AttributeType="{attr_type}" Center="Cell">
        <DataItem Format="HDF" Dimensions="{dims_str}">
          {h5file.name}:/{name}
        </DataItem>
      </Attribute>
""")

        f.write(f"""
    </Grid>

    <!-- ======================= -->
    <!-- Patch centers (vectors) -->
    <!-- ======================= -->
    <Grid Name="PatchCenters" GridType="Uniform">

      <Topology TopologyType="Polyvertex" NumberOfElements="{N}">
        <DataItem Format="HDF" Dimensions="{N}">
          {h5file.name}:/Polyvertex
        </DataItem>
      </Topology>

      <Geometry GeometryType="XYZ">
        <DataItem Format="HDF" Dimensions="{N} 3">
          {h5file.name}:/Centers
        </DataItem>
      </Geometry>

      <Attribute Name="dir_strike" AttributeType="Vector" Center="Node">
        <DataItem Format="HDF" Dimensions="{N} 3">
          {h5file.name}:/dir_strike
        </DataItem>
      </Attribute>

      <Attribute Name="dir_dip" AttributeType="Vector" Center="Node">
        <DataItem Format="HDF" Dimensions="{N} 3">
          {h5file.name}:/dir_dip
        </DataItem>
      </Attribute>

      <Attribute Name="dir_normal" AttributeType="Vector" Center="Node">
        <DataItem Format="HDF" Dimensions="{N} 3">
          {h5file.name}:/dir_normal
        </DataItem>
      </Attribute>

    </Grid>

  </Domain>
</Xdmf>
""")


def grid2paraview(fname, grid, fields=[None], fieldnames=None, path='./', precision='Float'):
    """
    precision : str
        "Float" or "Double"
    """

    if path[-1] != '/':
        path += '/'
    
    path2file = path + fname

    path2file = Path(path2file)
    h5name = path2file.with_suffix(".h5").name
    xdmfname = path2file.with_suffix(".xdmf")

    nx, ny, nz = grid.nx, grid.ny, grid.nz
    X, Y, Z = grid.xaxis, grid.yaxis, grid.zaxis

    # Empty input field
    if len(fields) == 1 and fields[0] is None:
        fields = [np.ones(grid.shape, dtype=np.int32)]
        fieldnames = ['empty']

    # -------------------------
    # Write HDF5
    # -------------------------
    with h5py.File(path2file.with_suffix(".h5"), "w") as h5:
        h5.create_dataset("X", data=X)
        h5.create_dataset("Y", data=Y)
        h5.create_dataset("Z", data=Z)

        for i, field in enumerate(fields):
            name = (
                fieldnames[i]
                if fieldnames is not None
                else f"field_{i}"
            )
            h5.create_dataset(name, data=field)

    # -------------------------
    # Write XDMF (manual XML)
    # -------------------------
    dtype = "Float" if precision == "Float" else "Double"

    with open(xdmfname, "w") as f:
        f.write('<?xml version="1.0" ?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        f.write('<Xdmf Version="3.0">\n')
        f.write(' <Domain>\n')
        f.write('  <Grid Name="UniformGrid" GridType="Uniform">\n')

        # Topology
        f.write(
            f'   <Topology TopologyType="3DRectMesh" '
            f'Dimensions="{nz} {ny} {nx}"/>\n'
        )

        # Geometry
        f.write('   <Geometry GeometryType="VXVYVZ">\n')
        for coord, name in zip(["X", "Y", "Z"], [nx, ny, nz]):
            f.write(
                f'    <DataItem Dimensions="{name}" '
                f'NumberType="{dtype}" Precision="8" Format="HDF">\n'
            )
            f.write(f'     {h5name}:/{coord}\n')
            f.write('    </DataItem>\n')
        f.write('   </Geometry>\n')

        # Attributes
        for i, field in enumerate(fields):
            name = (
                fieldnames[i]
                if fieldnames is not None
                else f"field_{i}"
            )

            if field.ndim == 3:
                atype = "Scalar"
                dims = f"{nz} {ny} {nx}"
            elif field.ndim == 4 and field.shape[-1] == 3:  # displacement vector
                atype = "Vector"
                dims = f"{nz} {ny} {nx} 3"
            elif field.ndim == 4 and field.shape[-1] == 4:  # stress/strain invariants
                atype = "Vector"
                dims = f"{nz} {ny} {nx} 4"
            elif field.ndim == 4 and field.shape[-1] == 6:  # optimal failure place
                atype = "Vector"
                dims = f"{nz} {ny} {nx} 6"
            elif field.ndim == 4 and field.shape[-1] == 9:
                atype = "Tensor"
                dims = f"{nz} {ny} {nx} 9"
            else:
                raise ValueError(
                    f"Unsupported field shape for {name}: {field.shape}"
                )

            f.write(
                f'   <Attribute Name="{name}" '
                f'AttributeType="{atype}" Center="Node">\n'
            )
            f.write(
                f'    <DataItem Dimensions="{dims}" '
                f'NumberType="{dtype}" Precision="8" Format="HDF">\n'
            )
            f.write(f'     {h5name}:/{name}\n')
            f.write('    </DataItem>\n')
            f.write('   </Attribute>\n')

        f.write('  </Grid>\n')
        f.write(' </Domain>\n')
        f.write('</Xdmf>\n')
