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
from copy import deepcopy

# Internal dependencies
from .generics import im
from .geotransform import displacement_sdt_to_xyz, normal2fault
from . import GridScalar, GridVector, GridTensor, \
                    GridDisplacement, GridStress, GridStrain, \
                    GridPrincipalStrainOrientation, \
                    GridOptimalFailurePlane, GridDisplacementGradient


# ----------------- FUNCTIONS -----------------


def plotFault2D(patches, patchIDs=None, group=None,\
                v=None, normal=True, refID=0,\
                view_vec=None, cmap="viridis",\
                vmin=None, vmax=None,\
                vec=None, vec_scale=1e2, vec_width=2e-3, vec_color='k',\
                edgecolor="k", linewidth=0.2,\
                figsize=(10, 4), cbar_title=None,\
                aspect='equal',invert_xaxis=False, invert_yaxis=False,\
                title='Fault patches (projected view)',\
                xlabel='Projected X', ylabel='Projected Y',
                verbose=True):
    """
    Fast 2D plot of fault patches projected onto a plane
    orthogonal to a given 3D vector.

    Args:
        patches (PatchCollection object):
            Fault patches
        v (array-like or None):
            Values to color patches
        normal (bool):
            Condition plot the collection of patches
            in a plane normal of the normal to the 1st patch.
            Normal computed with the function fault2normal.
        view_vec (array-like, shape (3,) or None):
            3D vector defining projection direction.
            If not None, take the priority over the argument 'normal'
        cmap (str):
            Matplotlib colormap
        edgecolor (str):
            Patch edge color
        linewidth (float):
            Patch edge linewidth
        figsize (tuple):
            Shape of the figure. Default: (7, 5)
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

    if vec is not None:

        X, Y, U, V = [], [], [], []

    j = 0
    v2plot = []
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

            polygons.append(patch2d)

            v2plot.append(v[j])

            if vec is not None:

                element_j = patches.patches[j]
                usj = vec.ss[j]
                udj = vec.ds[j]
                utj = vec.ts[j]

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

    im('Plot the displacement field', pName, verbose)
    ax.quiver(X, Y, U, V,
            scale=vec_scale,
            width=vec_width,
            color=vec_color,
            zorder=1)

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

    plt.show()



def plotFault3D(patches, centers=True, refpts= True, displayGroup=True, cmap=plt.cm.viridis, figsize=(9,7), camera=(30,-45,0), aspect='equal'):
    """
    Plot 3D fault subpatches with arbitrary strike and dip.
    
    Parameters:
        patches : list of dicts
            Output from discretize_fault function
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
    Write XDMF + HDF5 files for ParaView visualization.

    Parameters
    ----------
    filename : str
        Base filename (without extension)
    precision : str
        "Float" or "Double"
    """

    if path[-1] != '/':
        path += '/'
    
    path2file = path + fname

    path2file = Path(path2file)
    h5name = path2file.with_suffix(".h5").name
    xdmfname = path2file.with_suffix(".xdmf")

    shape = grid.shape
    nx, ny, nz = grid.nx, grid.ny, grid.nz
    X, Y, Z = grid.xaxis, grid.yaxis, grid.zaxis

    # empty fields
    if len(fields) == 1 and fields[0] is None:
        fields = [np.zeros(int(grid.nx * grid.ny * grid.nz), dtype=np.int32).reshape(grid.shape)]
        fieldnames = ['empty']
    
    # deepcopy of the fields
    fields = deepcopy(fields)

    # -------------------------
    # Restructure fields
    # -------------------------

    # The idea is to homogeneise the type: only GridScalar, GridVector, GridTensor
    for i in range(len(fields)):
        field = fields[i]
        old_shape = field.shape
        # Vectorial fields
        if isinstance(field, GridDisplacement):
            field = GridVector(field.solution)
        # Tensor fields
        elif isinstance(field, GridStress) or isinstance(field, GridStrain) or \
             isinstance(field, GridPrincipalStrainOrientation) or \
             isinstance(field, GridOptimalFailurePlane) or \
             isinstance(field, GridDisplacementGradient):
            field = GridTensor(field.solution)
        # Scalar fields
        elif isinstance(field, GridScalar):
            pass
        else:
            raise TypeError('Unknown input type for the field number %s'%str(i))
        # restore the original shape
        field.reshape(old_shape)
        fields[i] = field # update
    
    for i in range(len(fields)):
        field = fields[i]
        # check shape and return array
        if field.shape == shape:
            myarray = field.array
            if len(myarray.shape) == 5: # tensor
                myarray = myarray.reshape(shape + (9,))
            fields[i] = myarray
        else:
            raise ValueError('The input field %s does not have a shape compatible with the input grid'%fieldnames[i])

        # transpose to have the same convention as Paraview (nz, ny, nx, ncomponent)
        if len(fields[i].shape) == 3: # scalair
            fields[i] = np.transpose(fields[i], (2, 1, 0))
        elif len(fields[i].shape) == 4: # vector or tensor
            fields[i] = np.transpose(fields[i], (2, 1, 0, 3))
        else:
            raise ValueError('The input field %s does not have a shape compatible with the input grid'%fieldnames[i])
        

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
            elif field.ndim == 4 and field.shape[-1] == 3:
                atype = "Vector"
                dims = f"{nz} {ny} {nx} 3"
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
