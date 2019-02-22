#!/usr/bin/env python2

import argparse
import h5py
import itertools
import json

XMF_HEADER = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">
  <Domain Name="Fluid">
    <Grid Name="FluidTimeSeries" GridType="Collection" CollectionType="Temporal">
"""

XMF_BODY = """
      <Grid Name="Fluid" GridType="Uniform">
        <Time Value="@TIMESTEP"/>
        <!-- Topology: orthonormal 3D grid -->
        <Topology TopologyType="3DCoRectMesh" Dimensions="@POINTS"></Topology>
        <!-- Geometry: Node positions derived implicitly, based on grid origin and cell size -->
        <Geometry GeometryType="Origin_DxDyDz">
          <DataItem Name="Origin" Dimensions="3" NumberType="Float" Precision="8" Format="XML">@GRID_ORIGIN</DataItem>
          <DataItem Name="Spacing" Dimensions="3" NumberType="Float" Precision="8" Format="XML">@GRID_SPACING</DataItem>
        </Geometry>
        <Attribute Name="pressure" AttributeType="Scalar" Center="Cell">
          <DataItem Dimensions="@CELLS" NumberType="Float" Precision="8" Format="HDF">@HDF_FILE:/pressure</DataItem>
        </Attribute>
        <Attribute Name="rho" AttributeTypen="Scalar" Center="Cell">
          <DataItem Dimensions="@CELLS" NumberType="Float" Precision="8" Format="HDF">@HDF_FILE:/rho</DataItem>
        </Attribute>
        <Attribute Name="temperature" AttributeType="Scalar" Center="Cell">
          <DataItem Dimensions="@CELLS" NumberType="Float" Precision="8" Format="HDF">@HDF_FILE:/temperature</DataItem>
        </Attribute>
        <Attribute Name="velocity" AttributeType="Vector" Center="Cell">
          <DataItem Dimensions="@CELLS 3" NumberType="Float" Precision="8" Format="HDF">@HDF_FILE:/velocity</DataItem>
        </Attribute>
      </Grid>
"""

XMF_FOOTER = """
    </Grid>
  </Domain>
</Xdmf>
"""

parser = argparse.ArgumentParser()
parser.add_argument('json_file',
                    help='original simulation configuration file')
parser.add_argument('-s', '--section', choices=['1','2'],
                    help='which section to visualize (if multi-section sim)')
parser.add_argument('hdf_file', nargs='+',
                    help='fluid restart file(s) to visualize')
args = parser.parse_args()

nx = None
ny = None
nz = None

hdf_out_basename = 'out_fluid_'
for (f, i) in zip(args.hdf_file, itertools.count()):
    hdf_in = h5py.File(f, 'r')
    # Extract domain size.
    if nx is None:
        nx = hdf_in['pressure'].shape[0]
        ny = hdf_in['pressure'].shape[1]
        nz = hdf_in['pressure'].shape[2]
    else:
        assert nx == hdf_in['pressure'].shape[0]
        assert ny == hdf_in['pressure'].shape[1]
        assert nz == hdf_in['pressure'].shape[2]
    hdf_out = h5py.File(hdf_out_basename+'%010d.hdf' % i, 'w')
    # Copy pressure over.
    hdf_out['pressure'] = hdf_in['pressure'][:]
    # Copy rho over.
    hdf_out['rho'] = hdf_in['rho'][:]
    # Copy temperature over.
    hdf_out['temperature'] = hdf_in['temperature'][:]
    # Convert velocity from an XxYxZ matrix of triples to an XxYxZx3 matrix.
    hdf_out['velocity'] = hdf_in['velocity'][:][:,:,:,:]
    hdf_out.close()
    hdf_in.close()

# NOTE: We flip the X and Z dimensions, because Legion dumps data in
# column-major order.
with open(args.json_file) as json_in:
    config = json.load(json_in)
    if args.section is not None:
        config = config['configs'][int(args.section)-1]
    # Compute number of boundary cells on each dimension.
    bx = nx - config['Grid']['zNum']
    by = ny - config['Grid']['yNum']
    bz = nz - config['Grid']['xNum']
    assert bx == 0 or bx == 2, 'Expected at most 1-cell boundary'
    assert by == 0 or by == 2, 'Expected at most 1-cell boundary'
    assert bz == 0 or bz == 2, 'Expected at most 1-cell boundary'
    # Compute cell width.
    dx = config['Grid']['zWidth'] / config['Grid']['zNum']
    dy = config['Grid']['yWidth'] / config['Grid']['yNum']
    dz = config['Grid']['xWidth'] / config['Grid']['xNum']
    # Compute grid origin (taking boundary cells into account).
    ox = config['Grid']['origin'][2] - dx * (bx / 2)
    oy = config['Grid']['origin'][1] - dy * (by / 2)
    oz = config['Grid']['origin'][0] - dz * (bz / 2)

# NOTE: The XMF format expects grid dimensions in points, not cells, so we have
# to add 1 on each dimension.
with open('out_fluid.xmf', 'w') as xmf_out:
    xmf_out.write(XMF_HEADER)
    for i in range(len(args.hdf_file)):
        xmf_out.write(XMF_BODY
                      .replace('@TIMESTEP', str(i))
                      .replace('@POINTS', '%s %s %s' % (nx+1,ny+1,nz+1))
                      .replace('@CELLS', '%s %s %s' % (nx,ny,nz))
                      .replace('@GRID_ORIGIN', '%s %s %s' % (ox,oy,oz))
                      .replace('@GRID_SPACING', '%s %s %s' % (dx,dy,dz))
                      .replace('@HDF_FILE', hdf_out_basename+'%010d.hdf' % i))
    xmf_out.write(XMF_FOOTER)
