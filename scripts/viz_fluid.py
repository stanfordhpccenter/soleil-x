#!/usr/bin/env python

import argparse
import h5py
import json

XMF_TEMPLATE = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">
  <Domain>
    <Grid Name="Fluid" GridType="Uniform">
      <!-- Topology: orthonormal 3D grid -->
      <Topology TopologyType="3DCoRectMesh" Dimensions="GRID_DIMS"></Topology>
      <!-- Geometry: Node positions derived implicitly, based on grid origin and cell size -->
      <Geometry GeometryType="Origin_DxDyDz">
        <DataItem Name="Origin" Dimensions="3" NumberType="Float" Precision="8" Format="XML">GRID_ORIGIN</DataItem>
        <DataItem Name="Spacing" Dimensions="3" NumberType="Float" Precision="8" Format="XML">GRID_SPACING</DataItem>
      </Geometry>
      <Attribute Name="pressure" AttributeType="Scalar" Center="Cell">
        <DataItem Dimensions="GRID_DIMS" NumberType="Float" Precision="8" Format="HDF">out.hdf:/pressure</DataItem>
      </Attribute>
      <Attribute Name="rho" AttributeTypen="Scalar" Center="Cell">
        <DataItem Dimensions="GRID_DIMS" NumberType="Float" Precision="8" Format="HDF">out.hdf:/rho</DataItem>
      </Attribute>
      <Attribute Name="temperature" AttributeType="Scalar" Center="Cell">
        <DataItem Dimensions="GRID_DIMS" NumberType="Float" Precision="8" Format="HDF">out.hdf:/temperature</DataItem>
      </Attribute>
      <Attribute Name="velocity" AttributeType="Vector" Center="Cell">
        <DataItem Dimensions="GRID_DIMS 3" NumberType="Float" Precision="8" Format="HDF">out.hdf:/velocity</DataItem>
      </Attribute>
    </Grid>
  </Domain>
</Xdmf>
"""

parser = argparse.ArgumentParser()
parser.add_argument('hdf_file', help='simulation results to visualize')
parser.add_argument('json_file', help='simulation config file')
parser.add_argument('section_num', nargs='?', choices=['1','2'],
                    help='which section to visualize (if multi-section sim)')
args = parser.parse_args()

hdf_in = h5py.File(args.hdf_file, 'r')
# Extract domain size.
nx = hdf_in['pressure'].shape[0]
ny = hdf_in['pressure'].shape[1]
nz = hdf_in['pressure'].shape[2]
hdf_out = h5py.File('out.hdf', 'w')
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
    if args.section_num is not None:
        config = config['configs'][int(args.section_num)-1]
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
with open('out.xmf', 'w') as xmf_out:
    xmf_out.write(XMF_TEMPLATE
                  .replace('GRID_DIMS', '%s %s %s' % (nx+1,ny+1,nz+1))
                  .replace('GRID_ORIGIN', '%s %s %s' % (ox,oy,oz))
                  .replace('GRID_SPACING', '%s %s %s' % (dx,dy,dz)))
