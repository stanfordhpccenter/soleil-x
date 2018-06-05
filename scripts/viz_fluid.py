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
parser.add_argument('hdf_file')
parser.add_argument('json_file')
args = parser.parse_args()

hdf_in = h5py.File(args.hdf_file, 'r')
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

with open(args.json_file) as json_in:
    config = json.load(json_in)
    # CAUTION: We flip the X and Z dimensions, because Legion dumps data in
    # column-major order.
    nx = config['Grid']['zNum']
    ny = config['Grid']['yNum']
    nz = config['Grid']['xNum']
    lx = config['Grid']['zWidth']
    ly = config['Grid']['yWidth']
    lz = config['Grid']['xWidth']
    ox = config['Grid']['origin'][2]
    oy = config['Grid']['origin'][1]
    oz = config['Grid']['origin'][0]

with open('out.xmf', 'w') as xmf_out:
    xmf_out.write(XMF_TEMPLATE
                  .replace('GRID_DIMS', '%s %s %s' % (nx,ny,nz))
                  .replace('GRID_ORIGIN', '%s %s %s' % (ox,oy,oz))
                  .replace('GRID_SPACING', '%s %s %s' % (lx/nx,ly/ny,lz/nz)))
