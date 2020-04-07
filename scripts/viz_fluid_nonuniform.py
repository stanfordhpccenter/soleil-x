#!/usr/bin/env python2

import argparse
import h5py
import numpy as np
import itertools
#import json

XMF_HEADER = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">
  <Domain Name="Fluid">
    <Grid Name="FluidTimeSeries" GridType="Collection" CollectionType="Temporal">
"""

XMF_BODY = """
      <Grid Name="Fluid" GridType="Uniform">
        <Time Value="@TIME"/>
        <!-- Topology: orthonormal 3D grid -->
        <Topology TopologyType="3DRectMesh" Dimensions="@POINTS"></Topology>
        <Geometry GeometryType="VXVYVZ">
          <!-- Z Points -->
          <DataItem Dimensions="@NZ_POINTS" NumberType="Float" Precision="8" Format="XML">
            @Z_POINTS
          </DataItem>
          <!-- Y Points -->
          <DataItem Dimensions="@NY_POINTS" NumberType="Float" Precision="8" Format="XML">
            @Y_POINTS
          </DataItem>
          <!-- X Points -->
          <DataItem Dimensions="@NX_POINTS" NumberType="Float" Precision="8" Format="XML">
            @X_POINTS
          </DataItem>
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
parser.add_argument('hdf_file', nargs='+',
                    help='fluid restart file(s) to visualize')
args = parser.parse_args()

nx = None
ny = None
nz = None

time = {}
for (f, i) in zip(args.hdf_file, itertools.count()):
    hdf_in = h5py.File(f, 'r')
    # Read simulation time at timestep
    time[i] = hdf_in.attrs['simTime']
    # Extract domain size.
    if nx is None:
        nx = hdf_in['pressure'].shape[0]
        ny = hdf_in['pressure'].shape[1]
        nz = hdf_in['pressure'].shape[2]
    else:
        assert nx == hdf_in['pressure'].shape[0]
        assert ny == hdf_in['pressure'].shape[1]
        assert nz == hdf_in['pressure'].shape[2]
    hdf_out = h5py.File('fluid%010d.hdf' % i, 'w')
    # Copy pressure over.
    hdf_out['pressure'] = hdf_in['pressure'][:]
    # Copy rho over.
    hdf_out['rho'] = hdf_in['rho'][:]
    # Copy temperature over.
    hdf_out['temperature'] = hdf_in['temperature'][:]
    # Convert velocity from an XxYxZ matrix of triples to an XxYxZx3 matrix.
    hdf_out['velocity'] = hdf_in['velocity'][:][:,:,:,:]

    # Add the gird points to the hdf file
    centerCoordinates = hdf_in['centerCoordinates'][:][:,:,:,:]
    cellWidth = hdf_in['cellWidth'][:][:,:,:,:]

    # Get the points that define the mesh
    x_points = np.squeeze(centerCoordinates[:,0,0,2] - cellWidth[:,0,0,2]/2.0)
    x_points = np.append(x_points, centerCoordinates[-1,0,0,2] + cellWidth[-1,0,0,2]/2.0)

    y_points = np.squeeze(centerCoordinates[0,:,0,1] - cellWidth[0,:,0,1]/2.0)
    y_points = np.append(y_points, centerCoordinates[0,-1,0,1] + cellWidth[0,-1,0,1]/2.0)

    z_points = np.squeeze(centerCoordinates[0,0,:,0] - cellWidth[0,0,:,0]/2.0)
    z_points = np.append(z_points, centerCoordinates[0,0,-1,0] + cellWidth[0,0,-1,0]/2.0)

    hdf_out['x_points'] = x_points
    hdf_out['y_points'] = y_points
    hdf_out['z_points'] = z_points

    hdf_out.close()
    hdf_in.close()

# NOTE: The XMF format expects grid dimensions in points, not cells, so we have
# to add 1 on each dimension.
with open('fluid.xmf', 'w') as xmf_out:
    xmf_out.write(XMF_HEADER)
    for i in range(len(args.hdf_file)):
        xmf_out.write(XMF_BODY
                      .replace('@TIME', str(time[i]))
                      .replace('@POINTS', '%s %s %s' % (nx+1,ny+1,nz+1))
                      .replace('@NX_POINTS', '%s' % (nx+1))
                      .replace('@NY_POINTS', '%s' % (ny+1))
                      .replace('@NZ_POINTS', '%s' % (nz+1))
                      .replace('@X_POINTS', '%s' % np.array2string(x_points)[1:-1])
                      .replace('@Y_POINTS', '%s' % np.array2string(y_points)[1:-1])
                      .replace('@Z_POINTS', '%s' % np.array2string(z_points)[1:-1])
                      .replace('@CELLS', '%s %s %s' % (nx,ny,nz))
                      .replace('@HDF_FILE', 'fluid%010d.hdf' % i))
    xmf_out.write(XMF_FOOTER)
