# Script to generate xdmf file from hdf file (readable with Paraview/Vtk)
# Usage: python generate_xmf.py name_of_hdf_file
import h5py
import numpy as np
import os
import sys

#Xdmf data for fluid (assumes uniform mesh)
xdmfdata = \
"""<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">
  <Domain>
    <Grid Name="MyGrid" GridType="Uniform">
      <Topology TopologyType="3DCORECTMesh" Dimensions="size_z size_y size_x"/>
      <Geometry GeometryType="ORIGIN_DXDYDZ">
        <DataItem Name="Origin" Dimensions="3" NumberType="Float" Format="XML">
                        0 0 0
        </DataItem>
        <DataItem Name="Spacing" Dimensions="3" NumberType="Float" Format="XML">
                        1 1 1
        </DataItem>
      </Geometry>
      <Attribute Name="pressure" AttributeType="Scalar" Center="Node">
        <DataItem Format="HDF" NumberType="Float" Dimensions="size_z size_y size_x">filename:/pressure</DataItem>
      </Attribute>
      <Attribute Name="rho" AttributeType="Scalar" Center="Node">
        <DataItem Format="HDF" NumberType="Float" Dimensions="size_z size_y size_x">filename:/rho</DataItem>
      </Attribute>

      <Attribute Name="Velocity_X" AttributeType="Scalar" Center="Node">
        <DataItem Format="HDF" NumberType="Float" Dimensions="size_z size_y size_x">filename:/Velocity_X</DataItem>
      </Attribute>

      <Attribute Name="Velocity_Y" AttributeType="Scalar" Center="Node">
        <DataItem Format="HDF" NumberType="Float" Dimensions="size_z size_y size_x">filename:/Velocity_Y</DataItem>
      </Attribute>

      <Attribute Name="Velocity_Z" AttributeType="Scalar" Center="Node">
        <DataItem Format="HDF" NumberType="Float" Dimensions="size_z size_y size_x">filename:/Velocity_Z</DataItem>
      </Attribute>

    </Grid>
  </Domain>
</Xdmf>
"""


#File reading
filename = sys.argv[1]
fic = h5py.File( filename )


shape = fic['pressure'].shape
print 'Grid size is ', shape

keys = fic.keys()

#Splitting Velocity in three components, should be possible to have it as a vector though
to_add = ['Velocity_X', 'Velocity_Y', 'Velocity_Z']
for k,key in enumerate(to_add):
  if not( key in keys):
    fic[key] =  fic['velocity'][:][:,:,:,k]
fic.close()

#Opening xmf file for writing
xmffile = open( filename.replace( '.hdf', '.xmf' ) , 'w' )

#Replacing with correct values for grid size and filename:
xdmfdata = xdmfdata.replace( 'size_z', str( shape[0] ) )
xdmfdata = xdmfdata.replace( 'size_y', str( shape[1] ) )
xdmfdata = xdmfdata.replace( 'size_x', str( shape[2] ) )
xdmfdata = xdmfdata.replace( 'filename', filename )

#Dumping xmf
xmffile.write( xdmfdata )
xmffile.close()
