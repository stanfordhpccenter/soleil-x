#!/usr/bin/env python2

import argparse
import h5py
import numpy

XMF_TEMPLATE = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">
  <Domain>
    <Grid Name="Particles" GridType="Uniform">
      <!-- Topology: One cell for every node, no connectivity -->
      <Topology TopologyType="Polyvertex" NumberOfElements="NUM_PARTICLES"></Topology>
      <!-- Geometry: Position of every node (particle) given explicitly -->
      <Geometry GeometryType="XYZ">
        <DataItem Dimensions="NUM_PARTICLES 3" NumberType="Float" Precision="8" Format="HDF">out.hdf:/position</DataItem>
      </Geometry>
      <Attribute Name="temperature" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="NUM_PARTICLES" NumberType="Float" Precision="8" Format="HDF">out.hdf:/temperature</DataItem>
      </Attribute>
      <Attribute Name="diameter" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="NUM_PARTICLES" NumberType="Float" Precision="8" Format="HDF">out.hdf:/diameter</DataItem>
      </Attribute>
      <Attribute Name="velocity" AttributeType="Vector" Center="Node">
        <DataItem Dimensions="NUM_PARTICLES 3" NumberType="Float" Precision="8" Format="HDF">out.hdf:/velocity</DataItem>
      </Attribute>
    </Grid>
  </Domain>
</Xdmf>
"""

parser = argparse.ArgumentParser()
parser.add_argument('hdf_file', help='simulation results to visualize')
args = parser.parse_args()

hdf_in = h5py.File(args.hdf_file, 'r')
hdf_out = h5py.File('out.hdf', 'w')
# Create a validity field, use it to filter the rest of the datasets.
valids = numpy.array(hdf_in['__valid'][:]) > 0
# Convert position from an array of N triples to an Nx3 matrix.
hdf_out['position'] = hdf_in['position'][valids][:]
# Copy temperature over.
hdf_out['temperature'] = hdf_in['temperature'][valids]
# Copy diameter over.
hdf_out['diameter'] = hdf_in['diameter'][valids]
# Convert velocity from an array of N triples to an Nx3 matrix.
hdf_out['velocity'] = hdf_in['velocity'][valids][:]
hdf_out.close()
hdf_in.close()

with open('out.xmf', 'w') as xmf_out:
    xmf_out.write(XMF_TEMPLATE.replace('NUM_PARTICLES', str(sum(valids))))
