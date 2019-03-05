#!/usr/bin/env python2

import argparse
import h5py
import itertools
import numpy

XMF_HEADER = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">
  <Domain Name="Particles">
    <Grid Name="ParticlesTimeSeries" GridType="Collection" CollectionType="Temporal">
"""

XMF_BODY = """
      <Grid Name="Particles" GridType="Uniform">
        <Time Value="@TIME"/>
        <!-- Topology: One cell for every node, no connectivity -->
        <Topology TopologyType="Polyvertex" NumberOfElements="@NUM_PARTICLES"></Topology>
        <!-- Geometry: Position of every node (particle) given explicitly -->
        <Geometry GeometryType="XYZ">
          <DataItem Dimensions="@NUM_PARTICLES 3" NumberType="Float" Precision="8" Format="HDF">@HDF_FILE:/position</DataItem>
        </Geometry>
        <Attribute Name="temperature" AttributeType="Scalar" Center="Node">
          <DataItem Dimensions="@NUM_PARTICLES" NumberType="Float" Precision="8" Format="HDF">@HDF_FILE:/temperature</DataItem>
        </Attribute>
        <Attribute Name="diameter" AttributeType="Scalar" Center="Node">
          <DataItem Dimensions="@NUM_PARTICLES" NumberType="Float" Precision="8" Format="HDF">@HDF_FILE:/diameter</DataItem>
        </Attribute>
        <Attribute Name="velocity" AttributeType="Vector" Center="Node">
          <DataItem Dimensions="@NUM_PARTICLES 3" NumberType="Float" Precision="8" Format="HDF">@HDF_FILE:/velocity</DataItem>
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
                    help='particle restart file(s) to visualize')
args = parser.parse_args()

size = {}
time = {}
for (f, i) in zip(args.hdf_file, itertools.count()):
    hdf_in = h5py.File(f, 'r')
    hdf_out = h5py.File('particles%010d.hdf' % i, 'w')
    # Read simulation time at timestep
    time[i] = hdf_in.attrs['simTime']
    # Create a validity field, use it to filter the rest of the datasets.
    valids = numpy.array(hdf_in['__valid'][:]) > 0
    size[i] = sum(valids)
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

with open('particles.xmf', 'w') as xmf_out:
    xmf_out.write(XMF_HEADER)
    for i in range(len(args.hdf_file)):
        if size[i] == 0:
            print 'Skipping timestep %s: Paraview cannot handle empty particle files' % i
            continue
        xmf_out.write(XMF_BODY
                      .replace('@TIME', str(time[i]))
                      .replace('@NUM_PARTICLES', str(size[i]))
                      .replace('@HDF_FILE', 'particles%010d.hdf' % i))
    xmf_out.write(XMF_FOOTER)
