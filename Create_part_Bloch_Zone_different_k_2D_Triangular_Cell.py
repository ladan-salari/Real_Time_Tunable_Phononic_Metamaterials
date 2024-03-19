# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Create a single unit cell of structure with bloch BC.
# The one with a bar in the middle
# that the stiffness changes with voltage.
# Units are mm,kg,s.
# Modeling 2D 6-fold symmetry triangular unit cell!
# Ladan Salari Sharif
# August 2017
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# IMPORT ABAQUS MODULUS
# -----------------------------------------------------------------------
from part import *
from material import *
from assembly import *
from step import *
from mesh import *
from sketch import *
from abaqusConstants import *
from caeModules import *
from interaction import *
from connectorBehavior import *
from section import *

import math
import os
import numpy as np

# --------------------------------------------------------------
# MODEL PARAMETERS
# --------------------------------------------------------------
# All the UNITS are in "mm"
length = 24.10 * 2.0
angle = 45.0  # Angle is in degree
length_winge = 13.41 + 18.2999  # internal wing
magnet_hieght = 12.7
width_magnet = 6.30
thickness_edge_magnet = 0.69
thickness_edge_magnet2 = 1.85
thickness_edge_magnet3 = 1.65
width_winge = 3.0
thickness_main = 1.78 * 2.0
length_extension = 35.0
Material_name = "Acrylic"  # Name the materials
E = 3.2e6  # Young's Modulus (kg/mms^2)
nu = 0.3  # poission Ratio
Dens = 1.18e-6  # Density (kg/mm^3)

Plate_thickess = 4.0


def Model_without_Bloch_BC(Negative_Stiffness):
    # --------------------------------------------------------------
    # SKETCH PART
    # --------------------------------------------------------------
    s = mdb.models["Model-1"].ConstrainedSketch(name="__profile__", sheetSize=0.2)
    s.Line(point1=(0.0, 0.0), point2=(length_winge, 0.0))
    ##        #Point 2->3
    ##    s.Line(point1=(length_winge,0.0),
    ##           point2=(length_winge,magnet_hieght/2.0))
    ##        #Point 3->4
    ##    s.Line(point1=(length_winge,magnet_hieght/2.0),
    ##           point2=(length_winge+width_magnet,magnet_hieght/2.0))
    ##        #Point 4->5
    ##    s.Line(point1=(length_winge+width_magnet,magnet_hieght/2.0),
    ##           point2=(length_winge+width_magnet,0.0))
    ##        #Point 5->6
    ##    s.Line(point1=(length_winge+width_magnet,0.0),
    ##           point2=(length_winge+width_magnet+thickness_edge_magnet,0.0))
    # Point 2->6
    s.Line(
        point1=(length_winge, 0.0),
        point2=(length_winge + width_magnet + thickness_edge_magnet, 0.0),
    )
    # point 6->7
    s.Line(
        point1=(length_winge + width_magnet + thickness_edge_magnet, 0.0),
        point2=(
            length_winge + width_magnet + thickness_edge_magnet,
            magnet_hieght / 2.0 + thickness_edge_magnet2,
        ),
    )
    # point 7->8
    s.Line(
        point1=(
            length_winge + width_magnet + thickness_edge_magnet,
            magnet_hieght / 2.0 + thickness_edge_magnet2,
        ),
        point2=(
            length_winge - thickness_edge_magnet3,
            magnet_hieght / 2.0 + thickness_edge_magnet2,
        ),
    )
    # point 8->9
    s.Line(
        point1=(
            length_winge - thickness_edge_magnet3,
            magnet_hieght / 2.0 + thickness_edge_magnet2,
        ),
        point2=(length_winge - thickness_edge_magnet3, width_winge),
    )
    # point 9->10
    s.Line(
        point1=(length_winge - thickness_edge_magnet3, width_winge),
        point2=(thickness_main * math.sqrt(2), width_winge),
    )
    # point 10->11
    s.Line(
        point1=(thickness_main * math.sqrt(2.0), width_winge),
        point2=(
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180),
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0,
        ),
    )
    # point 11->12
    ##    s.Line(point1=(thickness_main*math.sqrt(2.0)/2.0+length*math.cos(angle*math.pi/180),width_winge+length*math.sin(angle*math.pi/180)-thickness_main*math.sqrt(2.0)/2.0),
    ##           point2=(thickness_main*math.sqrt(2.0)/2.0+length*math.cos(angle*math.pi/180),width_winge+length*math.sin(angle*math.pi/180)-thickness_main*math.sqrt(2.0)/2.0-width_winge))
    # point 11->13
    s.Line(
        point1=(
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180),
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0,
        ),
        point2=(
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180)
            + width_winge,
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0,
        ),
    )
    # point 13->14
    s.Line(
        point1=(
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180)
            + width_winge,
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0,
        ),
        point2=(
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180)
            + width_winge,
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0
            + thickness_main * math.sqrt(2.0),
        ),
    )
    # point 14->17
    s.Line(
        point1=(
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180)
            + width_winge,
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0
            + thickness_main * math.sqrt(2.0),
        ),
        point2=(
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180)
            + width_winge,
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0
            + thickness_main * math.sqrt(2.0)
            + length_extension,
        ),
    )
    # point 17->18
    s.Line(
        point1=(
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180)
            + width_winge,
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0
            + thickness_main * math.sqrt(2.0)
            + length_extension,
        ),
        point2=(
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180),
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0
            + thickness_main * math.sqrt(2.0)
            + length_extension,
        ),
    )
    # point 18->15
    s.Line(
        point1=(
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180),
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0
            + thickness_main * math.sqrt(2.0)
            + length_extension,
        ),
        point2=(
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180),
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0
            + thickness_main * math.sqrt(2.0),
        ),
    )
    # point 15->16
    s.Line(
        point1=(
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180),
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0
            + thickness_main * math.sqrt(2.0),
        ),
        point2=(0.0, width_winge),
    )
    # point 16->1
    s.Line(point1=(0.0, width_winge), point2=(0.0, 0.0))
    mdb.models["Model-1"].Part(
        dimensionality=TWO_D_PLANAR, name="Part-1", type=DEFORMABLE_BODY
    )
    mdb.models["Model-1"].parts["Part-1"].BaseShell(
        sketch=mdb.models["Model-1"].sketches["__profile__"]
    )
    p = mdb.models["Model-1"].parts["Part-1"]
    # partition the two magnets
    Sketch = (
        mdb.models["Model-1"]
        .sketches["__profile__"]
        .Line(
            point1=(length_winge, -magnet_hieght / 2.0),
            point2=(length_winge, magnet_hieght / 2.0),
        )
    )
    Sketch = (
        mdb.models["Model-1"]
        .sketches["__profile__"]
        .Line(
            point1=(length_winge, magnet_hieght / 2.0),
            point2=(length_winge + width_magnet, magnet_hieght / 2.0),
        )
    )
    Sketch = (
        mdb.models["Model-1"]
        .sketches["__profile__"]
        .Line(
            point1=(length_winge, -magnet_hieght / 2.0),
            point2=(length_winge + width_magnet, -magnet_hieght / 2.0),
        )
    )
    Sketch = (
        mdb.models["Model-1"]
        .sketches["__profile__"]
        .Line(
            point1=(length_winge + width_magnet, magnet_hieght / 2.0),
            point2=(length_winge + width_magnet, -magnet_hieght / 2.0),
        )
    )
    p.PartitionFaceBySketch(
        faces=p.faces.findAt(
            (0.0, 0.0, 0.0),
        ),
        sketch=s,
    )
    # mirror the part along y axias
    mydatum = (
        mdb.models["Model-1"]
        .parts["Part-1"]
        .DatumPlaneByPrincipalPlane(
            principalPlane=YZPLANE,
            offset=thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180)
            + width_winge,
        )
        .id
    )
    Datummm = mdb.models["Model-1"].parts["Part-1"].datums
    mdb.models["Model-1"].parts["Part-1"].Mirror(
        mirrorPlane=Datummm[mydatum], keepOriginal=ON
    )
    # mirror the part along x axias
    mydatum = (
        mdb.models["Model-1"]
        .parts["Part-1"]
        .DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
        .id
    )
    Datummm = mdb.models["Model-1"].parts["Part-1"].datums
    mdb.models["Model-1"].parts["Part-1"].Mirror(
        mirrorPlane=Datummm[mydatum], keepOriginal=ON
    )

    s = mdb.models["Model-1"].ConstrainedSketch(name="__profile2__", sheetSize=0.2)
    s.Line(point1=(0.0, 0.0), point2=(2.0 * width_winge, 0.0))
    s.Line(
        point1=(2.0 * width_winge, 0.0),
        point2=(width_winge, width_winge * math.sqrt(3.0)),
    )
    s.Line(
        point1=(width_winge, width_winge * math.sqrt(3.0)),
        point2=(-width_winge, width_winge * math.sqrt(3.0)),
    )
    s.Line(point1=(-width_winge, width_winge * math.sqrt(3.0)), point2=(0.0, 0.0))

    mdb.models["Model-1"].Part(
        dimensionality=TWO_D_PLANAR, name="Part-2", type=DEFORMABLE_BODY
    )
    mdb.models["Model-1"].parts["Part-2"].BaseShell(
        sketch=mdb.models["Model-1"].sketches["__profile2__"]
    )
    p2 = mdb.models["Model-1"].parts["Part-2"]

    # --------------------------------------------------------------
    # DEFINE MATERIAL
    # --------------------------------------------------------------
    mdb.models["Model-1"].Material(name=Material_name)
    mdb.models["Model-1"].materials[Material_name].Elastic(table=((E, nu),))
    mdb.models["Model-1"].materials[Material_name].Density(table=((Dens,),))
    mdb.models["Model-1"].HomogeneousSolidSection(
        name="Solid-section", material=Material_name, thickness=Plate_thickess
    )
    mdb.models["Model-1"].Material(name="magnet")
    mdb.models["Model-1"].materials["magnet"].Elastic(table=((210e6, nu),))
    mdb.models["Model-1"].materials["magnet"].Density(table=((7.5804e-06,),))
    mdb.models["Model-1"].HomogeneousSolidSection(
        name="Solid-section_magnet", material="magnet", thickness=width_magnet
    )
    p = mdb.models["Model-1"].parts["Part-1"]
    f = p.faces.findAt(((0.0, 0.0, 0.0),))
    region = regionToolset.Region(faces=f)
    p.SectionAssignment(region=region, sectionName="Solid-section")
    f = p.faces.findAt(((length_winge + width_magnet / 2.0, 0.0, 0.0),))
    region = regionToolset.Region(faces=f)
    p.SectionAssignment(region=region, sectionName="Solid-section_magnet")
    f = p.faces.findAt(((length_winge + width_magnet / 2.0 + 7.0, 0.0, 0.0),))
    region = regionToolset.Region(faces=f)
    p.SectionAssignment(region=region, sectionName="Solid-section_magnet")
    # --------------------------------------------------------------
    # ASSEMBLY
    # --------------------------------------------------------------
    a = mdb.models["Model-1"].rootAssembly
    session.viewports["Viewport: 1"].setValues(displayedObject=a)
    session.viewports["Viewport: 1"].assemblyDisplay.setValues(mesh=ON)
    session.viewports["Viewport: 1"].assemblyDisplay.meshOptions.setValues(
        meshTechnique=OFF
    )
    a = mdb.models["Model-1"].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mdb.models["Model-1"].parts["Part-1"]
    a.Instance(name="Part1-1", part=p, dependent=ON)
    a.Instance(name="Part1-2", part=p, dependent=ON, autoOffset=OFF)
    a.Instance(name="Part1-3", part=p, dependent=ON, autoOffset=OFF)
    # rotate the instance
    a.rotate(
        instanceList=("Part1-1",),
        axisPoint=(0.0, 0.0, 0.0),
        axisDirection=(0.0, 0.0, 1.0),
        angle=90.0,
    )
    a.rotate(
        instanceList=("Part1-2",),
        axisPoint=(0.0, 0.0, 0.0),
        axisDirection=(0.0, 0.0, 1.0),
        angle=-30.0,
    )
    a.rotate(
        instanceList=("Part1-3",),
        axisPoint=(0.0, 0.0, 0.0),
        axisDirection=(0.0, 0.0, 1.0),
        angle=30.0,
    )
    # translate the instance to the right possition
    a.translate(instanceList=("Part1-2",), vector=(-74.192473, 128.505133, 0.0))
    a.translate(instanceList=("Part1-3",), vector=(5.603526, 88.905286, 0.0))

    a.Instance(name="Part2-1", part=p2, dependent=ON, autoOffset=OFF)
    a.Instance(name="Part2-2", part=p2, dependent=ON, autoOffset=OFF)
    a.Instance(name="Part2-3", part=p2, dependent=ON, autoOffset=OFF)
    # rotate the instance
    a.rotate(
        instanceList=("Part2-1",),
        axisPoint=(0.0, 0.0, 0.0),
        axisDirection=(0.0, 0.0, 1.0),
        angle=90.0,
    )
    a.rotate(
        instanceList=("Part2-2",),
        axisPoint=(0.0, 0.0, 0.0),
        axisDirection=(0.0, 0.0, 1.0),
        angle=30.0,
    )
    a.rotate(
        instanceList=("Part2-3",),
        axisPoint=(0.0, 0.0, 0.0),
        axisDirection=(0.0, 0.0, 1.0),
        angle=-30.0,
    )
    a.translate(instanceList=("Part2-2",), vector=(0.0, 171.810572, 0.0))
    a.translate(instanceList=("Part2-3",), vector=(-79.795999, 39.599847, 0.0))
    a.translate(instanceList=("Part2-1",), vector=(79.795999, 39.599847, 0.0))
    # merge:
    SingleInstances_List = a.instances.keys()
    a.InstanceFromBooleanMerge(
        name="Unit_cell",
        instances=(
            [
                a.instances[SingleInstances_List[i]]
                for i in range(len(SingleInstances_List))
            ]
        ),
    )

    p3 = mdb.models["Model-1"].parts["Unit_cell"]
    a.Instance(name="Unit_cell-1", part=p3, dependent=ON, autoOffset=ON)
    a.Instance(name="Unit_cell-2", part=p3, dependent=ON, autoOffset=ON)
    # --------------------------------------------------------------
    # MESH
    # --------------------------------------------------------------
    # partition the edges magnets
    p3.PartitionFaceByShortestPath(
        faces=p3.faces.findAt(
            (0.0, 0.0, 0.0),
        ),
        point1=(
            -(
                width_winge
                + length * math.sin(angle * math.pi / 180)
                - thickness_main * math.sqrt(2.0) / 2.0
                + thickness_main * math.sqrt(2.0)
                + length_extension
            ),
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180),
            0.0,
        ),
        point2=(
            -(
                width_winge
                + length * math.sin(angle * math.pi / 180)
                - thickness_main * math.sqrt(2.0) / 2.0
                + thickness_main * math.sqrt(2.0)
                + length_extension
            ),
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180)
            + 2.0 * width_winge,
            0.0,
        ),
    )
    p3.PartitionFaceByShortestPath(
        faces=p3.faces.findAt(
            (0.0, 0.0, 0.0),
        ),
        point1=(
            -(
                width_winge
                + length * math.sin(angle * math.pi / 180)
                - thickness_main * math.sqrt(2.0) / 2.0
                + thickness_main * math.sqrt(2.0)
            ),
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180),
            0.0,
        ),
        point2=(
            -(
                width_winge
                + length * math.sin(angle * math.pi / 180)
                - thickness_main * math.sqrt(2.0) / 2.0
                + thickness_main * math.sqrt(2.0)
            ),
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180)
            + 2.0 * width_winge,
            0.0,
        ),
    )
    p3.PartitionFaceByShortestPath(
        faces=p3.faces.findAt(
            (0.0, 0.0, 0.0),
        ),
        point1=(
            (
                width_winge
                + length * math.sin(angle * math.pi / 180)
                - thickness_main * math.sqrt(2.0) / 2.0
                + thickness_main * math.sqrt(2.0)
                + length_extension
            ),
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180),
            0.0,
        ),
        point2=(
            (
                width_winge
                + length * math.sin(angle * math.pi / 180)
                - thickness_main * math.sqrt(2.0) / 2.0
                + thickness_main * math.sqrt(2.0)
                + length_extension
            ),
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180)
            + 2.0 * width_winge,
            0.0,
        ),
    )
    p3.PartitionFaceByShortestPath(
        faces=p3.faces.findAt(
            (0.0, 0.0, 0.0),
        ),
        point1=(
            (
                width_winge
                + length * math.sin(angle * math.pi / 180)
                - thickness_main * math.sqrt(2.0) / 2.0
                + thickness_main * math.sqrt(2.0)
            ),
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180),
            0.0,
        ),
        point2=(
            (
                width_winge
                + length * math.sin(angle * math.pi / 180)
                - thickness_main * math.sqrt(2.0) / 2.0
                + thickness_main * math.sqrt(2.0)
            ),
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180)
            + 2.0 * width_winge,
            0.0,
        ),
    )
    p3.PartitionFaceByAuto(
        face=p3.faces.findAt(
            (0.0, 0.0, 0.0),
        )
    )
    p3.PartitionFaceByShortestPath(
        faces=p3.faces.findAt(
            (-79.795999, 39.599847, 0.0),
        ),
        point1=(74.599847, 42.599847, 0.0),
        point2=(79.795999, 45.599847, 0.0),
    )
    p3.PartitionFaceByShortestPath(
        faces=p3.faces.findAt(
            (-79.795999, 39.599847, 0.0),
        ),
        point1=(57.099847, 72.910736, 0.0),
        point2=(62.295999, 75.910736, 0.0),
    )
    p3.PartitionFaceByShortestPath(
        faces=p3.faces.findAt(
            (-79.795999, 39.599847, 0.0),
        ),
        point1=(17.5, 141.499683, 0.0),
        point2=(22.696152, 144.499683, 0.0),
    )
    p3.PartitionFaceByShortestPath(
        faces=p3.faces.findAt(
            (-79.795999, 39.599847, 0.0),
        ),
        point1=(5.196152, 174.810572, 0.0),
        point2=(0.0, 171.810572, 0.0),
    )
    p3.PartitionFaceByShortestPath(
        faces=p3.faces.findAt(
            (-79.795999, 39.599847, 0.0),
        ),
        point1=(-5.196152, 174.810572, 0.0),
        point2=(0.0, 171.810572, 0.0),
    )
    p3.PartitionFaceByShortestPath(
        faces=p3.faces.findAt(
            (-79.795999, 39.599847, 0.0),
        ),
        point1=(-17.5, 141.499683, 0.0),
        point2=(-22.696152, 144.499683, 0.0),
    )
    p3.PartitionFaceByShortestPath(
        faces=p3.faces.findAt(
            (-79.795999, 39.599847, 0.0),
        ),
        point1=(-62.295999, 75.910736, 0.0),
        point2=(-57.099847, 72.910736, 0.0),
    )
    p3.PartitionFaceByShortestPath(
        faces=p3.faces.findAt(
            (-79.795999, 39.599847, 0.0),
        ),
        point1=(-79.795999, 45.599847, 0.0),
        point2=(-74.599847, 42.599847, 0.0),
    )
    p3.PartitionFaceByAuto(
        face=p3.faces.findAt(
            (-75.692473, 125.907057, 0.0),
        )
    )
    p3.PartitionFaceByAuto(
        face=p3.faces.findAt(
            (75.692473, 125.907057, 0.0),
        )
    )
    p3.seedPart(deviationFactor=0.1, size=0.5)
    p3.setMeshControls(
        elemShape=QUAD_DOMINATED,
        regions=p3.faces.findAt(((0.0, 0.0, 0.0),)),
        technique=FREE,
    )
    p3.setElementType(
        elemTypes=(
            ElemType(elemCode=CPS4R, elemLibrary=STANDARD),
            ElemType(elemCode=CPS3, elemLibrary=STANDARD),
        ),
        regions=(p3.faces.findAt(((0.0, 0.0, 0.0),)),),
    )
    p3.generateMesh()
    # --------------------------------------------------------------
    # step
    # --------------------------------------------------------------
    mdb.models["Model-1"].FrequencyStep(
        name="Frequency",
        previous="Initial",
        eigensolver=LANCZOS,
        numEigen=20,
        blockSize=DEFAULT,
    )
    # --------------------------------------------------------------
    # Define Constrain between plate and Ref nodes
    # --------------------------------------------------------------
    ##    #leg 1:
    edges1 = (
        mdb.models["Model-1"]
        .parts["Unit_cell"]
        .edges.getByBoundingCylinder(
            center1=(magnet_hieght / 2.0, length_winge + width_magnet, 0.0),
            center2=(-magnet_hieght / 2.0, length_winge + width_magnet, 0.0),
            radius=1e-4,
        )
    )
    p3.Set(edges=edges1, name="Leg1_right_Plate")

    edges2 = (
        mdb.models["Model-1"]
        .parts["Unit_cell"]
        .edges.getByBoundingCylinder(
            center1=(magnet_hieght / 2.0, length_winge + width_magnet + 3.179894, 0.0),
            center2=(-magnet_hieght / 2.0, length_winge + width_magnet + 3.179894, 0.0),
            radius=1e-4,
        )
    )
    p3.Set(edges=edges2, name="Leg1_Left_Plate")

    # leg 2:
    edges3 = (
        mdb.models["Model-1"]
        .parts["Unit_cell"]
        .edges.getByBoundingCylinder(
            center1=(-41.696065, 102.410975, 0.0),
            center2=(-35.346065, 113.409497, 0.0),
            radius=1e-4,
        )
    )
    p3.Set(edges=edges3, name="Leg2_right_Plate")
    edges4 = (
        mdb.models["Model-1"]
        .parts["Unit_cell"]
        .edges.getByBoundingCylinder(
            center1=(-44.449934, 104.000922, 0.0),
            center2=(-38.099934, 114.999444, 0.0),
            radius=1e-4,
        )
    )
    p3.Set(edges=edges4, name="Leg2_Left_Plate")

    # leg 3:
    edges5 = (
        mdb.models["Model-1"]
        .parts["Unit_cell"]
        .edges.getByBoundingCylinder(
            center1=(41.696065, 102.410975, 0.0),
            center2=(35.346065, 113.409497, 0.0),
            radius=1e-4,
        )
    )
    p3.Set(edges=edges5, name="Leg3_Left_Plate")
    edges6 = (
        mdb.models["Model-1"]
        .parts["Unit_cell"]
        .edges.getByBoundingCylinder(
            center1=(44.449934, 104.000922, 0.0),
            center2=(38.099934, 114.999444, 0.0),
            radius=1e-4,
        )
    )
    p3.Set(edges=edges6, name="Leg3_right_Plate")
    # Ref node of the real leg 1:
    id_number1 = a.ReferencePoint(point=(95.755199, 38.0099, 0.0)).id
    ref1 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints1 = (ref1[id_number1],)
    a.Set(referencePoints=refPoints1, name="Leg1_Right_Node_Real")

    id_number2 = a.ReferencePoint(point=(95.755199, 41.189794, 0.0)).id
    ref1 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints1 = (ref1[id_number2],)
    a.Set(referencePoints=refPoints1, name="Leg1_Left_Node_Real")

    # Ref node of the real leg 2:
    id_number3 = a.ReferencePoint(point=(57.234134, 107.910236, 0.0)).id
    ref1 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints1 = (ref1[id_number3],)
    a.Set(referencePoints=refPoints1, name="Leg2_Right_Node_Real")

    id_number4 = a.ReferencePoint(point=(54.480265, 109.500183, 0.0)).id
    ref1 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints1 = (ref1[id_number4],)
    a.Set(referencePoints=refPoints1, name="Leg2_Left_Node_Real")

    # Ref node of the real leg 3:
    id_number5 = a.ReferencePoint(point=(137.030133, 109.500183, 0.0)).id
    ref1 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints1 = (ref1[id_number5],)
    a.Set(referencePoints=refPoints1, name="Leg3_Right_Node_Real")

    id_number6 = a.ReferencePoint(point=(134.276264, 107.910236, 0.0)).id
    ref1 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints1 = (ref1[id_number6],)
    a.Set(referencePoints=refPoints1, name="Leg3_Left_Node_Real")
    # Ref node of the Imag leg 1:
    id_number7 = a.ReferencePoint(point=(271.306397, 38.0099, 0.0)).id
    ref1 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints1 = (ref1[id_number7],)
    a.Set(referencePoints=refPoints1, name="Leg1_Right_Node_Imag")

    id_number8 = a.ReferencePoint(point=(271.306397, 41.189794, 0.0)).id
    ref1 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints1 = (ref1[id_number8],)
    a.Set(referencePoints=refPoints1, name="Leg1_Left_Node_Imag")

    # Ref node of the Imag leg 2:
    id_number9 = a.ReferencePoint(point=(232.785332, 107.910236, 0.0)).id
    ref1 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints1 = (ref1[id_number9],)
    a.Set(referencePoints=refPoints1, name="Leg2_Right_Node_Imag")

    id_number10 = a.ReferencePoint(point=(230.031463, 109.500183, 0.0)).id
    ref1 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints1 = (ref1[id_number10],)
    a.Set(referencePoints=refPoints1, name="Leg2_Left_Node_Imag")

    # Ref node of the Imag leg 3:
    id_number11 = a.ReferencePoint(point=(312.581331, 109.500183, 0.0)).id
    ref1 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints1 = (ref1[id_number11],)
    a.Set(referencePoints=refPoints1, name="Leg3_Right_Node_Imag")

    id_number12 = a.ReferencePoint(point=(309.827462, 107.910236, 0.0)).id
    ref1 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints1 = (ref1[id_number12],)
    a.Set(referencePoints=refPoints1, name="Leg3_Left_Node_Imag")

    # Coupling Real part
    mdb.models["Model-1"].Coupling(
        name="Leg1_Left-Real",
        controlPoint=a.sets["Leg1_Left_Node_Real"],
        surface=a.sets["Unit_cell-1.Leg1_Left_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )

    mdb.models["Model-1"].Coupling(
        name="Leg1_right-Real",
        controlPoint=a.sets["Leg1_Right_Node_Real"],
        surface=a.sets["Unit_cell-1.Leg1_right_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )

    mdb.models["Model-1"].Coupling(
        name="Leg2_Left-Real",
        controlPoint=a.sets["Leg2_Left_Node_Real"],
        surface=a.sets["Unit_cell-1.Leg2_Left_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )

    mdb.models["Model-1"].Coupling(
        name="Leg2_right-Real",
        controlPoint=a.sets["Leg2_Right_Node_Real"],
        surface=a.sets["Unit_cell-1.Leg2_right_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )

    mdb.models["Model-1"].Coupling(
        name="Leg3_Left-Real",
        controlPoint=a.sets["Leg3_Left_Node_Real"],
        surface=a.sets["Unit_cell-1.Leg3_Left_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )

    mdb.models["Model-1"].Coupling(
        name="Leg3_right-Real",
        controlPoint=a.sets["Leg3_Right_Node_Real"],
        surface=a.sets["Unit_cell-1.Leg3_right_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )

    mdb.models["Model-1"].Coupling(
        name="Leg1_Left-Imag",
        controlPoint=a.sets["Leg1_Left_Node_Imag"],
        surface=a.sets["Unit_cell-2.Leg1_Left_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )

    mdb.models["Model-1"].Coupling(
        name="Leg1_right-Imag",
        controlPoint=a.sets["Leg1_Right_Node_Imag"],
        surface=a.sets["Unit_cell-2.Leg1_right_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )

    mdb.models["Model-1"].Coupling(
        name="Leg2_Left-Imag",
        controlPoint=a.sets["Leg2_Left_Node_Imag"],
        surface=a.sets["Unit_cell-2.Leg2_Left_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )

    mdb.models["Model-1"].Coupling(
        name="Leg2_right-Imag",
        controlPoint=a.sets["Leg2_Right_Node_Imag"],
        surface=a.sets["Unit_cell-2.Leg2_right_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )

    mdb.models["Model-1"].Coupling(
        name="Leg3_Left-Imag",
        controlPoint=a.sets["Leg3_Left_Node_Imag"],
        surface=a.sets["Unit_cell-2.Leg3_Left_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )

    mdb.models["Model-1"].Coupling(
        name="Leg3_right-Imag",
        controlPoint=a.sets["Leg3_Right_Node_Imag"],
        surface=a.sets["Unit_cell-2.Leg3_right_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )

    # --------------------------------------------------------------
    # Define Wire
    # --------------------------------------------------------------
    # UNFORTUNATLY EACH TIME I NEED TO FIND OUT THE DISTSNCE BETWEEN THE TWO PART AND PUT IT ISTEAD OF 62.22 LINE!
    mdb.models["Model-1"].rootAssembly.WirePolyLine(
        points=(
            (
                mdb.models["Model-1"].rootAssembly.referencePoints[id_number1],
                mdb.models["Model-1"].rootAssembly.referencePoints[id_number2],
            )
        ),
        meshable=OFF,
    )
    edges1 = mdb.models["Model-1"].rootAssembly.edges.getSequenceFromMask(
        mask=("[#1 ]",),
    )
    mdb.models["Model-1"].rootAssembly.Set(edges=edges1, name="Square1_wire_real")

    mdb.models["Model-1"].rootAssembly.WirePolyLine(
        points=(
            (
                mdb.models["Model-1"].rootAssembly.referencePoints[id_number3],
                mdb.models["Model-1"].rootAssembly.referencePoints[id_number4],
            )
        ),
        meshable=OFF,
    )
    edges1 = mdb.models["Model-1"].rootAssembly.edges.getSequenceFromMask(
        mask=("[#1 ]",),
    )
    mdb.models["Model-1"].rootAssembly.Set(edges=edges1, name="Square2_wire_real")

    mdb.models["Model-1"].rootAssembly.WirePolyLine(
        points=(
            (
                mdb.models["Model-1"].rootAssembly.referencePoints[id_number5],
                mdb.models["Model-1"].rootAssembly.referencePoints[id_number6],
            )
        ),
        meshable=OFF,
    )
    edges1 = mdb.models["Model-1"].rootAssembly.edges.getSequenceFromMask(
        mask=("[#1 ]",),
    )
    mdb.models["Model-1"].rootAssembly.Set(edges=edges1, name="Square3_wire_real")

    mdb.models["Model-1"].rootAssembly.WirePolyLine(
        points=(
            (
                mdb.models["Model-1"].rootAssembly.referencePoints[id_number7],
                mdb.models["Model-1"].rootAssembly.referencePoints[id_number8],
            )
        ),
        meshable=OFF,
    )
    edges1 = mdb.models["Model-1"].rootAssembly.edges.getSequenceFromMask(
        mask=("[#1 ]",),
    )
    mdb.models["Model-1"].rootAssembly.Set(edges=edges1, name="Square1_wire_Imag")

    mdb.models["Model-1"].rootAssembly.WirePolyLine(
        points=(
            (
                mdb.models["Model-1"].rootAssembly.referencePoints[id_number9],
                mdb.models["Model-1"].rootAssembly.referencePoints[id_number10],
            )
        ),
        meshable=OFF,
    )
    edges1 = mdb.models["Model-1"].rootAssembly.edges.getSequenceFromMask(
        mask=("[#1 ]",),
    )
    mdb.models["Model-1"].rootAssembly.Set(edges=edges1, name="Square2_wire_Imag")

    mdb.models["Model-1"].rootAssembly.WirePolyLine(
        points=(
            (
                mdb.models["Model-1"].rootAssembly.referencePoints[id_number11],
                mdb.models["Model-1"].rootAssembly.referencePoints[id_number12],
            )
        ),
        meshable=OFF,
    )
    edges1 = mdb.models["Model-1"].rootAssembly.edges.getSequenceFromMask(
        mask=("[#1 ]",),
    )
    mdb.models["Model-1"].rootAssembly.Set(edges=edges1, name="Square3_wire_Imag")

    # --------------------------------------------------------------
    # Define Connector
    # --------------------------------------------------------------

    mdb.models["Model-1"].ConnectorSection(
        name="Negative_conector", translationalType=AXIAL
    )
    stop_0 = connectorBehavior.ConnectorElasticity(
        behavior=NONLINEAR,
        coupling=UNCOUPLED,
        components=(1,),
        table=((-Negative_Stiffness, 1.0), (Negative_Stiffness, -1.0)),
    )
    mdb.models["Model-1"].sections["Negative_conector"].setValues(
        behaviorOptions=(stop_0,)
    )
    mdb.models["Model-1"].rootAssembly.SectionAssignment(
        sectionName="Negative_conector", region=a.sets["Square1_wire_real"]
    )
    mdb.models["Model-1"].rootAssembly.SectionAssignment(
        sectionName="Negative_conector", region=a.sets["Square1_wire_Imag"]
    )
    mdb.models["Model-1"].rootAssembly.SectionAssignment(
        sectionName="Negative_conector", region=a.sets["Square2_wire_real"]
    )
    mdb.models["Model-1"].rootAssembly.SectionAssignment(
        sectionName="Negative_conector", region=a.sets["Square2_wire_Imag"]
    )
    mdb.models["Model-1"].rootAssembly.SectionAssignment(
        sectionName="Negative_conector", region=a.sets["Square3_wire_real"]
    )
    mdb.models["Model-1"].rootAssembly.SectionAssignment(
        sectionName="Negative_conector", region=a.sets["Square3_wire_Imag"]
    )
    # --------------------------------------------------------------
    # Define node set
    # --------------------------------------------------------------
    # Left Face in a2 direction
    delta = 1e-4
    Face1 = (
        mdb.models["Model-1"]
        .parts["Unit_cell"]
        .nodes.getByBoundingCylinder(
            center1=(-79.795999, 39.599847, 0.0),
            center2=(-74.599847, 36.599847, 0.0),
            radius=delta,
        )
    )

    # mdb.models['Model-1'].parts['Unit_cell'].NodeSet(name='P1', nodes=
    Face1_sort = sorted(Face1, key=lambda Face1: Face1.coordinates[1])
    Face1_Node = []
    mylist = list(xrange(len(Face1_sort)))
    for index in mylist:
        Face1_Node.append(Face1_sort[index].label)
        mdb.models["Model-1"].parts["Unit_cell"].SetFromNodeLabels(
            nodeLabels=(Face1_sort[index].label,), name="Left_a2-" + str(index)
        )

    # Right Face in a2 direction
    delta = 1e-4
    Face1 = (
        mdb.models["Model-1"]
        .parts["Unit_cell"]
        .nodes.getByBoundingCylinder(
            center1=(5.196152, 174.810572, 0.0),
            center2=(0.0, 177.810572, 0.0),
            radius=delta,
        )
    )

    # mdb.models['Model-1'].parts['Unit_cell'].NodeSet(name='P1', nodes=
    Face1_sort = sorted(Face1, key=lambda Face1: Face1.coordinates[1])
    Face1_Node = []
    mylist = list(xrange(len(Face1_sort)))
    for index in mylist:
        Face1_Node.append(Face1_sort[index].label)
        mdb.models["Model-1"].parts["Unit_cell"].SetFromNodeLabels(
            nodeLabels=(Face1_sort[index].label,), name="right_a2-" + str(index)
        )

    # Left Face in a1
    delta = 1e-4
    Face1 = (
        mdb.models["Model-1"]
        .parts["Unit_cell"]
        .nodes.getByBoundingCylinder(
            center1=(-79.795999, 39.599847, 0.0),
            center2=(-79.795999, 45.599847, 0.0),
            radius=delta,
        )
    )

    Face1_sort = sorted(Face1, key=lambda Face1: Face1.coordinates[1])
    Face1_Node = []
    mylist = list(xrange(len(Face1_sort)))
    for index in mylist:
        Face1_Node.append(Face1_sort[index].label)
        mdb.models["Model-1"].parts["Unit_cell"].SetFromNodeLabels(
            nodeLabels=(Face1_sort[index].label,), name="Left_a1-" + str(index)
        )

    # Right Face in a1 direction
    delta = 1e-4
    Face1 = (
        mdb.models["Model-1"]
        .parts["Unit_cell"]
        .nodes.getByBoundingCylinder(
            center1=(79.795999, 39.599847, 0.0),
            center2=(79.795999, 45.599847, 0.0),
            radius=delta,
        )
    )

    # mdb.models['Model-1'].parts['Unit_cell'].NodeSet(name='P1', nodes=
    Face1_sort = sorted(Face1, key=lambda Face1: Face1.coordinates[1])
    Face1_Node = []
    mylist = list(xrange(len(Face1_sort)))
    for index in mylist:
        Face1_Node.append(Face1_sort[index].label)
        mdb.models["Model-1"].parts["Unit_cell"].SetFromNodeLabels(
            nodeLabels=(Face1_sort[index].label,), name="right_a1-" + str(index)
        )

    # Left Face in a2-a1 direction
    delta = 1e-4
    Face1 = (
        mdb.models["Model-1"]
        .parts["Unit_cell"]
        .nodes.getByBoundingCylinder(
            center1=(-5.196152, 174.810572, 0.0),
            center2=(0.0, 177.810572, 0.0),
            radius=delta,
        )
    )

    # mdb.models['Model-1'].parts['Unit_cell'].NodeSet(name='P1', nodes=
    Face1_sort = sorted(Face1, key=lambda Face1: Face1.coordinates[1])
    Face1_Node = []
    mylist = list(xrange(len(Face1_sort)))
    for index in mylist:
        Face1_Node.append(Face1_sort[index].label)
        mdb.models["Model-1"].parts["Unit_cell"].SetFromNodeLabels(
            nodeLabels=(Face1_sort[index].label,), name="Left_a2-a1-" + str(index)
        )

    # right Face in a2-a1 direction
    delta = 1e-4
    Face1 = (
        mdb.models["Model-1"]
        .parts["Unit_cell"]
        .nodes.getByBoundingCylinder(
            center1=(74.599847, 36.599847, 0.0),
            center2=(79.795999, 39.599847, 0.0),
            radius=delta,
        )
    )

    # mdb.models['Model-1'].parts['Unit_cell'].NodeSet(name='P1', nodes=
    Face1_sort = sorted(Face1, key=lambda Face1: Face1.coordinates[1])
    Face1_Node = []
    mylist = list(xrange(len(Face1_sort)))
    for index in mylist:
        Face1_Node.append(Face1_sort[index].label)
        mdb.models["Model-1"].parts["Unit_cell"].SetFromNodeLabels(
            nodeLabels=(Face1_sort[index].label,), name="right_a2-a1-" + str(index)
        )

    return mylist


# -------------------------------------------------------
# The loop
# -------------------------------------------------------

Step_size = 0.02
Frame_Stiffness = (12 * E * Plate_thickess * (thickness_main) ** 3) / (
    12 * (length) ** 3 * (math.cos(angle * math.pi / 180)) ** 2
)
Number_K = 50.0
K_range = np.arange(
    Frame_Stiffness * 25.0 / Number_K, Frame_Stiffness * 1.3, Frame_Stiffness / Number_K
)
Negative_Stiffness = 0.0
mylist = Model_without_Bloch_BC(Negative_Stiffness)
for Negative_Stiffness in K_range:
    directory = "/share/amg/lsalaris/Brelium_Zone/Bloch_Zone_2D_6fold_sym_Triangular/k_"
    New_folder = directory + str(Negative_Stiffness)
    os.makedirs(New_folder)
    os.chdir(New_folder)
    # --------------------------------------------------------------
    # Delete the negative Connector
    # --------------------------------------------------------------
    del mdb.models["Model-1"].sections["Negative_conector"]
    # --------------------------------------------------------------
    # Define Connector
    # --------------------------------------------------------------
    # in kg/mm.s^2 (it is force)
    mdb.models["Model-1"].ConnectorSection(
        name="Negative_conector", translationalType=AXIAL
    )
    stop_0 = connectorBehavior.ConnectorElasticity(
        behavior=NONLINEAR,
        coupling=UNCOUPLED,
        components=(1,),
        table=(
            (-Negative_Stiffness * 100.0, 100.0),
            (Negative_Stiffness * 100.0, -100.0),
        ),
    )
    mdb.models["Model-1"].sections["Negative_conector"].setValues(
        behaviorOptions=(stop_0,)
    )
    # --------------------------------------------------------------
    # Define Equation connector
    # --------------------------------------------------------------
    my_steps = np.arange(0.0, 1.02, Step_size).tolist()
    for index_x in my_steps:
        for i in mylist:
            mdb.models["Model-1"].Equation(
                name="Pa1_real_i" + str(i),
                terms=(
                    (1.0, "Unit_cell-1.right_a1-" + str(i), 1, 1),
                    (
                        -math.cos(0.0 * index_x * math.pi),
                        "Unit_cell-1.Left_a1-" + str(i),
                        1,
                        1,
                    ),
                    (
                        math.sin(0.0 * index_x * math.pi),
                        "Unit_cell-2.Left_a1-" + str(i),
                        1,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa1_imag_i" + str(i),
                terms=(
                    (1.0, "Unit_cell-2.right_a1-" + str(i), 1, 1),
                    (
                        -math.sin(0.0 * index_x * math.pi),
                        "Unit_cell-1.Left_a1-" + str(i),
                        1,
                        1,
                    ),
                    (
                        -math.cos(0.0 * index_x * math.pi),
                        "Unit_cell-2.Left_a1-" + str(i),
                        1,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa1_real_j" + str(i),
                terms=(
                    (1.0, "Unit_cell-1.right_a1-" + str(i), 2, 1),
                    (
                        -math.cos(0.0 * index_x * math.pi),
                        "Unit_cell-1.Left_a1-" + str(i),
                        2,
                        1,
                    ),
                    (
                        math.sin(0.0 * index_x * math.pi),
                        "Unit_cell-2.Left_a1-" + str(i),
                        2,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa1_imag_j" + str(i),
                terms=(
                    (1.0, "Unit_cell-2.right_a1-" + str(i), 2, 1),
                    (
                        -math.sin(0.0 * index_x * math.pi),
                        "Unit_cell-1.Left_a1-" + str(i),
                        2,
                        1,
                    ),
                    (
                        -math.cos(0.0 * index_x * math.pi),
                        "Unit_cell-2.Left_a1-" + str(i),
                        2,
                        1,
                    ),
                ),
            )

            mdb.models["Model-1"].Equation(
                name="Pa2_real_i" + str(i),
                terms=(
                    (1.0, "Unit_cell-1.right_a2-" + str(i), 1, 1),
                    (
                        -math.cos(index_x * math.pi),
                        "Unit_cell-1.Left_a2-" + str(i),
                        1,
                        1,
                    ),
                    (
                        math.sin(index_x * math.pi),
                        "Unit_cell-2.Left_a2-" + str(i),
                        1,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa2_imag_i" + str(i),
                terms=(
                    (1.0, "Unit_cell-2.right_a2-" + str(i), 1, 1),
                    (
                        -math.sin(index_x * math.pi),
                        "Unit_cell-1.Left_a2-" + str(i),
                        1,
                        1,
                    ),
                    (
                        -math.cos(index_x * math.pi),
                        "Unit_cell-2.Left_a2-" + str(i),
                        1,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa2_real_j" + str(i),
                terms=(
                    (1.0, "Unit_cell-1.right_a2-" + str(i), 2, 1),
                    (
                        -math.cos(index_x * math.pi),
                        "Unit_cell-1.Left_a2-" + str(i),
                        2,
                        1,
                    ),
                    (
                        math.sin(index_x * math.pi),
                        "Unit_cell-2.Left_a2-" + str(i),
                        2,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa2_imag_j" + str(i),
                terms=(
                    (1.0, "Unit_cell-2.right_a2-" + str(i), 2, 1),
                    (
                        -math.sin(index_x * math.pi),
                        "Unit_cell-1.Left_a2-" + str(i),
                        2,
                        1,
                    ),
                    (
                        -math.cos(index_x * math.pi),
                        "Unit_cell-2.Left_a2-" + str(i),
                        2,
                        1,
                    ),
                ),
            )
        # for i in range (0, mylist[-1]):
        #     mdb.models['Model-1'].Equation(name='Pa2-a1_real_i'+str(i), terms=((1.0,'Unit_cell-1.Left_a2-a1-'+str(i),1,1),(-math.cos(index_x*math.pi),'Unit_cell-1.right_a2-a1-'+str(i),1,1),(math.sin(index_x*math.pi),'Unit_cell-2.right_a2-a1-'+str(i),1,1)))
        #     mdb.models['Model-1'].Equation(name='Pa2-a1_imag_i'+str(i), terms=((1.0,'Unit_cell-2.Left_a2-a1-'+str(i),1,1),(-math.sin(index_x*math.pi),'Unit_cell-1.right_a2-a1-'+str(i),1,1),(-math.cos(index_x*math.pi),'Unit_cell-2.right_a2-a1-'+str(i),1,1)))
        #     mdb.models['Model-1'].Equation(name='Pa2-a1_real_j'+str(i), terms=((1.0,'Unit_cell-1.Left_a2-a1-'+str(i),2,1),(-math.cos(index_x*math.pi),'Unit_cell-1.right_a2-a1-'+str(i),2,1),(math.sin(index_x*math.pi),'Unit_cell-2.right_a2-a1-'+str(i),2,1)))
        #     mdb.models['Model-1'].Equation(name='Pa2-a1_imag_j'+str(i), terms=((1.0,'Unit_cell-2.Left_a2-a1-'+str(i),2,1),(-math.sin(index_x*math.pi),'Unit_cell-1.right_a2-a1-'+str(i),2,1),(-math.cos(index_x*math.pi),'Unit_cell-2.right_a2-a1-'+str(i),2,1)))

        #######################################
        #######DEFINE JOB AND SUBMIT
        #######################################
        job_name = (
            "from_O_A"
            + str(my_steps.index(index_x))
            + "k_"
            + str(int(Negative_Stiffness))
        )
        mdb.Job(
            name=job_name,
            model="Model-1",
            description="",
            type=ANALYSIS,
            atTime=None,
            waitMinutes=0,
            waitHours=0,
            queue=None,
            memory=90,
            memoryUnits=PERCENTAGE,
            getMemoryFromAnalysis=True,
            explicitPrecision=SINGLE,
            nodalOutputPrecision=SINGLE,
            echoPrint=OFF,
            modelPrint=OFF,
            contactPrint=OFF,
            historyPrint=OFF,
            userSubroutine="",
            scratch="",
            multiprocessingMode=DEFAULT,
            numCpus=4,
            numDomains=4,
        )
        mdb.jobs[job_name].submit(consistencyChecking=OFF)
        mdb.jobs[job_name].waitForCompletion()
    # --------------------------------------------------------------
    # Define Equation connector
    # --------------------------------------------------------------
    my_steps2 = np.arange(0.0, 1.02, Step_size).tolist()
    for index_y in my_steps2:
        for i in mylist:
            mdb.models["Model-1"].Equation(
                name="Pa1_real_i" + str(i),
                terms=(
                    (1.0, "Unit_cell-1.right_a1-" + str(i), 1, 1),
                    (
                        -math.cos(2.0 * index_y * math.pi / 3.0),
                        "Unit_cell-1.Left_a1-" + str(i),
                        1,
                        1,
                    ),
                    (
                        math.sin(2.0 * index_y * math.pi / 3.0),
                        "Unit_cell-2.Left_a1-" + str(i),
                        1,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa1_imag_i" + str(i),
                terms=(
                    (1.0, "Unit_cell-2.right_a1-" + str(i), 1, 1),
                    (
                        -math.sin(2.0 * index_y * math.pi / 3.0),
                        "Unit_cell-1.Left_a1-" + str(i),
                        1,
                        1,
                    ),
                    (
                        -math.cos(2.0 * index_y * math.pi / 3.0),
                        "Unit_cell-2.Left_a1-" + str(i),
                        1,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa1_real_j" + str(i),
                terms=(
                    (1.0, "Unit_cell-1.right_a1-" + str(i), 2, 1),
                    (
                        -math.cos(2.0 * index_y * math.pi / 3.0),
                        "Unit_cell-1.Left_a1-" + str(i),
                        2,
                        1,
                    ),
                    (
                        math.sin(2.0 * index_y * math.pi / 3.0),
                        "Unit_cell-2.Left_a1-" + str(i),
                        2,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa1_imag_j" + str(i),
                terms=(
                    (1.0, "Unit_cell-2.right_a1-" + str(i), 2, 1),
                    (
                        -math.sin(2.0 * index_y * math.pi / 3.0),
                        "Unit_cell-1.Left_a1-" + str(i),
                        2,
                        1,
                    ),
                    (
                        -math.cos(2.0 * index_y * math.pi / 3.0),
                        "Unit_cell-2.Left_a1-" + str(i),
                        2,
                        1,
                    ),
                ),
            )

            mdb.models["Model-1"].Equation(
                name="Pa2_real_i" + str(i),
                terms=(
                    (1.0, "Unit_cell-1.right_a2-" + str(i), 1, 1),
                    (
                        -math.cos(math.pi + (index_y * math.pi / 3.0)),
                        "Unit_cell-1.Left_a2-" + str(i),
                        1,
                        1,
                    ),
                    (
                        math.sin(math.pi + (index_y * math.pi / 3.0)),
                        "Unit_cell-2.Left_a2-" + str(i),
                        1,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa2_imag_i" + str(i),
                terms=(
                    (1.0, "Unit_cell-2.right_a2-" + str(i), 1, 1),
                    (
                        -math.sin(math.pi + (index_y * math.pi / 3.0)),
                        "Unit_cell-1.Left_a2-" + str(i),
                        1,
                        1,
                    ),
                    (
                        -math.cos(math.pi + (index_y * math.pi / 3.0)),
                        "Unit_cell-2.Left_a2-" + str(i),
                        1,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa2_real_j" + str(i),
                terms=(
                    (1.0, "Unit_cell-1.right_a2-" + str(i), 2, 1),
                    (
                        -math.cos(math.pi + (index_y * math.pi / 3.0)),
                        "Unit_cell-1.Left_a2-" + str(i),
                        2,
                        1,
                    ),
                    (
                        math.sin(math.pi + (index_y * math.pi / 3.0)),
                        "Unit_cell-2.Left_a2-" + str(i),
                        2,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa2_imag_j" + str(i),
                terms=(
                    (1.0, "Unit_cell-2.right_a2-" + str(i), 2, 1),
                    (
                        -math.sin(math.pi + (index_y * math.pi / 3.0)),
                        "Unit_cell-1.Left_a2-" + str(i),
                        2,
                        1,
                    ),
                    (
                        -math.cos(math.pi + (index_y * math.pi / 3.0)),
                        "Unit_cell-2.Left_a2-" + str(i),
                        2,
                        1,
                    ),
                ),
            )

        # for i in range (0, mylist[-1]):
        #     mdb.models['Model-1'].Equation(name='Pa2-a1_real_i'+str(i), terms=((1.0,'Unit_cell-1.Left_a2-a1-'+str(i),1,1),(-math.cos(math.pi-(index_y*math.pi/3.0)),'Unit_cell-1.right_a2-a1-'+str(i),1,1),(math.sin(math.pi-(index_y*math.pi/3.0)),'Unit_cell-2.right_a2-a1-'+str(i),1,1)))
        #     mdb.models['Model-1'].Equation(name='Pa2-a1_imag_i'+str(i), terms=((1.0,'Unit_cell-2.Left_a2-a1-'+str(i),1,1),(-math.sin(math.pi-(index_y*math.pi/3.0)),'Unit_cell-1.right_a2-a1-'+str(i),1,1),(-math.cos(math.pi-(index_y*math.pi/3.0)),'Unit_cell-2.right_a2-a1-'+str(i),1,1)))
        #     mdb.models['Model-1'].Equation(name='Pa2-a1_real_j'+str(i), terms=((1.0,'Unit_cell-1.Left_a2-a1-'+str(i),2,1),(-math.cos(math.pi-(index_y*math.pi/3.0)),'Unit_cell-1.right_a2-a1-'+str(i),2,1),(math.sin(math.pi-(index_y*math.pi/3.0)),'Unit_cell-2.right_a2-a1-'+str(i),2,1)))
        #     mdb.models['Model-1'].Equation(name='Pa2-a1_imag_j'+str(i), terms=((1.0,'Unit_cell-2.Left_a2-a1-'+str(i),2,1),(-math.sin(math.pi-(index_y*math.pi/3.0)),'Unit_cell-1.right_a2-a1-'+str(i),2,1),(-math.cos(math.pi-(index_y*math.pi/3.0)),'Unit_cell-2.right_a2-a1-'+str(i),2,1)))

        #####################################
        #####DEFINE JOB AND SUBMIT
        #####################################
        job_name = (
            "from_A_B"
            + str(my_steps2.index(index_y))
            + "k_"
            + str(int(Negative_Stiffness))
        )
        mdb.Job(
            name=job_name,
            model="Model-1",
            description="",
            type=ANALYSIS,
            atTime=None,
            waitMinutes=0,
            waitHours=0,
            queue=None,
            memory=90,
            memoryUnits=PERCENTAGE,
            getMemoryFromAnalysis=True,
            explicitPrecision=SINGLE,
            nodalOutputPrecision=SINGLE,
            echoPrint=OFF,
            modelPrint=OFF,
            contactPrint=OFF,
            historyPrint=OFF,
            userSubroutine="",
            scratch="",
            multiprocessingMode=DEFAULT,
            numCpus=4,
            numDomains=4,
        )
        mdb.jobs[job_name].submit(consistencyChecking=OFF)
        mdb.jobs[job_name].waitForCompletion()
    # --------------------------------------------------------------
    # Define Equation connector
    # --------------------------------------------------------------
    my_steps3 = np.arange(0.0, 1.02, Step_size).tolist()
    for index_z in my_steps3:
        for i in mylist:
            mdb.models["Model-1"].Equation(
                name="Pa1_real_i" + str(i),
                terms=(
                    (1.0, "Unit_cell-1.right_a1-" + str(i), 1, 1),
                    (
                        -math.cos(2.0 * index_z * math.pi / 3.0),
                        "Unit_cell-1.Left_a1-" + str(i),
                        1,
                        1,
                    ),
                    (
                        math.sin(2.0 * index_z * math.pi / 3.0),
                        "Unit_cell-2.Left_a1-" + str(i),
                        1,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa1_imag_i" + str(i),
                terms=(
                    (1.0, "Unit_cell-2.right_a1-" + str(i), 1, 1),
                    (
                        -math.sin(2.0 * index_z * math.pi / 3.0),
                        "Unit_cell-1.Left_a1-" + str(i),
                        1,
                        1,
                    ),
                    (
                        -math.cos(2.0 * index_z * math.pi / 3.0),
                        "Unit_cell-2.Left_a1-" + str(i),
                        1,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa1_real_j" + str(i),
                terms=(
                    (1.0, "Unit_cell-1.right_a1-" + str(i), 2, 1),
                    (
                        -math.cos(2.0 * index_z * math.pi / 3.0),
                        "Unit_cell-1.Left_a1-" + str(i),
                        2,
                        1,
                    ),
                    (
                        math.sin(2.0 * index_z * math.pi / 3.0),
                        "Unit_cell-2.Left_a1-" + str(i),
                        2,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa1_imag_j" + str(i),
                terms=(
                    (1.0, "Unit_cell-2.right_a1-" + str(i), 2, 1),
                    (
                        -math.sin(2.0 * index_z * math.pi / 3.0),
                        "Unit_cell-1.Left_a1-" + str(i),
                        2,
                        1,
                    ),
                    (
                        -math.cos(2.0 * index_z * math.pi / 3.0),
                        "Unit_cell-2.Left_a1-" + str(i),
                        2,
                        1,
                    ),
                ),
            )

            mdb.models["Model-1"].Equation(
                name="Pa2_real_i" + str(i),
                terms=(
                    (1.0, "Unit_cell-1.right_a2-" + str(i), 1, 1),
                    (
                        -math.cos(4.0 * index_z * math.pi / 3.0),
                        "Unit_cell-1.Left_a2-" + str(i),
                        1,
                        1,
                    ),
                    (
                        math.sin(4.0 * index_z * math.pi / 3.0),
                        "Unit_cell-2.Left_a2-" + str(i),
                        1,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa2_imag_i" + str(i),
                terms=(
                    (1.0, "Unit_cell-2.right_a2-" + str(i), 1, 1),
                    (
                        -math.sin(4.0 * index_z * math.pi / 3.0),
                        "Unit_cell-1.Left_a2-" + str(i),
                        1,
                        1,
                    ),
                    (
                        -math.cos(4.0 * index_z * math.pi / 3.0),
                        "Unit_cell-2.Left_a2-" + str(i),
                        1,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa2_real_j" + str(i),
                terms=(
                    (1.0, "Unit_cell-1.right_a2-" + str(i), 2, 1),
                    (
                        -math.cos(4.0 * index_z * math.pi / 3.0),
                        "Unit_cell-1.Left_a2-" + str(i),
                        2,
                        1,
                    ),
                    (
                        math.sin(4.0 * index_z * math.pi / 3.0),
                        "Unit_cell-2.Left_a2-" + str(i),
                        2,
                        1,
                    ),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa2_imag_j" + str(i),
                terms=(
                    (1.0, "Unit_cell-2.right_a2-" + str(i), 2, 1),
                    (
                        -math.sin(4.0 * index_z * math.pi / 3.0),
                        "Unit_cell-1.Left_a2-" + str(i),
                        2,
                        1,
                    ),
                    (
                        -math.cos(4.0 * index_z * math.pi / 3.0),
                        "Unit_cell-2.Left_a2-" + str(i),
                        2,
                        1,
                    ),
                ),
            )

        # for i in range (0, mylist[-1]):
        #     mdb.models['Model-1'].Equation(name='Pa2-a1_real_i'+str(i), terms=((1.0,'Unit_cell-1.Left_a2-a1-'+str(i),1,1),(-math.cos(2.0*index_z*math.pi/3.0),'Unit_cell-1.right_a2-a1-'+str(i),1,1),(math.sin(2.0*index_z*math.pi/3.0),'Unit_cell-2.right_a2-a1-'+str(i),1,1)))
        #     mdb.models['Model-1'].Equation(name='Pa2-a1_imag_i'+str(i), terms=((1.0,'Unit_cell-2.Left_a2-a1-'+str(i),1,1),(-math.sin(2.0*index_z*math.pi/3.0),'Unit_cell-1.right_a2-a1-'+str(i),1,1),(-math.cos(2.0*index_z*math.pi/3.0),'Unit_cell-2.right_a2-a1-'+str(i),1,1)))
        #     mdb.models['Model-1'].Equation(name='Pa2-a1_real_j'+str(i), terms=((1.0,'Unit_cell-1.Left_a2-a1-'+str(i),2,1),(-math.cos(2.0*index_z*math.pi/3.0),'Unit_cell-1.right_a2-a1-'+str(i),2,1),(math.sin(2.0*index_z*math.pi/3.0),'Unit_cell-2.right_a2-a1-'+str(i),2,1)))
        #     mdb.models['Model-1'].Equation(name='Pa2-a1_imag_j'+str(i), terms=((1.0,'Unit_cell-2.Left_a2-a1-'+str(i),2,1),(-math.sin(2.0*index_z*math.pi/3.0),'Unit_cell-1.right_a2-a1-'+str(i),2,1),(-math.cos(2.0*index_z*math.pi/3.0),'Unit_cell-2.right_a2-a1-'+str(i),2,1)))
        #####################################
        #####DEFINE JOB AND SUBMIT
        #####################################
        job_name = (
            "from_B_O"
            + str(my_steps3.index(index_z))
            + "k_"
            + str(int(Negative_Stiffness))
        )
        mdb.Job(
            name=job_name,
            model="Model-1",
            description="",
            type=ANALYSIS,
            atTime=None,
            waitMinutes=0,
            waitHours=0,
            queue=None,
            memory=90,
            memoryUnits=PERCENTAGE,
            getMemoryFromAnalysis=True,
            explicitPrecision=SINGLE,
            nodalOutputPrecision=SINGLE,
            echoPrint=OFF,
            modelPrint=OFF,
            contactPrint=OFF,
            historyPrint=OFF,
            userSubroutine="",
            scratch="",
            multiprocessingMode=DEFAULT,
            numCpus=4,
            numDomains=4,
        )
        mdb.jobs[job_name].submit(consistencyChecking=OFF)
        mdb.jobs[job_name].waitForCompletion()
