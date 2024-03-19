# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Create a single unit cell of structure with bloch BC.
# 1D structure
# that the stiffness changes with voltage.
# Units are mm,kg,s.
# Ladan Salari Sharif
# June 2017
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# IMPORT ABAQUS MODULUS
# -----------------------------------------------------------------------
from abaqus import *
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
length_winge = 13.41 + 18.2999 + 3.0  # internal wing
magnet_hieght = 12.7
width_magnet = 6.30
thickness_edge_magnet = 0.69
thickness_edge_magnet2 = 1.85
thickness_edge_magnet3 = 1.65
width_winge = 3.0 * 2.0
thickness_main = 1.78 * 2.0
length_extension = 25.0
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
    # Point 1->2
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
            + length * math.cos(angle * math.pi / 180),
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0
            - width_winge
            + thickness_main * math.sqrt(2.0),
        ),
    )
    # point 12->13
    s.Line(
        point1=(
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180),
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0
            - width_winge
            + thickness_main * math.sqrt(2.0),
        ),
        point2=(
            thickness_main * math.sqrt(2.0) / 2.0
            + length * math.cos(angle * math.pi / 180)
            + width_winge,
            width_winge
            + length * math.sin(angle * math.pi / 180)
            - thickness_main * math.sqrt(2.0) / 2.0
            - width_winge
            + thickness_main * math.sqrt(2.0),
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
            - thickness_main * math.sqrt(2.0) / 2.0
            - width_winge
            + thickness_main * math.sqrt(2.0),
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
    # MESH
    # --------------------------------------------------------------
    num_nodes = int(length / 0.5)  # 0.5mm is the size of the nodes
    p.PartitionFaceByAuto(
        face=p.faces.findAt(
            (0.0, 0.0, 0.0),
        )
    )
    faces1 = (
        mdb.models["Model-1"]
        .parts["Part-1"]
        .faces.getByBoundingCylinder(
            center1=(0.0, 0.0, 0.0), center2=(45.0, 0.0, 0.0), radius=45.0
        )
    )
    p.setMeshControls(elemShape=QUAD_DOMINATED, regions=faces1, technique=STRUCTURED)
    p.setElementType(
        elemTypes=(
            ElemType(elemCode=CPS4R, elemLibrary=STANDARD),
            ElemType(elemCode=CPS3, elemLibrary=STANDARD),
        ),
        regions=(
            mdb.models["Model-1"].parts["Part-1"].faces.findAt(((0.0, 0.0, 0.0),)),
        ),
    )
    p.seedPart(deviationFactor=0.1, size=0.4 * 2.0)
    p.generateMesh()
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
    a.Instance(name="Part1-2", part=p, dependent=ON, autoOffset=ON)
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
    edges2 = (
        mdb.models["Model-1"]
        .parts["Part-1"]
        .edges.getByBoundingCylinder(
            center1=(length_winge + width_magnet, magnet_hieght / 2.0, 0.0),
            center2=(length_winge + width_magnet, -magnet_hieght / 2.0, 0.0),
            radius=1e-4,
        )
    )
    mdb.models["Model-1"].parts["Part-1"].Set(edges=edges2, name="Left_Plate")
    edges3 = (
        mdb.models["Model-1"]
        .parts["Part-1"]
        .edges.getByBoundingCylinder(
            center1=(3.179847 + length_winge + width_magnet, magnet_hieght / 2.0, 0.0),
            center2=(3.179847 + length_winge + width_magnet, -magnet_hieght / 2.0, 0.0),
            radius=1e-4,
        )
    )
    mdb.models["Model-1"].parts["Part-1"].Set(edges=edges3, name="right_Plate")
    mdb.models["Model-1"].rootAssembly.ReferencePoint(
        point=(length_winge + width_magnet, 0.0, 0.0)
    )
    ref1 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints1 = (ref1[6],)
    a.Set(referencePoints=refPoints1, name="Left_Node_Real")

    mdb.models["Model-1"].rootAssembly.ReferencePoint(
        point=(3.179847 + length_winge + width_magnet, 0.0, 0.0)
    )
    ref2 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints2 = (ref2[8],)
    a.Set(referencePoints=refPoints2, name="Right_Node_Real")
    mdb.models["Model-1"].rootAssembly.ReferencePoint(
        point=(103.09163 + length_winge + width_magnet, 0.0, 0.0)
    )
    ref3 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints3 = (ref3[10],)
    a.Set(referencePoints=refPoints3, name="Left_Node_Imag")
    mdb.models["Model-1"].rootAssembly.ReferencePoint(
        point=(103.09163 + 3.179847 + length_winge + width_magnet, 0.0, 0.0)
    )
    ref4 = mdb.models["Model-1"].rootAssembly.referencePoints
    refPoints4 = (ref4[12],)
    a.Set(referencePoints=refPoints4, name="Right_Node_Imag")
    mdb.models["Model-1"].Coupling(
        name="Right-Real",
        controlPoint=a.sets["Right_Node_Real"],
        surface=a.sets["Part1-1.right_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )
    mdb.models["Model-1"].Coupling(
        name="Right-Imag",
        controlPoint=a.sets["Right_Node_Imag"],
        surface=a.sets["Part1-2.right_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )
    mdb.models["Model-1"].Coupling(
        name="Left-Real",
        controlPoint=a.sets["Left_Node_Real"],
        surface=a.sets["Part1-1.Left_Plate"],
        couplingType=KINEMATIC,
        influenceRadius=WHOLE_SURFACE,
        localCsys=None,
        u1=ON,
        u2=ON,
        ur3=ON,
    )
    mdb.models["Model-1"].Coupling(
        name="Left-Imag",
        controlPoint=a.sets["Left_Node_Imag"],
        surface=a.sets["Part1-2.Left_Plate"],
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
    ###UNFORTUNATLY EACH TIME I NEED TO FIND OUT THE DISTSNCE BETWEEN THE TWO PART AND PUT IT ISTEAD OF 62.22 LINE!
    mdb.models["Model-1"].rootAssembly.WirePolyLine(
        points=(
            (
                mdb.models["Model-1"].rootAssembly.referencePoints[6],
                mdb.models["Model-1"].rootAssembly.referencePoints[8],
            )
        ),
        meshable=OFF,
    )
    edges1 = mdb.models["Model-1"].rootAssembly.edges.getSequenceFromMask(
        mask=("[#1 ]",),
    )
    mdb.models["Model-1"].rootAssembly.Set(edges=edges1, name="wire_set_real")

    mdb.models["Model-1"].rootAssembly.WirePolyLine(
        points=(
            (
                mdb.models["Model-1"].rootAssembly.referencePoints[10],
                mdb.models["Model-1"].rootAssembly.referencePoints[12],
            )
        ),
        meshable=OFF,
    )
    edges2 = mdb.models["Model-1"].rootAssembly.edges.getSequenceFromMask(
        mask=("[#1 ]",),
    )
    mdb.models["Model-1"].rootAssembly.Set(edges=edges2, name="wire_set_Imag")
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
        table=((-Negative_Stiffness, 100.0), (Negative_Stiffness, -100.0)),
    )
    mdb.models["Model-1"].sections["Negative_conector"].setValues(
        behaviorOptions=(stop_0,)
    )
    mdb.models["Model-1"].rootAssembly.SectionAssignment(
        sectionName="Negative_conector", region=a.sets["wire_set_real"]
    )
    mdb.models["Model-1"].rootAssembly.SectionAssignment(
        sectionName="Negative_conector", region=a.sets["wire_set_Imag"]
    )
    # --------------------------------------------------------------
    # Define node set
    # --------------------------------------------------------------
    # Bottom Face
    delta = 1e-4
    Face1 = (
        mdb.models["Model-1"]
        .parts["Part-1"]
        .nodes.getByBoundingCylinder(
            center1=(
                thickness_main * math.sqrt(2.0) / 2.0
                + length * math.cos(angle * math.pi / 180),
                width_winge
                + length * math.sin(angle * math.pi / 180)
                - thickness_main * math.sqrt(2.0) / 2.0
                + thickness_main * math.sqrt(2.0)
                + length_extension,
                0.0,
            ),
            center2=(
                thickness_main * math.sqrt(2.0) / 2.0
                + length * math.cos(angle * math.pi / 180)
                + 2.0 * width_winge,
                width_winge
                + length * math.sin(angle * math.pi / 180)
                - thickness_main * math.sqrt(2.0) / 2.0
                + thickness_main * math.sqrt(2.0)
                + length_extension,
                0.0,
            ),
            radius=delta,
        )
    )
    Face1_sort = sorted(Face1, key=lambda Face1: Face1.coordinates[0])
    Face1_Node = []
    mylist = list(xrange(len(Face1_sort)))
    for index in mylist:
        Face1_Node.append(Face1_sort[index].label)
        mdb.models["Model-1"].parts["Part-1"].SetFromNodeLabels(
            nodeLabels=(Face1_sort[index].label,), name="P1-" + str(index)
        )
    # Top Face
    Face2 = (
        mdb.models["Model-1"]
        .parts["Part-1"]
        .nodes.getByBoundingCylinder(
            center1=(
                thickness_main * math.sqrt(2.0) / 2.0
                + length * math.cos(angle * math.pi / 180),
                -(
                    width_winge
                    + length * math.sin(angle * math.pi / 180)
                    - thickness_main * math.sqrt(2.0) / 2.0
                    + thickness_main * math.sqrt(2.0)
                    + length_extension
                ),
                0.0,
            ),
            center2=(
                thickness_main * math.sqrt(2.0) / 2.0
                + length * math.cos(angle * math.pi / 180)
                + 2.0 * width_winge,
                -(
                    width_winge
                    + length * math.sin(angle * math.pi / 180)
                    - thickness_main * math.sqrt(2.0) / 2.0
                    + thickness_main * math.sqrt(2.0)
                    + length_extension
                ),
                0.0,
            ),
            radius=delta,
        )
    )
    Face2_sort = sorted(Face2, key=lambda Face2: Face2.coordinates[0])
    Face2_Node = []
    mylist2 = list(xrange(len(Face2_sort)))
    for index2 in mylist2:
        Face2_Node.append(Face2_sort[index2].label)
        mdb.models["Model-1"].parts["Part-1"].SetFromNodeLabels(
            nodeLabels=(Face2_sort[index2].label,), name="P2-" + str(index2)
        )
    return mylist


# -------------------------------------------------------
# The loop
# -------------------------------------------------------
Step_size = 0.02
Frame_Stiffness = (E * Plate_thickess * (thickness_main) ** 3.0) / (
    (length) ** 3 * (math.cos(angle * math.pi / 180)) ** 2
)
Number_K = 3.0
K_range = np.arange(0.0, Frame_Stiffness * 1.5, Frame_Stiffness / Number_K)
Negative_Stiffness = 0.0
mylist = Model_without_Bloch_BC(Negative_Stiffness)
for Negative_Stiffness in K_range:
    directory = "/share/amg/lsalaris/Brelium_Zone/Experiment_Magnet_Acrylic/With_wings_length_thickness_2times_shorter_wings_thicker_wings_25mm/k_"
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
    my_steps = np.arange(0.0, 1.01, Step_size).tolist()
    for index_y in my_steps:
        for i in mylist:
            mdb.models["Model-1"].Equation(
                name="Pa2_real_a" + str(i),
                terms=(
                    (1.0, "Part1-1.P1-" + str(i), 2, 1),
                    (-math.cos(index_y * math.pi), "Part1-1.P2-" + str(i), 2, 1),
                    (math.sin(index_y * math.pi), "Part1-2.P2-" + str(i), 2, 1),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa2_imag_a" + str(i),
                terms=(
                    (1.0, "Part1-2.P1-" + str(i), 2, 1),
                    (-math.sin(index_y * math.pi), "Part1-1.P2-" + str(i), 2, 1),
                    (-math.cos(index_y * math.pi), "Part1-2.P2-" + str(i), 2, 1),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa1_real_a" + str(i),
                terms=(
                    (1.0, "Part1-1.P1-" + str(i), 1, 1),
                    (-math.cos(index_y * math.pi), "Part1-1.P2-" + str(i), 1, 1),
                    (math.sin(index_y * math.pi), "Part1-2.P2-" + str(i), 1, 1),
                ),
            )
            mdb.models["Model-1"].Equation(
                name="Pa1_imag_a" + str(i),
                terms=(
                    (1.0, "Part1-2.P1-" + str(i), 1, 1),
                    (-math.sin(index_y * math.pi), "Part1-1.P2-" + str(i), 1, 1),
                    (-math.cos(index_y * math.pi), "Part1-2.P2-" + str(i), 1, 1),
                ),
            )
        ########################################
        #########DEFINE JOB AND SUBMIT
        #########################################
        job_name = (
            "fromA_B_x"
            + str(my_steps.index(index_y))
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
