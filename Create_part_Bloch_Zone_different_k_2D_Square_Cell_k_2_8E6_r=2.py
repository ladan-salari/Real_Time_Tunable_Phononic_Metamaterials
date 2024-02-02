# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Create a single unit cell of structure with bloch BC.
#  The square structure with a hole inside of it which has a negative stiffness
# that the stiffness changes with voltage.
# Units are mm,kg,s.
# Ladan Salari Sharif
# October 2017
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
import interaction
import math
import part
import mesh
import copy
import re
import decimal
import math
import os
import numpy
import numpy as np
import section
from connectorBehavior import *
from operator import itemgetter


# --------------------------------------------------------------
# MODEL PARAMETERS
# --------------------------------------------------------------
# All the UNITS are in "mm"

Material_name = "Acrylic"  # Name the materials
E = 3.2e6  # Young's Modulus (kg/mms^2)
nu = 0.3  # poission Ratio
Dens = 1.18e-6  # Density (kg/mm^3)
Plate_thickess = 5.0
Radious = 2.0
number_of_Neg_st_element = 36.0
Square = 8.0
Gap = 2.0
length = Square + 2 * Gap + 2 * Radious


def Model_without_Bloch_BC(Negative_Stiffness):

    # --------------------------------------------------------------
    # SKETCH PART
    # --------------------------------------------------------------
    s = mdb.models["Model-1"].ConstrainedSketch(name="__profile__", sheetSize=0.2)
    # Point 1->2
    s.Line(point1=(Square, 0.0), point2=(length, 0.0))
    # Point 2->3
    s.Line(point1=(length, 0.0), point2=(length, length))
    # Point 3->4
    s.Line(point1=(length, length), point2=(0.0, length))
    # Point 4->5
    s.Line(point1=(0.0, length), point2=(0.0, Square))
    # Point 5->6
    s.Line(point1=(0.0, Square), point2=(Square, Square))
    # Point 5->6
    s.Line(point1=(Square, Square), point2=(Square, 0.0))

    s.CircleByCenterPerimeter(
        center=(Square + Gap + Radious, Square / 2.0),
        point1=(Square + Gap, Square / 2.0),
    )

    s.CircleByCenterPerimeter(
        center=(Square / 2.0, Square + Gap + Radious),
        point1=(Square / 2.0, Square + Gap),
    )

    mdb.models["Model-1"].Part(
        dimensionality=TWO_D_PLANAR, name="Part-1", type=DEFORMABLE_BODY
    )
    mdb.models["Model-1"].parts["Part-1"].BaseShell(
        sketch=mdb.models["Model-1"].sketches["__profile__"]
    )
    p = mdb.models["Model-1"].parts["Part-1"]
    ###--------------------------------------------------------------
    ### DEFINE MATERIAL
    ###--------------------------------------------------------------
    mdb.models["Model-1"].Material(name=Material_name)
    mdb.models["Model-1"].materials[Material_name].Elastic(table=((E, nu),))
    mdb.models["Model-1"].materials[Material_name].Density(table=((Dens,),))
    mdb.models["Model-1"].HomogeneousSolidSection(
        name="Solid-section", material=Material_name, thickness=Plate_thickess
    )
    f = p.faces.findAt(((length, length, 0.0),))
    region = regionToolset.Region(faces=f)
    p.SectionAssignment(region=region, sectionName="Solid-section")
    ###--------------------------------------------------------------
    ##### ASSEMBLY
    #####--------------------------------------------------------------
    a = mdb.models["Model-1"].rootAssembly
    session.viewports["Viewport: 1"].setValues(displayedObject=a)
    session.viewports["Viewport: 1"].assemblyDisplay.setValues(mesh=ON)
    session.viewports["Viewport: 1"].assemblyDisplay.meshOptions.setValues(
        meshTechnique=OFF
    )
    a = mdb.models["Model-1"].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mdb.models["Model-1"].parts["Part-1"]
    a.Instance(name="Real", part=p, dependent=ON)
    a.Instance(name="Imag", part=p, dependent=ON, autoOffset=OFF)
    # translate the instance to the right possition
    a.translate(instanceList=("Imag",), vector=(1.5 * length, 0.0, 0.0))
    #####--------------------------------------------------------------
    ##### MESH
    #####--------------------------------------------------------------
    # partition the edges magnets
    p.PartitionFaceByShortestPath(
        faces=p.faces.findAt(
            (length, length, 0.0),
        ),
        point1=(Square, 0.0, 0.0),
        point2=(Square, length, 0.0),
    )
    p.PartitionFaceByShortestPath(
        faces=p.faces.findAt(
            (length, length, 0.0),
        ),
        point1=(0.0, Square, 0.0),
        point2=(length, Square, 0.0),
    )
    p.PartitionFaceByShortestPath(
        faces=p.faces.findAt(
            (0.0, Square, 0.0),
        ),
        point1=(Square / 2.0, Square, 0.0),
        point2=(Square / 2.0, length, 0.0),
    )
    p.PartitionFaceByShortestPath(
        faces=p.faces.findAt(
            (0.0, Square, 0.0),
        ),
        point1=(0.0, Square + Gap + Radious, 0.0),
        point2=(Square, Square + Gap + Radious, 0.0),
    )
    p.PartitionFaceByShortestPath(
        faces=p.faces.findAt(
            (Square * 2 / 3, Square, 0.0),
        ),
        point1=(0.0, Square + Gap + Radious, 0.0),
        point2=(Square, Square + Gap + Radious, 0.0),
    )
    p.PartitionFaceByShortestPath(
        faces=p.faces.findAt(
            (Square, 0.0, 0.0),
        ),
        point1=(Square + Gap + Radious, 0.0, 0.0),
        point2=(Square + Gap + Radious, Square, 0.0),
    )
    p.PartitionFaceByShortestPath(
        faces=p.faces.findAt(
            (Square, 0.0, 0.0),
        ),
        point1=(Square, Square / 2.0, 0.0),
        point2=(length, Square / 2.0, 0.0),
    )
    p.PartitionFaceByShortestPath(
        faces=p.faces.findAt(
            (length, 0.0, 0.0),
        ),
        point1=(Square, Square / 2.0, 0.0),
        point2=(length, Square / 2.0, 0.0),
    )
    # Create a dutumn Plane Top Circle
    center_plane = p.DatumPlaneByPrincipalPlane(
        offset=Square + Gap + Radious, principalPlane=YZPLANE
    ).id
    Z_axi = p.DatumAxisByPrincipalAxis(principalAxis=ZAXIS).id

    center_axi = p.DatumAxisByParToEdge(
        edge=p.datums[Z_axi], point=(Square + Gap + Radious, Square / 2.0, 0.0)
    ).id
    angle_steps = numpy.arange(
        math.pi / number_of_Neg_st_element,
        math.pi / 2.0,
        math.pi / number_of_Neg_st_element,
    ).tolist()
    for i in angle_steps:
        partition_plane = p.DatumPlaneByRotation(
            angle=i * 180.0 / math.pi,
            axis=p.datums[center_axi],
            plane=p.datums[center_plane],
        ).id
        p.PartitionEdgeByDatumPlane(
            datumPlane=p.datums[partition_plane],
            edges=p.edges.findAt(
                (
                    (
                        Square
                        + Gap
                        + Radious
                        - Radious * math.cos(math.pi / number_of_Neg_st_element),
                        Square / 2.0
                        + Radious * math.sin(math.pi / number_of_Neg_st_element),
                        0.0,
                    ),
                )
            ),
        )
        p.PartitionEdgeByDatumPlane(
            datumPlane=p.datums[partition_plane],
            edges=p.edges.findAt(
                (
                    (
                        Square
                        + Gap
                        + Radious
                        + Radious * math.cos(math.pi / number_of_Neg_st_element),
                        Square / 2.0
                        - Radious * math.sin(math.pi / number_of_Neg_st_element),
                        0.0,
                    ),
                )
            ),
        )

    for i in angle_steps:
        partition_plane = p.DatumPlaneByRotation(
            angle=180 - i * 180.0 / math.pi,
            axis=p.datums[center_axi],
            plane=p.datums[center_plane],
        ).id
        p.PartitionEdgeByDatumPlane(
            datumPlane=p.datums[partition_plane],
            edges=p.edges.findAt(
                (
                    (
                        Square
                        + Gap
                        + Radious
                        - Radious * math.cos(math.pi / number_of_Neg_st_element),
                        Square / 2.0
                        - Radious * math.sin(math.pi / number_of_Neg_st_element),
                        0.0,
                    ),
                )
            ),
        )
        p.PartitionEdgeByDatumPlane(
            datumPlane=p.datums[partition_plane],
            edges=p.edges.findAt(
                (
                    (
                        Square
                        + Gap
                        + Radious
                        + Radious * math.cos(math.pi / number_of_Neg_st_element),
                        Square / 2.0
                        + Radious * math.sin(math.pi / number_of_Neg_st_element),
                        0.0,
                    ),
                )
            ),
        )

    # Create a dutumn Plane Bottom Circle
    center_plane = p.DatumPlaneByPrincipalPlane(
        offset=Square / 2.0, principalPlane=YZPLANE
    ).id
    Z_axi = p.DatumAxisByPrincipalAxis(principalAxis=ZAXIS).id

    center_axi = p.DatumAxisByParToEdge(
        edge=p.datums[Z_axi], point=(Square / 2.0, Square + Gap + Radious, 0.0)
    ).id
    angle_steps = numpy.arange(
        math.pi / number_of_Neg_st_element,
        math.pi / 2.0,
        math.pi / number_of_Neg_st_element,
    ).tolist()
    for i in angle_steps:
        partition_plane = p.DatumPlaneByRotation(
            angle=i * 180.0 / math.pi,
            axis=p.datums[center_axi],
            plane=p.datums[center_plane],
        ).id
        p.PartitionEdgeByDatumPlane(
            datumPlane=p.datums[partition_plane],
            edges=p.edges.findAt(
                (
                    (
                        Square / 2.0
                        - Radious * math.cos(math.pi / number_of_Neg_st_element),
                        Square
                        + Gap
                        + Radious
                        + Radious * math.sin(math.pi / number_of_Neg_st_element),
                        0.0,
                    ),
                )
            ),
        )
        p.PartitionEdgeByDatumPlane(
            datumPlane=p.datums[partition_plane],
            edges=p.edges.findAt(
                (
                    (
                        Square / 2.0
                        + Radious * math.cos(math.pi / number_of_Neg_st_element),
                        Square
                        + Gap
                        + Radious
                        - Radious * math.sin(math.pi / number_of_Neg_st_element),
                        0.0,
                    ),
                )
            ),
        )

    for i in angle_steps:
        partition_plane = p.DatumPlaneByRotation(
            angle=180 - i * 180.0 / math.pi,
            axis=p.datums[center_axi],
            plane=p.datums[center_plane],
        ).id
        p.PartitionEdgeByDatumPlane(
            datumPlane=p.datums[partition_plane],
            edges=p.edges.findAt(
                (
                    (
                        Square / 2.0
                        - Radious * math.cos(math.pi / number_of_Neg_st_element),
                        Square
                        + Gap
                        + Radious
                        - Radious * math.sin(math.pi / number_of_Neg_st_element),
                        0.0,
                    ),
                )
            ),
        )
        p.PartitionEdgeByDatumPlane(
            datumPlane=p.datums[partition_plane],
            edges=p.edges.findAt(
                (
                    (
                        Square / 2.0
                        + Radious * math.cos(math.pi / number_of_Neg_st_element),
                        Square
                        + Gap
                        + Radious
                        + Radious * math.sin(math.pi / number_of_Neg_st_element),
                        0.0,
                    ),
                )
            ),
        )

    p.seedPart(deviationFactor=0.1, size=0.2)
    region = p.faces.getByBoundingCylinder(
        center1=(length / 2.0, 0.0, 0.0),
        center2=(length / 2.0, length, 0.0),
        radius=length,
    )
    p.setMeshControls(elemShape=TRI, regions=region, technique=STRUCTURED)
    ###    p.setMeshControls(elemShape=QUAD, regions=region, technique=STRUCTURED)
    p.seedEdgeByNumber(
        constraint=FINER,
        edges=p.edges.findAt(((0.0, Square + Gap + Radious + 1, 0.0),)),
        number=30,
    )
    p.seedEdgeByNumber(
        constraint=FINER,
        edges=p.edges.findAt(((0.0, Square + Gap + Radious - 1, 0.0),)),
        number=30,
    )
    p.seedEdgeByNumber(
        constraint=FINER,
        edges=p.edges.findAt(((length, Square + Gap + Radious, 0.0),)),
        number=60,
    )
    p.seedEdgeByNumber(
        constraint=FINER,
        edges=p.edges.findAt(((Square + Gap + Radious - 1, 0.0, 0.0),)),
        number=30,
    )
    p.seedEdgeByNumber(
        constraint=FINER,
        edges=p.edges.findAt(((Square + Gap + Radious + 1, 0.0, 0.0),)),
        number=30,
    )
    p.seedEdgeByNumber(
        constraint=FINER,
        edges=p.edges.findAt(((Square + Gap + Radious, length, 0.0),)),
        number=60,
    )
    p.seedEdgeByNumber(
        constraint=FINER,
        edges=p.edges.findAt(((Square / 2.0 + 1, length, 0.0),)),
        number=20,
    )
    p.seedEdgeByNumber(
        constraint=FINER,
        edges=p.edges.findAt(((Square / 2.0 - 1, length, 0.0),)),
        number=20,
    )
    p.seedEdgeByNumber(
        constraint=FINER,
        edges=p.edges.findAt(((Square / 2.0 + 1, Square, 0.0),)),
        number=20,
    )
    p.seedEdgeByNumber(
        constraint=FINER,
        edges=p.edges.findAt(((Square / 2.0 - 1, Square, 0.0),)),
        number=20,
    )
    p.seedEdgeByNumber(
        constraint=FINER,
        edges=p.edges.findAt(((length, Square / 2.0 - 1.0, 0.0),)),
        number=20,
    )
    p.seedEdgeByNumber(
        constraint=FINER,
        edges=p.edges.findAt(((length, Square / 2.0 + 1.0, 0.0),)),
        number=20,
    )
    p.seedEdgeByNumber(
        constraint=FINER,
        edges=p.edges.findAt(((Square, Square / 2.0 - 1.0, 0.0),)),
        number=20,
    )
    p.seedEdgeByNumber(
        constraint=FINER,
        edges=p.edges.findAt(((Square, Square / 2.0 + 1.0, 0.0),)),
        number=20,
    )
    p.setElementType(
        elemTypes=(
            ElemType(elemCode=CPS3, elemLibrary=STANDARD),
            ElemType(elemCode=CPS3, elemLibrary=STANDARD),
        ),
        regions=(p.faces.findAt(((length, length, 0.0),)),),
    )
    p.generateMesh()
    # --------------------------------------------------------------
    # step
    # --------------------------------------------------------------
    mdb.models["Model-1"].FrequencyStep(
        name="Frequency",
        previous="Initial",
        eigensolver=LANCZOS,
        numEigen=100,
        blockSize=DEFAULT,
    )
    ###--------------------------------------------------------------
    ## Define Constrain between edges and middle node and Ref nodes
    ##--------------------------------------------------------------
    # Top Vertex Real Part Top Circle
    angle_steps = numpy.arange(
        0.0, math.pi, math.pi / number_of_Neg_st_element
    ).tolist()
    vertex_1_Top_C = []
    for angle in angle_steps:
        vertex_1_Top_C.append(
            a.instances["Real"].vertices.findAt(
                (
                    Square / 2.0 + Radious * math.cos(angle),
                    Square + Gap + Radious + Radious * math.sin(angle),
                    0.0,
                ),
            )
        )
    # Bottom Vertex Real Part
    angle_steps = numpy.arange(
        math.pi / number_of_Neg_st_element,
        math.pi + math.pi / number_of_Neg_st_element,
        math.pi / number_of_Neg_st_element,
    ).tolist()
    vertex_2_Top_C = []
    for angle in angle_steps:
        vertex_2_Top_C.append(
            a.instances["Real"].vertices.findAt(
                (
                    Square / 2.0 + Radious * math.cos(angle),
                    Square + Gap + Radious - Radious * math.sin(angle),
                    0.0,
                ),
            )
        )

    # Top Vertex Imag Part Top Circle
    angle_steps = numpy.arange(
        0.0, math.pi, math.pi / number_of_Neg_st_element
    ).tolist()
    vertex_3_Top_C = []
    for angle in angle_steps:
        vertex_3_Top_C.append(
            a.instances["Imag"].vertices.findAt(
                (
                    Square / 2.0 + 1.5 * length + Radious * math.cos(angle),
                    Square + Gap + Radious + Radious * math.sin(angle),
                    0.0,
                ),
            )
        )
    # Bottom Vertex Imag Part
    angle_steps = numpy.arange(
        math.pi / number_of_Neg_st_element,
        math.pi + math.pi / number_of_Neg_st_element,
        math.pi / number_of_Neg_st_element,
    ).tolist()
    vertex_4_Top_C = []
    for angle in angle_steps:
        vertex_4_Top_C.append(
            a.instances["Imag"].vertices.findAt(
                (
                    Square / 2.0 + 1.5 * length + Radious * math.cos(angle),
                    Square + Gap + Radious - Radious * math.sin(angle),
                    0.0,
                ),
            )
        )

    # Top Vertex Real Part Bottom Circle
    angle_steps = numpy.arange(
        0.0, math.pi, math.pi / number_of_Neg_st_element
    ).tolist()
    vertex_1_Bottom_C = []
    for angle in angle_steps:
        vertex_1_Bottom_C.append(
            a.instances["Real"].vertices.findAt(
                (
                    Square + Gap + Radious + Radious * math.cos(angle),
                    Square / 2.0 + Radious * math.sin(angle),
                    0.0,
                ),
            )
        )
    # Bottom Vertex Real Part
    angle_steps = numpy.arange(
        math.pi / number_of_Neg_st_element,
        math.pi + math.pi / number_of_Neg_st_element,
        math.pi / number_of_Neg_st_element,
    ).tolist()
    vertex_2_Bottom_C = []
    for angle in angle_steps:
        vertex_2_Bottom_C.append(
            a.instances["Real"].vertices.findAt(
                (
                    Square + Gap + Radious + Radious * math.cos(angle),
                    Square / 2.0 - Radious * math.sin(angle),
                    0.0,
                ),
            )
        )

    # Top Vertex Imag Part Bottom Circle
    angle_steps = numpy.arange(
        0.0, math.pi, math.pi / number_of_Neg_st_element
    ).tolist()
    vertex_3_Bottom_C = []
    for angle in angle_steps:
        vertex_3_Bottom_C.append(
            a.instances["Imag"].vertices.findAt(
                (
                    Square + Gap + Radious + 1.5 * length + Radious * math.cos(angle),
                    Square / 2.0 + Radious * math.sin(angle),
                    0.0,
                ),
            )
        )
    # Bottom Vertex Imag Part
    angle_steps = numpy.arange(
        math.pi / number_of_Neg_st_element,
        math.pi + math.pi / number_of_Neg_st_element,
        math.pi / number_of_Neg_st_element,
    ).tolist()
    vertex_4_Bottom_C = []
    for angle in angle_steps:
        vertex_4_Bottom_C.append(
            a.instances["Imag"].vertices.findAt(
                (
                    Square + Gap + Radious + 1.5 * length + Radious * math.cos(angle),
                    Square / 2.0 - Radious * math.sin(angle),
                    0.0,
                ),
            )
        )
    ##
    #########--------------------------------------------------------------
    ######### Define Wire
    #########--------------------------------------------------------------
    # Define wires between points in Real Part
    for i in range(0, len(angle_steps)):

        mdb.models["Model-1"].rootAssembly.WirePolyLine(
            points=((vertex_1_Top_C[i], vertex_2_Top_C[len(angle_steps) - i - 1])),
            meshable=OFF,
        )
        edges1 = mdb.models["Model-1"].rootAssembly.edges.getSequenceFromMask(
            mask=("[#1 ]",),
        )
        mdb.models["Model-1"].rootAssembly.Set(edges=edges1, name="Wire_Real" + str(i))
    for i in range(0, len(angle_steps)):

        mdb.models["Model-1"].rootAssembly.WirePolyLine(
            points=(
                (vertex_1_Bottom_C[i], vertex_2_Bottom_C[len(angle_steps) - i - 1])
            ),
            meshable=OFF,
        )
        edges1 = mdb.models["Model-1"].rootAssembly.edges.getSequenceFromMask(
            mask=("[#1 ]",),
        )
        mdb.models["Model-1"].rootAssembly.Set(
            edges=edges1, name="Wire_Real" + str(i + len(angle_steps))
        )
    for i in range(0, len(angle_steps)):
        mdb.models["Model-1"].rootAssembly.WirePolyLine(
            points=((vertex_3_Top_C[i], vertex_4_Top_C[len(angle_steps) - i - 1])),
            meshable=OFF,
        )
        edges1 = mdb.models["Model-1"].rootAssembly.edges.getSequenceFromMask(
            mask=("[#1 ]",),
        )
        mdb.models["Model-1"].rootAssembly.Set(edges=edges1, name="Wire_Imag" + str(i))
    for i in range(0, len(angle_steps)):
        mdb.models["Model-1"].rootAssembly.WirePolyLine(
            points=(
                (vertex_3_Bottom_C[i], vertex_4_Bottom_C[len(angle_steps) - i - 1])
            ),
            meshable=OFF,
        )
        edges1 = mdb.models["Model-1"].rootAssembly.edges.getSequenceFromMask(
            mask=("[#1 ]",),
        )
        mdb.models["Model-1"].rootAssembly.Set(
            edges=edges1, name="Wire_Imag" + str(i + len(angle_steps))
        )

    ####
    #########--------------------------------------------------------------
    ######### Define Connector
    #########--------------------------------------------------------------
    mdb.models["Model-1"].ConnectorSection(
        name="Negative_conector", translationalType=AXIAL
    )
    stop_0 = connectorBehavior.ConnectorElasticity(
        behavior=NONLINEAR,
        coupling=UNCOUPLED,
        components=(1,),
        table=((-Negative_Stiffness * 100, 100.0), (Negative_Stiffness * 100, -100.0)),
    )
    mdb.models["Model-1"].sections["Negative_conector"].setValues(
        behaviorOptions=(stop_0,)
    )
    for i in range(0, len(angle_steps)):
        mdb.models["Model-1"].rootAssembly.SectionAssignment(
            sectionName="Negative_conector", region=a.sets["Wire_Real" + str(i)]
        )
        mdb.models["Model-1"].rootAssembly.SectionAssignment(
            sectionName="Negative_conector",
            region=a.sets["Wire_Real" + str(i + len(angle_steps))],
        )
    for i in range(0, len(angle_steps)):
        mdb.models["Model-1"].rootAssembly.SectionAssignment(
            sectionName="Negative_conector",
            region=a.sets["Wire_Imag" + str(i + len(angle_steps))],
        )
        mdb.models["Model-1"].rootAssembly.SectionAssignment(
            sectionName="Negative_conector", region=a.sets["Wire_Imag" + str(i)]
        )
    ##
    #####--------------------------------------------------------------
    ##### Define node set
    #####--------------------------------------------------------------
    # Right Face in a1 direction
    delta = 1e-4
    Face1 = (
        mdb.models["Model-1"]
        .parts["Part-1"]
        .nodes.getByBoundingCylinder(
            center1=(length, Square, 0.0), center2=(length, length, 0.0), radius=delta
        )
    )
    Face1_sort = sorted(Face1, key=lambda Face1: Face1.coordinates[1])
    Face1_Node = []
    mylist_1 = list(xrange(len(Face1_sort)))
    for index in mylist_1:
        Face1_Node.append(Face1_sort[index].label)
        mdb.models["Model-1"].parts["Part-1"].SetFromNodeLabels(
            nodeLabels=(Face1_sort[index].label,), name="Right-" + str(index)
        )

    # Left Face in a1 direction
    delta = 1e-4
    Face1 = (
        mdb.models["Model-1"]
        .parts["Part-1"]
        .nodes.getByBoundingCylinder(
            center1=(0.0, Square, 0.0), center2=(0.0, length, 0.0), radius=delta
        )
    )

    # mdb.models['Model-1'].parts['Unit_cell'].NodeSet(name='P1', nodes=
    Face1_sort = sorted(Face1, key=lambda Face1: Face1.coordinates[1])
    Face1_Node = []
    mylist_2 = list(xrange(len(Face1_sort)))
    for index in mylist_2:
        Face1_Node.append(Face1_sort[index].label)
        mdb.models["Model-1"].parts["Part-1"].SetFromNodeLabels(
            nodeLabels=(Face1_sort[index].label,), name="Left-" + str(index)
        )

    # Bottom in a2 direction
    delta = 1e-4
    Face1 = (
        mdb.models["Model-1"]
        .parts["Part-1"]
        .nodes.getByBoundingCylinder(
            center1=(Square, 0.0, 0.0), center2=(length, 0.0, 0.0), radius=delta
        )
    )

    Face1_sort = sorted(Face1, key=lambda Face1: Face1.coordinates[0])
    Face1_Node = []
    mylist_3 = list(xrange(len(Face1_sort)))
    for index in mylist_3:
        Face1_Node.append(Face1_sort[index].label)
        mdb.models["Model-1"].parts["Part-1"].SetFromNodeLabels(
            nodeLabels=(Face1_sort[index].label,), name="Bottom-" + str(index)
        )

    # Top Face in a2 direction
    delta = 1e-4
    Face1 = (
        mdb.models["Model-1"]
        .parts["Part-1"]
        .nodes.getByBoundingCylinder(
            center1=(Square, length, 0.0), center2=(length, length, 0.0), radius=delta
        )
    )

    Face1_sort = sorted(Face1, key=lambda Face1: Face1.coordinates[0])
    Face1_Node = []
    mylist_4 = list(xrange(len(Face1_sort)))
    for index in mylist_4:
        Face1_Node.append(Face1_sort[index].label)
        mdb.models["Model-1"].parts["Part-1"].SetFromNodeLabels(
            nodeLabels=(Face1_sort[index].label,), name="Top-" + str(index)
        )
    return mylist_1, mylist_4


###-------------------------------------------------------
#####The loop
#####-------------------------------------------------------
##
Step_size = 0.04
fraction = (math.pi * Radious**2) / (length**2)
Frame_Stiffness = (E * Plate_thickess * (1 - fraction) * 2 * math.pi) / (
    (1 + nu) + fraction * (1 - nu)
)
Stiffness_of_each_wire = Frame_Stiffness / (2 * number_of_Neg_st_element)
Negative_Stiffness = Stiffness_of_each_wire / 9.0
(mylist_1, mylist_4) = Model_without_Bloch_BC(Negative_Stiffness)
######for Negative_Stiffness in K_range:
directory = "/share/amg/lsalaris/Brelium_Zone/Bloch_Zone_Negative_stiffness_Materials/Square_With_2_Circle_and_Wholes/R=2/k_"
New_folder = directory + str(Negative_Stiffness)
os.makedirs(New_folder)
os.chdir(New_folder)
##--------------------------------------------------------------
## Delete the negative Connector
##--------------------------------------------------------------
del mdb.models["Model-1"].sections["Negative_conector"]
# --------------------------------------------------------------
## Define Connector
# --------------------------------------------------------------
# in kg/mm.s^2 (it is force)
mdb.models["Model-1"].ConnectorSection(
    name="Negative_conector", translationalType=AXIAL
)
stop_0 = connectorBehavior.ConnectorElasticity(
    behavior=NONLINEAR,
    coupling=UNCOUPLED,
    components=(1,),
    table=((-Negative_Stiffness * 100.0, 100.0), (Negative_Stiffness * 100.0, -100.0)),
)
mdb.models["Model-1"].sections["Negative_conector"].setValues(behaviorOptions=(stop_0,))
##--------------------------------------------------------------
## Define Equation connector
##--------------------------------------------------------------
my_steps = numpy.arange(0.0, 1.04, Step_size).tolist()
for index_x in my_steps:
    for i in mylist_1:
        mdb.models["Model-1"].Equation(
            name="r1_Real_i" + str(i),
            terms=(
                (1.0, "Real.Right-" + str(i), 1, 1),
                (-math.cos(index_x * math.pi), "Real.Left-" + str(i), 1, 1),
                (math.sin(index_x * math.pi), "Imag.Left-" + str(i), 1, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r1_Imag_i" + str(i),
            terms=(
                (1.0, "Imag.Right-" + str(i), 1, 1),
                (-math.sin(index_x * math.pi), "Real.Left-" + str(i), 1, 1),
                (-math.cos(index_x * math.pi), "Imag.Left-" + str(i), 1, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r1_Real_j" + str(i),
            terms=(
                (1.0, "Real.Right-" + str(i), 2, 1),
                (-math.cos(index_x * math.pi), "Real.Left-" + str(i), 2, 1),
                (math.sin(index_x * math.pi), "Imag.Left-" + str(i), 2, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r1_Imag_j" + str(i),
            terms=(
                (1.0, "Imag.Right-" + str(i), 2, 1),
                (-math.sin(index_x * math.pi), "Real.Left-" + str(i), 2, 1),
                (-math.cos(index_x * math.pi), "Imag.Left-" + str(i), 2, 1),
            ),
        )
    for i in range(0, len(mylist_4) - 1):
        mdb.models["Model-1"].Equation(
            name="r2_Real_i" + str(i),
            terms=(
                (1.0, "Real.Top-" + str(i), 1, 1),
                (-math.cos(0.0 * index_x * math.pi), "Real.Bottom-" + str(i), 1, 1),
                (math.sin(0.0 * index_x * math.pi), "Imag.Bottom-" + str(i), 1, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r2_Imag_i" + str(i),
            terms=(
                (1.0, "Imag.Top-" + str(i), 1, 1),
                (-math.sin(0.0 * index_x * math.pi), "Real.Bottom-" + str(i), 1, 1),
                (-math.cos(0.0 * index_x * math.pi), "Imag.Bottom-" + str(i), 1, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r2_Real_j" + str(i),
            terms=(
                (1.0, "Real.Top-" + str(i), 2, 1),
                (-math.cos(0.0 * index_x * math.pi), "Real.Bottom-" + str(i), 2, 1),
                (math.sin(0.0 * index_x * math.pi), "Imag.Bottom-" + str(i), 2, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r2_Imag_j" + str(i),
            terms=(
                (1.0, "Imag.Top-" + str(i), 2, 1),
                (-math.sin(0.0 * index_x * math.pi), "Real.Bottom-" + str(i), 2, 1),
                (-math.cos(0.0 * index_x * math.pi), "Imag.Bottom-" + str(i), 2, 1),
            ),
        )

    #############################################
    #############DEFINE JOB AND SUBMIT
    ############################################
    job_name = (
        "from_O_A" + str(my_steps.index(index_x)) + "k_" + str(int(Negative_Stiffness))
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
###--------------------------------------------------------------
#### Define Equation connector
###--------------------------------------------------------------
my_steps2 = numpy.arange(0.0, 1.04, Step_size).tolist()
for index_y in my_steps2:
    for i in mylist_1:
        mdb.models["Model-1"].Equation(
            name="r1_Real_i" + str(i),
            terms=(
                (1.0, "Real.Right-" + str(i), 1, 1),
                (-math.cos(math.pi), "Real.Left-" + str(i), 1, 1),
                (math.sin(math.pi), "Imag.Left-" + str(i), 1, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r1_Imag_i" + str(i),
            terms=(
                (1.0, "Imag.Right-" + str(i), 1, 1),
                (-math.sin(math.pi), "Real.Left-" + str(i), 1, 1),
                (-math.cos(math.pi), "Imag.Left-" + str(i), 1, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r1_Real_j" + str(i),
            terms=(
                (1.0, "Real.Right-" + str(i), 2, 1),
                (-math.cos(math.pi), "Real.Left-" + str(i), 2, 1),
                (math.sin(math.pi), "Imag.Left-" + str(i), 2, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r1_Imag_j" + str(i),
            terms=(
                (1.0, "Imag.Right-" + str(i), 2, 1),
                (-math.sin(math.pi), "Real.Left-" + str(i), 2, 1),
                (-math.cos(math.pi), "Imag.Left-" + str(i), 2, 1),
            ),
        )
    for i in range(0, len(mylist_4) - 1):
        mdb.models["Model-1"].Equation(
            name="r2_Real_i" + str(i),
            terms=(
                (1.0, "Real.Top-" + str(i), 1, 1),
                (-math.cos(index_y * math.pi), "Real.Bottom-" + str(i), 1, 1),
                (math.sin(index_y * math.pi), "Imag.Bottom-" + str(i), 1, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r2_Imag_i" + str(i),
            terms=(
                (1.0, "Imag.Top-" + str(i), 1, 1),
                (-math.sin(index_y * math.pi), "Real.Bottom-" + str(i), 1, 1),
                (-math.cos(index_y * math.pi), "Imag.Bottom-" + str(i), 1, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r2_Real_j" + str(i),
            terms=(
                (1.0, "Real.Top-" + str(i), 2, 1),
                (-math.cos(index_y * math.pi), "Real.Bottom-" + str(i), 2, 1),
                (math.sin(index_y * math.pi), "Imag.Bottom-" + str(i), 2, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r2_Imag_j" + str(i),
            terms=(
                (1.0, "Imag.Top-" + str(i), 2, 1),
                (-math.sin(index_y * math.pi), "Real.Bottom-" + str(i), 2, 1),
                (-math.cos(index_y * math.pi), "Imag.Bottom-" + str(i), 2, 1),
            ),
        )
        #
    #########################################
    #########DEFINE JOB AND SUBMIT
    ########################################
    job_name = (
        "from_A_B" + str(my_steps2.index(index_y)) + "k_" + str(int(Negative_Stiffness))
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
####--------------------------------------------------------------
### Define Equation connector
###--------------------------------------------------------------
my_steps3 = numpy.arange(0.0, 1.04, Step_size).tolist()
for index_z in my_steps3:
    for i in mylist_1:
        mdb.models["Model-1"].Equation(
            name="r1_Real_i" + str(i),
            terms=(
                (1.0, "Real.Right-" + str(i), 1, 1),
                (-math.cos(index_z * math.pi), "Real.Left-" + str(i), 1, 1),
                (math.sin(index_z * math.pi), "Imag.Left-" + str(i), 1, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r1_Imag_i" + str(i),
            terms=(
                (1.0, "Imag.Right-" + str(i), 1, 1),
                (-math.sin(index_z * math.pi), "Real.Left-" + str(i), 1, 1),
                (-math.cos(index_z * math.pi), "Imag.Left-" + str(i), 1, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r1_Real_j" + str(i),
            terms=(
                (1.0, "Real.Right-" + str(i), 2, 1),
                (-math.cos(index_z * math.pi), "Real.Left-" + str(i), 2, 1),
                (math.sin(index_z * math.pi), "Imag.Left-" + str(i), 2, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r1_Imag_j" + str(i),
            terms=(
                (1.0, "Imag.Right-" + str(i), 2, 1),
                (-math.sin(index_z * math.pi), "Real.Left-" + str(i), 2, 1),
                (-math.cos(index_z * math.pi), "Imag.Left-" + str(i), 2, 1),
            ),
        )
    for i in range(0, len(mylist_4) - 1):
        mdb.models["Model-1"].Equation(
            name="r2_Real_i" + str(i),
            terms=(
                (1.0, "Real.Top-" + str(i), 1, 1),
                (-math.cos(index_z * math.pi), "Real.Bottom-" + str(i), 1, 1),
                (math.sin(index_z * math.pi), "Imag.Bottom-" + str(i), 1, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r2_Imag_i" + str(i),
            terms=(
                (1.0, "Imag.Top-" + str(i), 1, 1),
                (-math.sin(index_z * math.pi), "Real.Bottom-" + str(i), 1, 1),
                (-math.cos(index_z * math.pi), "Imag.Bottom-" + str(i), 1, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r2_Real_j" + str(i),
            terms=(
                (1.0, "Real.Top-" + str(i), 2, 1),
                (-math.cos(index_z * math.pi), "Real.Bottom-" + str(i), 2, 1),
                (math.sin(index_z * math.pi), "Imag.Bottom-" + str(i), 2, 1),
            ),
        )
        mdb.models["Model-1"].Equation(
            name="r2_Imag_j" + str(i),
            terms=(
                (1.0, "Imag.Top-" + str(i), 2, 1),
                (-math.sin(index_z * math.pi), "Real.Bottom-" + str(i), 2, 1),
                (-math.cos(index_z * math.pi), "Imag.Bottom-" + str(i), 2, 1),
            ),
        )

    #############################################
    #############DEFINE JOB AND SUBMIT
    #############################################
    job_name = (
        "from_B_O" + str(my_steps3.index(index_z)) + "k_" + str(int(Negative_Stiffness))
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
