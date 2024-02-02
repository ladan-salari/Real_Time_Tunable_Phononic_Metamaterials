# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# extract Eginvalues of Modeling 2D 6-fold symmetry triangular unit cell
# with negative stiffness
# Ladan Salari Sharif
# October 2017
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------

print("ODB extraction")


# ODB extraction
from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
import math
import numpy
import os

# The number of the steps
Index_number = 25


# Dimension of the part
Material_name = "Acrylic"  # Name the materials
E = 3.2e6  # Young's Modulus (kg/mms^2)
nu = 0.3  # poission Ratio
Dens = 1.18e-6  # Density (kg/mm^3)

Plate_thickess = 5.0
Radious = 3.0
number_of_Neg_st_element = 36.0
Square = 8.0
Gap = 2.0
length = Square + 2 * Gap + 2 * Radious


# number of requested evalue
Evalue_Number = 100

# This while loop will run through each odb writing the appropriate eigen value to a text file
fraction = (math.pi * Radious**2) / (length**2)
Frame_Stiffness = (E * Plate_thickess * (1 - fraction) * 2 * math.pi) / (
    (1 + nu) + fraction * (1 - nu)
)
Stiffness_of_each_wire = Frame_Stiffness / (2 * number_of_Neg_st_element)
Negative_Stiffness = Stiffness_of_each_wire / 7.0
# Negative_Stiffness=0


directory = "/share/amg/lsalaris/Brelium_Zone/Bloch_Zone_Negative_stiffness_Materials/Square_With_2_Circle_and_Wholes/k_"
New_folder = directory + str(Negative_Stiffness)
os.chdir(New_folder)
counter_2 = 0
while counter_2 <= Index_number:
    odb = openOdb(
        "from_O_A" + str(counter_2) + "k_" + str(int(Negative_Stiffness)) + ".odb"
    )
    text_information1 = (
        "\nfrom_O_A" + str(counter_2) + "k_" + str(int(Negative_Stiffness))
    )
    file_evalue = open(
        "Eigen_Values_k_" + str(int(Negative_Stiffness)) + "-2D.txt", "a"
    )
    file_evalue.write(text_information1)

    # This pulls out the first eigen value from the buckling step and assigns it the File_eigenvalue variable
    try:
        for j in range(1, Evalue_Number):
            firstFrame = odb.steps["Frequency"].frames[j]
            Framedescription = firstFrame.description
            descriptionparts = Framedescription.split("=")
            File_e = descriptionparts[2]
            File_eigenvalue = File_e.split()[0]
            text_information = "," + str(File_eigenvalue)
            file_evalue.write(text_information)

        counter_2 = counter_2 + 1

    except:
        counter_2 = counter_2 + 1

counter_2 = 0
while counter_2 <= Index_number:
    odb = openOdb(
        "from_A_B" + str(counter_2) + "k_" + str(int(Negative_Stiffness)) + ".odb"
    )
    text_information1 = (
        "\nfrom_A_B" + str(counter_2) + "k_" + str(int(Negative_Stiffness))
    )
    file_evalue = open(
        "Eigen_Values_k_" + str(int(Negative_Stiffness)) + "-2D.txt", "a"
    )
    file_evalue.write(text_information1)
    # This pulls out the first eigen value from the buckling step and assigns it the File_eigenvalue variable
    try:
        for j in range(1, Evalue_Number):
            firstFrame = odb.steps["Frequency"].frames[j]
            Framedescription = firstFrame.description
            descriptionparts = Framedescription.split("=")
            File_e = descriptionparts[2]
            File_eigenvalue = File_e.split()[0]
            text_information = "," + str(File_eigenvalue)

            file_evalue.write(text_information)

        counter_2 = counter_2 + 1

    except:
        counter_2 = counter_2 + 1


counter_2 = Index_number
while counter_2 >= 0:
    odb = openOdb(
        "from_B_O" + str(counter_2) + "k_" + str(int(Negative_Stiffness)) + ".odb"
    )
    text_information1 = (
        "\nfrom_B_O" + str(counter_2) + "k_" + str(int(Negative_Stiffness))
    )
    file_evalue = open(
        "Eigen_Values_k_" + str(int(Negative_Stiffness)) + "-2D.txt", "a"
    )
    file_evalue.write(text_information1)
    # This pulls out the first eigen value from the buckling step and assigns it the File_eigenvalue variable
    try:
        for j in range(1, Evalue_Number):
            firstFrame = odb.steps["Frequency"].frames[j]
            Framedescription = firstFrame.description
            descriptionparts = Framedescription.split("=")
            File_e = descriptionparts[2]
            File_eigenvalue = File_e.split()[0]
            text_information = "," + str(File_eigenvalue)

            file_evalue.write(text_information)

        counter_2 = counter_2 - 1

    except:
        counter_2 = counter_2 - 1
counter_2 = 0
while counter_2 <= Index_number:
    odb = openOdb(
        "from_O_C" + str(counter_2) + "k_" + str(int(Negative_Stiffness)) + ".odb"
    )
    text_information1 = (
        "\nfrom_O_C" + str(counter_2) + "k_" + str(int(Negative_Stiffness))
    )
    file_evalue = open(
        "Eigen_Values_k_" + str(int(Negative_Stiffness)) + "-2D.txt", "a"
    )
    file_evalue.write(text_information1)
    # This pulls out the first eigen value from the buckling step and assigns it the File_eigenvalue variable
    try:
        for j in range(1, Evalue_Number):
            firstFrame = odb.steps["Frequency"].frames[j]
            Framedescription = firstFrame.description
            descriptionparts = Framedescription.split("=")
            File_e = descriptionparts[2]
            File_eigenvalue = File_e.split()[0]
            text_information = "," + str(File_eigenvalue)

            file_evalue.write(text_information)

        counter_2 = counter_2 + 1

    except:
        counter_2 = counter_2 + 1

counter_2 = 0
while counter_2 <= Index_number:
    odb = openOdb(
        "from_C_B" + str(counter_2) + "k_" + str(int(Negative_Stiffness)) + ".odb"
    )
    text_information1 = (
        "\nfrom_C_B" + str(counter_2) + "k_" + str(int(Negative_Stiffness))
    )
    file_evalue = open(
        "Eigen_Values_k_" + str(int(Negative_Stiffness)) + "-2D.txt", "a"
    )
    file_evalue.write(text_information1)
    # This pulls out the first eigen value from the buckling step and assigns it the File_eigenvalue variable
    try:
        for j in range(1, Evalue_Number):
            firstFrame = odb.steps["Frequency"].frames[j]
            Framedescription = firstFrame.description
            descriptionparts = Framedescription.split("=")
            File_e = descriptionparts[2]
            File_eigenvalue = File_e.split()[0]
            text_information = "," + str(File_eigenvalue)

            file_evalue.write(text_information)

        counter_2 = counter_2 + 1

    except:
        counter_2 = counter_2 + 1
