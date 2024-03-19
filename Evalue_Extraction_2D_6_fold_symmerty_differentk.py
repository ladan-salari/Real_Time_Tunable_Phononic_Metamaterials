# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# extract Eginvalues of Modeling 2D 6-fold symmetry triangular unit cell
# with negative stiffness
# element and wings with Bloch BC
# Ladan Salari Sharif
# August 2017
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------


print("ODB extraction")



import math
import numpy as np
import os
from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *

# The number of the steps in each edge of brillouin zone
Index_number = 25

# number of requested evalue
Evalue_Number = 22

# Dimension of the bar
Plate_thickess = 4.0
length = 24.10 * 2.0
thickness_main = 1.78 * 2.0
angle = 45.0

# Young's moduls
E = 3.2e6

# Define the range of stiffness of electromagnet
Frame_Stiffness = (12 * E * Plate_thickess * (thickness_main) ** 3) / (
    12 * (length) ** 3 * (math.cos(angle * math.pi / 180)) ** 2)
Number_K = 50.0
K_range = np.arange(
    Frame_Stiffness * 54.0 / Number_K, Frame_Stiffness * 1.3, Frame_Stiffness / Number_K
)
for Negative_Stiffness in K_range:
    directory = "/share/amg/lsalaris/Brelium_Zone/Bloch_Zone_2D_6fold_sym_Triangular/k_"
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

        # This pulls out the first eigen value from the buckling step
        # and assigns it the File_eigenvalue variable
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
        # This pulls out the first eigen value from the buckling step
        # and assigns it the File_eigenvalue variable
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
        # This pulls out the first eigen value from the buckling step
        # and assigns it the File_eigenvalue variable
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
