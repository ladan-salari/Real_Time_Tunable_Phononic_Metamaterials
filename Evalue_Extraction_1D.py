# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Plot the bandgap vs stiffness for one specific geometry
# with negative stiffness
# Ladan Salari Sharif
# June 2017
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
import math
import numpy
import decimal
import os


# find the range of stiffness for this model
E = 3.2e6
length = 24.10 * 2.0
Plate_thickess = 4.0
thickness_main = 1.78 * 2.0
angle = 45.0
Frame_Stiffness = (E * Plate_thickess * (thickness_main) ** 3.0) / (
    (length) ** 3 * (math.cos(angle * math.pi / 180)) ** 2
)
Number_K = 50.0
K_range = numpy.arange(0.0, Frame_Stiffness * 1.2, Frame_Stiffness / Number_K)
for Negative_Stiffness in K_range:
    directory = "/share/amg/lsalaris/Brelium_Zone/Experiment_Magnet_Acrylic/With_wings_length_thickness_2times_shorter_wings_thicker_wings_25mm/Smaller_Step_size/k_"
    New_folder = directory + str(Negative_Stiffness)
    os.chdir(New_folder)
    input_Eval = open(
        "Eigen_Values_k_" + str(int(Negative_Stiffness)) + "-magnet.txt", "r"
    )
    lines = input_Eval.readlines()
    list_x = []
    list_H = []
    list_G = []
    for num in range(1, 154):
        split_line = lines[num].split(",")
        G = decimal.Decimal(split_line[4])
        H = decimal.Decimal(split_line[5])
        x = H - G
        list_H.append(H)
        list_G.append(G)
        list_x.append(x)
    directory2 = "/share/amg/lsalaris/Brelium_Zone/Experiment_Magnet_Acrylic/With_wings_length_thickness_2times_shorter_wings_thicker_wings_25mm/Smaller_Step_size"
    os.chdir(directory2)
    text_information1 = "\n" + str(int(Negative_Stiffness))
    file_evalue = open("bandgap_vs_k2.txt", "a")
    file_evalue.write(text_information1)
    Bandgap = min(list_x)
    G_max = max(list_G)
    H_min = min(list_H)
    text_information = "," + str(Bandgap) + "," + str(G_max) + "," + str(H_min)
    file_evalue.write(text_information)
