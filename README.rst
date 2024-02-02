======================
Real Time Tunable Phoonic Metamaterials
======================

:Author: Ladan Salari Sharif
:Maintainer: Ladan Salari Sharif
:Language: Python 3.7
:Project: Negative-Stiffness Inclusions as a Platform for Real-Time Tunable Phononic Metamaterials



Installation
++++++++++++
**Windows:**

**Mac:**

Developement
++++++++++++
pip install -e .

Contributing
++++++++++++
1. Clone the repository and use `requirementstxt` to set up your virtual environment.

2. Set up the pre-commit hook:

.. code-block:: bash

    pre-commit install




Installing
++++++++++++
**Use the following line to pip install this code**

pip install git+https://github.com/ladan-salari/Real_Time_Tunable_Phononic_Metamaterials.git


Files
++++++++++++
**Create_part_Bloch_zone_1D_experimental_model.py**

    Create input files to model Blcoh wave simulation for 1D structure shown in the paper, Fig 3D.

**Create_part_Bloch_Zone_different_k_2D_Square_Cell_k_2_8E6_r=2.py**

    Create input files to model Blcoh wave simulation for 2D structures with negative stiffness materials, Fig 1. in the paper.
**Create_part_Bloch_Zone_different_k_2D_Triangular_Cell.py**

    Create input files to model Blcoh wave simulation for 2D infinite periodic structures with electromagnet, Fig 6. in the paper.

**Evalue_Extraction_1D.py**

    Eigenvalue extraction for different stiffness at different location of brillouin zone for 1D structure
**Evalue_Extraction_2D_differentk.py**

    Eigenvalue extraction for different stiffness at different location of brillouin zone for 2D structure with negative stiffness material

**Evalue_Extraction_2D_6_fold_symmerty_differentk.py**

    Eigenvalue extraction for different stiffness at different location of brillouin zone for 2D infinite periodic structures
