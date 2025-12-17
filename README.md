# VRES-T-Cell-Modelling
Brief description: Agent based modelling of T cells in new cancer treatment where T cells are activated outside the body using IL-2 secreting micro-rods. Code written in PhysiCell, C++.

To run the system save VRES under "sample_projects" in PhysiCell. Use "make-VRES-sample" command in Anaconda or similiar. Follow with command "make". Run executable name "VRES". Output will be in output folder. Can create a gif using "make gif" command. Will be put in output folder as "animation".
Save "substrate_and_cell_total_plotting.m" under matlab folder in PhysiCell. Run code after running the project in Anaconda or similiar.

# VRES
  Folder containing all sample PhysiCell sample project code.
  # Config 
    Folder contains PhysiCell_settings.xml.  
      - PhysiCell_settings.xml is where parameters for microenvironment, cells, initial conditions, boundary conditions, etc.       are set
        - Line 66 Micro environment
            - Line 67 - IL-2 set up (can set diffusion, decay rate, inital condition)
            - Line 83 - Chemokine set up (chemokine is what attracts T cells to rods)
        - Line 111 Cells Set Up 
            - Line 112 - T Cells 
                - Line 121 Death rate 
                - Line 165 Migration bias for chemokine attractant 
                - Line 171 Migration bias for chemokine activated (true)
                - Line 183 Secretion rules (IL-2 and chemokine secretion/uptake)
            - Line 229 - Rod Cells 
                - Line 300 Secretion rules (IL-2 and chemokine secretion/uptake)
        - Line 356 User Parameters (inital number of cells and parameters values)
  
  # custom_modules
    Folder contains custom.cpp and custom.h.
      - custom.cpp is where custom function for cell/ environment behaviour are made and stored. I have highlighted the two custom functions made for this project. 
        - Line 334 cell_proliferation_based_on_IL2
          Updates the T cell proliferation rate based on the amount of nearby IL-2
        - Line 334 secretion_rate_rod_cells
          Updates the secretion rate of rod cells based on the amount of time passed
      - custom.h is where the function headers live. All new functions must also be named in this file.
  
  # scripts
    empty 
  
  # Makefile
    Makes the project on Anaconda prompt.
  # main.cpp
    Main file that runs the code.

# substrate_and_cell_total_plotting.m
  Plots the total amount of T cells, rod cells, substrate and substrate diffusion over time. 
