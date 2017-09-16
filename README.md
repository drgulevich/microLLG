# microLLG

The program simulates evolution of the sublattice magnetization vector in ferromagnets or antiferromagnets on a square lattice by solving the Landau-Lifshitz-Gilbert equation in presence of spin-transfer torques (STT, homogeneous current distribution along x is implemented for now). Stationary solutions in absence of STT, such as skyrmions and domain walls, are found by starting from a seed solution and evolving the time-dependent system in presence of the Gilbert damping until the stationary state is reached.

At the moment the code is not well documented and lacks proper commentaries in the code.

**Examples:**

![Alt text](/lattice-deformation.png?raw=true "Skyrmion Lattice Deformation")

Deformation of a skyrmion lattice in presence of linearly polarized off-resonant pumping in Ref. D. Yudin, D. R. Gulevich, M. Titov, "Light-induced anisotropic skyrmion and stripe phases in Rashba ferromagnet", Phys. Rev. Lett., in press (2017). arXiv:1705.02261


![Alt text](/skyrmion-nucleation.gif?raw=true "Skyrmion Nucleation")

Test of model of a notch structure studied in Ref. J. Iwasaki, M. Mochizuki and N. Nagaosa, "Current-induced skyrmion dynamics in constricted geometries", Nan. Nanotech. 8, 743 (2013)

**Basic usage:**

Compile the C library (produces the dynamic library "libsky.so"):

$ make

Run iPython interactive interface:

$ ipython

Run the command in iPython

[1] run sky.py

to execute the simulation with the Python code sky.py linked to the compiled C library "libsky.so",

[2] run -i display

to display the last frame, or,

[3] run -i anim

to animate the simulation.

**Documentation:**

Documentation and commentaries are coming soon.
