import numpy as np
import math as math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

class HalfSinePulseForce:
    def __init__(self, pulse_forces):
        self.__pulse_forces = pulse_forces

    def get_amplitude_at_time(self, t):
        total_force = 0.0
        for force_amplitude, start_time, end_time in self.__pulse_forces:
            if start_time <= t <= end_time:
                total_force += force_amplitude * np.sin(np.pi * (t - start_time) / (end_time - start_time))
        return total_force



def mdof_simple_harmonic_motion(mass_M, stiff_K, inl_displ, inl_velo, time_range, pulse_force_list):
    # Find the angular frequency


# Usage:
# Mass
mass_M1 = 2.0 # Mass M1
mass_M2 = 5.0 # Mass M2
mass_M = [[mass_M1,0.0]
          [0.0,mass_M2]]

# Stiffness
stiff_K1 = 4.0 * 2.0 * (math.pi**2)  # Stiffness K1
stiff_K2 = 4.0 * 5.0 * (math.pi**2)  # Stiffness K2
stiff_K3 = 4.0 * 1.0 * (math.pi**2)  # Stiffness K2
stiff_K = [[stiff_K1+stiff_K2,-stiff_K2]
           [-stiff_K2,stiff_K2+stiff_K3]]


# Initial condition
inl_displ_1 = 0.0 # initial displacement Node 1
inl_velo_1 = 0.0 # initial velocity Node 1
#_______________________________________________
inl_displ_2 = 0.0 # initial displacement Node 2
inl_velo_2 = 0.0 # initial velocity Node 2

inl_displ = [[inl_displ_1]
             [inl_displ_2]]

inl_velo = [[inl_velo_1]
            [inl_velo_2]]

# Time range
time_range = (0, 10)  # Time range for the simulation (start and end time)

# Create an array of pulse forces with each element as (force_amplitude, start_time, end_time)
# Pulse force list 1
pulse_force_list_1 = [
    (0.0, 1.0, 2.5),
    (100, 1.0, 3.0),
    (200, 2.0, 4.5)
]

# Pulse force list 2
pulse_force_list_2 = [
    (0.0, 1.0, 2.5),
    (100, 1.0, 3.0),
    (200, 2.0, 4.5)
]

pulse_force_list = [[pulse_force_list_1]
                    [pulse_force_list_2]]
#______________________________________________________________________________
