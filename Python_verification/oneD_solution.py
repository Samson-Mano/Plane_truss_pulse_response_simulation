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


def pulse_force_displ_response(force_ampl_p0,force_start_time,force_end_time,stiff_K,mass_M,time_val):
    
    omega_n = math.sqrt(stiff_K/mass_M) # Omega n - angular natrual frequency
    T_n =  (2 * math.pi) / omega_n # Natural period
    t_d = force_end_time - force_start_time # Force period

    if(time_val>=force_start_time):
        t_at = time_val - force_start_time
        if(time_val<force_end_time):
            if(abs((t_d/T_n) - 0.5) == 0.000001):
                k_fact = (force_ampl_p0 / (2*stiff_K))
                displ = k_fact *(math.sin(omega_n*t_at) - (omega_n * t_at*math.cos(omega_n*t_at)))
                velo = k_fact * (omega_n**2) * t_at * (math.sin(omega_n*t_at))
                accl = k_fact * (omega_n**2) * (math.sin(omega_n * t_at) + (omega_n * t_at * math.cos(omega_n * t_at)))
            else:
                const1 = math.pi / (omega_n * t_d)
                const2 = 1 - const1**2
                k_fact = (force_ampl_p0 / stiff_K) *(1/const2)
                displ = k_fact * (math.sin((math.pi/t_d)*t_at) - const1 * math.sin(omega_n * t_at))
                velo = k_fact * ((math.pi/t_d)*math.cos((math.pi/t_d)*t_at) - (const1 * omega_n )* math.cos(omega_n * t_at))
                accl = k_fact * ((-1)*((math.pi/t_d)**2)*math.sin((math.pi/t_d)*t_at) + (const1 * omega_n * omega_n)* math.sin(omega_n * t_at))
        elif(time_val>force_end_time):
            if(abs((t_d/T_n) - 0.5) == 0.000001):
                k_fact =  ((force_ampl_p0 * math.pi) / (2*stiff_K))
                displ = k_fact * math.cos((omega_n * t_at) - math.pi)
                velo = k_fact * omega_n * math.sin(omega_n * t_at)
                accl = k_fact * omega_n * omega_n * math.cos(omega_n * t_at)
            else:
                const1 = math.pi / (omega_n * t_d)
                const2 = ((const1**2) - 1)
                k_fact =  (force_ampl_p0 / stiff_K)* ((2*const1) / const2) * math.cos(omega_n*t_d * 0.5)
                displ = k_fact * math.sin(omega_n * (t_at - (t_d * 0.5))) 
                velo = k_fact * omega_n * math.cos(omega_n * (t_at - (t_d * 0.5))) 
                accl = (-1) * k_fact * omega_n * omega_n * math.sin(omega_n * (t_at - (t_d * 0.5))) 

    return displ,velo,accl



def simple_harmonic_motion(mass_M, stiff_K, inl_displ, inl_velo, time_range, pulse_force_list):
    # Find the angular frequency
    omega_n = math.sqrt(stiff_K/ mass_M)
    t_n = (2 * math.pi) / omega_n

    # Generate time values for the simulation
    time_values = np.linspace(time_range[0], time_range[1], 10000)

    # Initialize arrays for displacement, velocity, and acceleration
    total_pulse_forces = np.zeros_like(time_values)
    displacement = np.zeros_like(time_values)
    velocity = np.zeros_like(time_values)
    acceleration = np.zeros_like(time_values)

    # Get the total force amplitude_at_time function with various time values for the pulse force system
    pulse_force_system = HalfSinePulseForce(pulse_force_list)
    total_pulse_forces = [pulse_force_system.get_amplitude_at_time(t) for t in time_values]

    # Find the individual pulse forces
    individual_pulse_force = np.zeros((len(pulse_force_list), time_values.size))
    # Generate random colors for each curve
    colors = list(mcolors.TABLEAU_COLORS.values())  # Use predefined Tableau colors

    for i,p_force in enumerate(pulse_force_list):
        individual_pulse_force_system = HalfSinePulseForce([pulse_force_list[i]])
        individual_pulse_force[i,:] = [individual_pulse_force_system.get_amplitude_at_time(t) for t in time_values]

    # Loop through time steps and calculate response
    for i, time_step in enumerate(time_values):
        # Reponse due to initial condition        
        inl_cond_displ_resp = (inl_displ * math.cos(omega_n * time_values[i])) + ((inl_velo / omega_n) * math.sin(omega_n * time_values[i]))
        inl_cond_velo_resp = ((-1)*inl_displ*omega_n*math.sin(omega_n * time_values[i]))+(inl_velo  * math.cos(omega_n * time_values[i]))
        inl_cond_accl_resp = ((-1)*inl_displ*omega_n*omega_n*math.cos(omega_n * time_values[i]))-(inl_velo *omega_n * math.sin(omega_n * time_values[i]))

        # Response due to pulse force
        plse_displ_resp = 0
        plse_velo_resp = 0
        plse_accl_resp = 0

        # Loop through individual pulse force
        for p_force in pulse_force_list:
            p0 = p_force[0] # get the force amplitude
            ft_start = p_force[1] # force start  
            ft_end = p_force[2] # force end
            t_d = ft_end - ft_start # force period
            
            if(time_values[i]>ft_start and time_values[i]<=ft_end):
                const1 = (p0/stiff_K)*(1/(1-(t_n/(2*t_d))**2))
                t_at = time_values[i] - ft_start
                p_0,v_0,a_0 = pulse_force_displ_response(p0,ft_start,ft_end,stiff_K,mass_M,time_values[i])
                # Add to the response
                plse_displ_resp = plse_displ_resp + p_0
                plse_velo_resp = plse_velo_resp + v_0
                plse_accl_resp = plse_accl_resp + a_0
            elif(time_values[i]>ft_end):
                const2 = (p0/stiff_K)*(((t_n/t_d)*math.cos((math.pi*t_d)/t_n))/(((t_n/(2*t_d))**2)-1))
                t_at = time_values[i] - ft_start
                p_0,v_0,a_0 = pulse_force_displ_response(p0,ft_start,ft_end,stiff_K,mass_M,time_values[i])
                # Add to the response
                plse_displ_resp = plse_displ_resp + p_0
                plse_velo_resp = plse_velo_resp + v_0
                plse_accl_resp = plse_accl_resp + a_0
    
        # Calculate acceleration at the current time step
        displacement[i] = inl_cond_displ_resp + plse_displ_resp
        velocity[i] = inl_cond_velo_resp + plse_velo_resp
        acceleration[i] = inl_cond_accl_resp + plse_accl_resp

    # Calculate the responses using numerical method
    displacement_numrl = np.zeros_like(time_values)
    velocity_numrl = np.zeros_like(time_values)
    acceleration_numrl = np.zeros_like(time_values)

    



    # Plot the results
    plt.figure(figsize=(10, 8))
    # Plot pulse forces
    plt.subplot(4, 1, 1)
    for i,p_force in enumerate(pulse_force_list):
        td1 = p_force[2] - p_force[1]
        plt.plot(time_values, individual_pulse_force[i,:], label=f'Pulse Force {i+1} td/tn = {td1/t_n}', color=colors[i])
    plt.xlabel('Time (s)')
    plt.ylabel('Pulse force')
    plt.legend()

    # Plot displacement
    plt.subplot(4, 1, 2)
    plt.plot(time_values, displacement, color='red', label='Displacement (x)')
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement')
    plt.legend()

    # Plot velocity
    plt.subplot(4, 1, 3)
    plt.plot(time_values, velocity, color='green', label='Velocity (v)')
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity')
    plt.legend()

    # Plot acceleration
    plt.subplot(4, 1, 4)
    plt.plot(time_values, acceleration, color='blue', label='Acceleration (a)')
    plt.xlabel('Time (s)')
    plt.ylabel('Acceleration')
    plt.legend()

    plt.tight_layout()
    plt.show()

# Example usage:
mass_M = 2.0        # Mass M
stiff_K = 4.0 * 2.0 * (math.pi**2)  # Stiffness K
inl_displ = 0.0 # initial displacement
inl_velo = 0.0 # initial velocity

time_range = (0, 10)  # Time range for the simulation (start and end time)

# Create an array of pulse forces with each element as (force_amplitude, start_time, end_time)
pulse_force_list = [
    (1.0, 1.0, 2.5),
    (5.0, 1.0, 5.0),
    (0.0, 5.0, 6.0)
]

simple_harmonic_motion(mass_M, stiff_K, inl_displ, inl_velo, time_range, pulse_force_list)