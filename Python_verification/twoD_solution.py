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
            if(abs((t_d/T_n) - 0.5) < 0.000001):
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
            if(abs((t_d/T_n) - 0.5) < 0.000001):
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

def mdof_simple_harmonic_motion_analytical(mass_M, stiff_K, inl_displ, inl_velo, time_range, pulse_force_list):
    # Find the angular frequency
    numDOF = 2 # number of defree of freedom
    # Calculate eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eig(np.linalg.inv(mass_M) @ stiff_K)

    # Sort the eigenvalues and eigenvectors in ascending order
    sorted_indices = np.argsort(eigenvalues)
    natural_frequencies = np.sqrt(np.abs(eigenvalues[sorted_indices]))
    mode_shapes = eigenvectors[:, sorted_indices]

    # Normalize the mode shapes
    normalized_mode_shapes = mode_shapes / np.abs(mode_shapes).max(axis=0)

    # Orthogonality of mass matrix
    norm_mass = normalized_mode_shapes.T @ mass_M @ normalized_mode_shapes 
    norm_stiff = normalized_mode_shapes.T @ stiff_K @ normalized_mode_shapes 

    # Generate time values for the simulation
    t_count = 1000 # time step count
    time_values = np.linspace(time_range[0], time_range[1], t_count)
    t_step = (time_range[1] - time_range[0])/ t_count # time step

    # Find the individual pulse forces
    individual_pulse_force = np.zeros((len(pulse_force_list), time_values.size))
    # Generate random colors for each curve
    colors = list(mcolors.TABLEAU_COLORS.values())  # Use predefined Tableau colors

    # Modal transform the initial displacement and initial velocity
    modal_inl_displ = np.zeros(numDOF) # modal initial displacement
    modal_inl_velo =  np.zeros(numDOF) # modal initial velocity
    omega_n =  np.zeros(numDOF) # modal angular frequency
    t_n = np.zeros(numDOF) # Time period

    # Initialize arrays for displacement, velocity, and acceleration
    total_pulse_forces = np.zeros((numDOF,t_count))
    displacement = np.zeros((numDOF,t_count))
    velocity = np.zeros((numDOF,t_count))
    acceleration = np.zeros((numDOF,t_count))

    for j in range(2):
        modal_inl_displ[j] = np.dot(normalized_mode_shapes[:,j],inl_displ)
        modal_inl_velo[j] = np.dot(normalized_mode_shapes[:,j],inl_velo)
        omega_n[j] = math.sqrt(norm_stiff[j,j]/norm_mass[j,j])
        t_n[j] = (2 * math.pi) / omega_n[j]

    # Loop through time steps and calculate response
    for i, time_step in enumerate(time_values):
        modal_displ_resp = np.zeros(numDOF) # Displacement response
        modal_velo_resp =  np.zeros(numDOF) # Velocity response
        modal_accl_resp =  np.zeros(numDOF) # Acceleration response
        for j in range(2):
            # Reponse due to initial condition        
            inl_cond_displ_resp = (modal_inl_displ[j] * math.cos(omega_n[j] * time_values[i])) + ((modal_inl_velo[j] / omega_n[j]) * math.sin(omega_n[j] * time_values[i]))
            inl_cond_velo_resp = ((-1)*modal_inl_displ[j]*omega_n[j]*math.sin(omega_n[j] * time_values[i]))+(modal_inl_velo[j]  * math.cos(omega_n[j] * time_values[i]))
            inl_cond_accl_resp = ((-1)*modal_inl_displ[j]*omega_n[j]*omega_n[j]*math.cos(omega_n[j] * time_values[i]))-(modal_inl_velo[j] *omega_n[j] * math.sin(omega_n[j] * time_values[i]))

            # Response due to pulse force
            plse_displ_resp = 0
            plse_velo_resp = 0
            plse_accl_resp = 0

            # Loop through individual nodes (pulse force)
            for k in range(2):
                # Loop through individual pulse force at a particular node
                for p_force_nd in pulse_force_list[k]:
                    for p_force in p_force_nd:
                        p0 = normalized_mode_shapes[k,j]  * p_force[0] # get the force amplitude
                        ft_start = p_force[1] # force start  
                        ft_end = p_force[2] # force end
                        t_d = ft_end - ft_start # force period
                    
                        if(time_values[i]>ft_start and time_values[i]<=ft_end):
                            p_0,v_0,a_0 = pulse_force_displ_response(p0,ft_start,ft_end,norm_stiff[k,k],norm_mass[k,k],time_values[i])
                            # Add to the response
                            plse_displ_resp = plse_displ_resp + p_0
                            plse_velo_resp = plse_velo_resp + v_0
                            plse_accl_resp = plse_accl_resp + a_0
                        elif(time_values[i]>ft_end):
                            p_0,v_0,a_0 = pulse_force_displ_response(p0,ft_start,ft_end,norm_stiff[k,k],norm_mass[k,k],time_values[i])
                            # Add to the response
                            plse_displ_resp = plse_displ_resp + p_0
                            plse_velo_resp = plse_velo_resp + v_0
                            plse_accl_resp = plse_accl_resp + a_0

            # Append to the displacement, velocity and acceleration response
            modal_displ_resp[j] = inl_cond_displ_resp + plse_displ_resp
            modal_velo_resp[j] =  inl_cond_velo_resp + plse_velo_resp
            modal_accl_resp[j] =  inl_cond_accl_resp + plse_accl_resp
        #_________________________________________________________________________
        for j in range(2):
            displacement[j,i] = np.dot( normalized_mode_shapes[j,:], modal_displ_resp)
            velocity[j,i] = np.dot( normalized_mode_shapes[j,:] ,modal_velo_resp)
            acceleration[j,i] = np.dot( normalized_mode_shapes[j,:] , modal_accl_resp)
        #_________________________________________________________________________


    # Print the results
    print("Normalized mass:")
    print(norm_mass)
    print("Normalized stiffness:")    
    print(norm_stiff)

    # Print the results
    print("Natural Frequencies (Hz):")
    print(natural_frequencies)

    print("\nNormalized Mode Shapes:")
    print(normalized_mode_shapes)


    # Plot the results
    plt.figure(figsize=(10, 8))
    # Plot pulse forces
    plt.subplot(4, 1, 1)
    # plt.plot(time_values, total_pulse_forces, color='darkred', label='Total pulse force (f)')
    for j in range(2):
        for p_force_nd in pulse_force_list[j]:
            for p_force in p_force_nd:
                individual_pulse_force_system = HalfSinePulseForce([(p_force[0], p_force[1], p_force[2])])
                individual_pulse_force = [individual_pulse_force_system.get_amplitude_at_time(t) for t in time_values]
                td1 = p_force[2] - p_force[1]
                plt.plot(time_values, individual_pulse_force, label=f'Node {j} Pulse Force {i+1} td/tn = {td1/t_n[j]}', color=colors[j])
    plt.xlabel('Time (s)')
    plt.ylabel('Pulse force')
    plt.legend()

    # Plot displacement
    plt.subplot(4, 1, 2)
    plt.plot(time_values, displacement[0,:], color='darkred', label='Displacement node 1 (x)')
    plt.plot(time_values, displacement[1,:], color='red', label='Displacement node 2 (x)')
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement')
    plt.legend()

    # Plot velocity
    plt.subplot(4, 1, 3)
    plt.plot(time_values, velocity[0,:], color='darkgreen', label='Velocity node 1 (v)')
    plt.plot(time_values, velocity[1,:], color='green', label='Velocity node 2 (v)')
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity')
    plt.legend()

    # Plot acceleration
    plt.subplot(4, 1, 4)
    plt.plot(time_values, acceleration[0,:], color='darkblue', label='Acceleration node 1 (a)')
    plt.plot(time_values, acceleration[1,:], color='blue', label='Acceleration node 2 (a)')
    plt.xlabel('Time (s)')
    plt.ylabel('Acceleration')
    plt.legend()

    plt.tight_layout()
    plt.show()


def mdof_simple_harmonic_solution_numerical(mass_M, stiff_K, inl_displ, inl_velo, time_range, pulse_force_list):
    # Time range
    t_start, t_end = time_range
    time_steps = 1000
    delta_t = (t_end - t_start) / time_steps

    # Function to calculate the pulse force at a given time
    def get_pulse_force(pulse_force_list, t):
        force = 0.0
        for f in pulse_force_list:
            if f[1] <= t <= f[2]:
                force += f[0]
        return force

    # Central difference time integration
    def central_difference_integration(displ_n, velo_n, t, pulse_force_list):
        accel_n = np.linalg.solve(mass_M, get_pulse_force(pulse_force_list, t) - stiff_K @ displ_n)
        displ_n_plus_1 = displ_n + velo_n * delta_t + 0.5 * accel_n * delta_t ** 2
        accel_n_plus_1 = np.linalg.solve(mass_M, get_pulse_force(pulse_force_list, t + delta_t) - stiff_K @ displ_n_plus_1)
        velo_n_plus_1 = velo_n + 0.5 * (accel_n + accel_n_plus_1) * delta_t
        return displ_n_plus_1, velo_n_plus_1, accel_n

    # Arrays to store the results
    displ_results = np.zeros((2, time_steps + 1))
    velo_results = np.zeros((2, time_steps + 1))
    accel_results = np.zeros((2, time_steps + 1))

    # Set initial conditions
    print(inl_displ)

    displ_results[:, 0] = np.array(inl_displ).flatten()
    velo_results[:, 0] =  np.array(inl_velo).flatten()
    accel_results[:, 0] = np.linalg.solve(mass_M, get_pulse_force(pulse_force_list[0], 0.0) - stiff_K @ displ_results[:, 0])

    # Perform time integration to calculate the response
    for i in range(time_steps):
        t = t_start + i * delta_t
        pulse_force_list_i = pulse_force_list[i % 2]  # Access the pulse force list for the current time step
        displ_results[:, i + 1], velo_results[:, i + 1], accel_results[:, i + 1] = central_difference_integration(
            displ_results[:, i], velo_results[:, i], t, pulse_force_list_i
        )

    ## Calculate acceleration from displacements and velocities
    #for i in range(time_steps + 1):
        #t = t_start + i * delta_t
        #pulse_force_list_i = pulse_force_list[i % 2]  # Access the pulse force list for the current time step
        #accel_results[:, i] = np.linalg.solve(mass_M, get_pulse_force(pulse_force_list_i, t) - stiff_K @ displ_results[:, i])

    # Time array for plotting
    time_array = np.linspace(t_start, t_end, time_steps + 1)

    # Plotting the results
    plt.figure(figsize=(10, 6))
    plt.subplot(3, 1, 1)
    plt.plot(time_array, displ_results[0],  color='darkred', label='Displacement node 1 (x)')
    plt.plot(time_array, displ_results[1],  color='red', label='Displacement node 2 (x)')
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement')
    plt.legend()

    plt.subplot(3, 1, 2)
    plt.plot(time_array, velo_results[0], color='darkgreen', label='Velocity node 1 (x)')
    plt.plot(time_array, velo_results[1],  color='green', label='Velocity node 2 (x)')
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity')
    plt.legend()

    plt.subplot(3, 1, 3)
    plt.plot(time_array, accel_results[0], color='darkblue', label='Acceleration node 1 (x)')
    plt.plot(time_array, accel_results[1],  color='blue', label='Acceleration node 2 (x)')
    plt.xlabel('Time (s)')
    plt.ylabel('Acceleration')
    plt.legend()

    plt.tight_layout()
    plt.show()


# Usage:
# Mass
mass_M1 = 20.0 # Mass M1
mass_M2 = 20.0 # Mass M2
mass_M = [[mass_M1,0.0],
          [0.0,mass_M2]]

# Stiffness
stiff_K1 = 4.0 * 2.0 * (math.pi**2)  # Stiffness K1
stiff_K2 = 4.0 * 2.0 * (math.pi**2)  # Stiffness K2
stiff_K3 = 4.0 * 2.0 * (math.pi**2)  # Stiffness K2
stiff_K = [[(stiff_K1+stiff_K2),(-stiff_K2)],
           [(-stiff_K2),(stiff_K2+stiff_K3)]]


# Initial condition
inl_displ_1 = 0.0 # initial displacement Node 1
inl_velo_1 = 0.0 # initial velocity Node 1
#_______________________________________________
inl_displ_2 = 0.0 # initial displacement Node 2
inl_velo_2 = 0.0 # initial velocity Node 2

inl_displ = [[inl_displ_1],
             [inl_displ_2]]

inl_velo = [[inl_velo_1],
            [inl_velo_2]]

# Time range
time_range = (0, 10)  # Time range for the simulation (start and end time)

# Create an array of pulse forces with each element as (force_amplitude, start_time, end_time)
# Pulse force list 1
pulse_force_list_1 = [
    (10.0, 1.0, 5.0),
    (0.0, 1.0, 3.0),
    (0.0, 2.0, 4.5)
]

# Pulse force list 2
pulse_force_list_2 = [
    (0.0, 1.0, 2.5),
    (0.0, 1.0, 3.0),
    (0.0, 2.0, 4.5)
]

pulse_force_list = [[pulse_force_list_1],
                    [pulse_force_list_2]]
#____________________________________________________________________________________________________________

mdof_simple_harmonic_motion_analytical(mass_M, stiff_K, inl_displ, inl_velo, time_range, pulse_force_list)
mdof_simple_harmonic_solution_numerical(mass_M, stiff_K, inl_displ, inl_velo, time_range, [pulse_force_list_1,pulse_force_list_2])

