from landlab.components import FlowAccumulator, LinearDiffuser
from landlab import FIXED_VALUE_BOUNDARY
from landlab import RasterModelGrid
from scipy import optimize
from mpmath import *
import numpy as np
import scipy

# Function for computing the erosion term using D_infinity flow direction method
def DinfEroder(mg, dt, K_sp, m_sp, n_sp):
    '''
    Function arguments:
    mg  :                    Raster grid used for the simulation
    dt  :                                  Time-step to be taken
    K_sp:                                    Erosion coefficient
    m_sp: Exponent of specific drainage area in the erosion term
    n_sp:            Exponent of local slope in the erosion term
    '''
    
    # Initializing local arrays for the function block
    nodes      =                                                      np.arange(mg.number_of_nodes)
    receiver   =                                                  mg.at_node['flow__receiver_node']
    A          =                                               np.copy(mg.at_node['drainage_area'])
    proportion =                                           mg.at_node['flow__receiver_proportions']
    ele        =                                      np.copy(mg.at_node['topographic__elevation'])
    neighbor   = np.append(mg.diagonal_adjacent_nodes_at_node, mg.adjacent_nodes_at_node, axis = 1)
    donor      =                                                            np.zeros_like(neighbor)
    
    # Listing the receivers of each node's neighbors
    receiver_of_neighbor        =                                         receiver[neighbor]
    receiver_of_neighbor[:,:,0] = np.where(neighbor == -1 , -1, receiver_of_neighbor[:,:,0])
    receiver_of_neighbor[:,:,1] = np.where(neighbor == -1 , -1, receiver_of_neighbor[:,:,1])
    
    # Determining the donors for each node in the domain
    for nei in range(8):
        donor[:, nei] = np.where((receiver_of_neighbor[:, nei, 0] == nodes) |
                                 (receiver_of_neighbor[:, nei, 1] == nodes) , neighbor[:, nei], -1)

    # Initializing the queue by adding sink nodes
    queue           = list(np.arange(mg.number_of_nodes)[mg.at_node['flow__sink_flag'] == 1])
    
    # Marking the sink nodes as processed
    processed_nodes        =   np.zeros_like(nodes, dtype=bool)
    processed_nodes[queue] =                               True
    
    # Traversing the remaining network
    while len(queue) > 0:
        
        # Popping the first element of the queue
        node = queue.pop(0)
        
        # Visiting its donors
        for nei in donor[node, donor[node] >= 0]:
            
            # Processing the donor if all the receivers of that donor are processesd
            if processed_nodes[nei] == False and np.all(processed_nodes[receiver[nei][receiver[nei] >= 0]]):
                processed_nodes[nei] = True
                
                # Finding the slope using the elevation of neighbors in cardinal and diagonal directions
                h_0 , h_1 = None, None
                for res, pro in zip(receiver[nei], proportion[nei]):
                    if res != -1:
                        if res // Nc == nei // Nc or res % Nc == nei % Nc:
                            h_0 = ele[res]
                        else:
                            h_1 = ele[res]
                h = ele[nei]
                width = dx
                if h_0 is not None and h_1 is not None:
                    f = lambda h_next: h_next - h + K_sp * (A[nei] / width)**m_sp *\
                                ((h_0 - h_1)**2 + (h_next - h_0)**2)**(n_sp / 2.) * dt / dx**n_sp
                    min_h = min(h_0, h_1)

                elif h_0 is not None:
                    f = lambda h_next: h_next - h + K_sp * (A[nei] / width)**m_sp *\
                                ((h_next - h_0)**2)**(n_sp / 2.) * dt / dx**n_sp
                    min_h = h_0
                elif h_1 is not None:
                    f = lambda h_next: h_next - h + K_sp * (A[nei] / width)**m_sp *\
                                (h_next - h_1)**n_sp * dt / (2**.5 * dx)**n_sp

                    min_h = h_1
                    
                # Updating the elevation implicitly using Brent’s method in Python
                ele[nei] = optimize.brenth(f, h, min_h)
                
                # Appending the processed node to the queue
                queue.append(nei)
    return ele
    
def implicit_diffusion(C, dt):
    '''
    Function arguments:
    C  :  RHS (known) column vector for the diffusion and uplift term
    dt :                                        Time-step to be taken
    '''
    
    kd   = D*dt/(dx**2)
    d1m1 = np.repeat(-kd, len(mg.core_nodes))
    d1m1[ncols-3:d1m1.size:ncols-2] = 0
    # Constructing the LHS (5-diagonal sparse matrix) for fixed boundary condtion in the rectangular domain)
    Am = scipy.sparse.diags([-kd, d1m1, 1 + 4*kd, d1m1, -kd], [-ncols+2, -1, 0, 1, ncols-2], shape=(len(mg.core_nodes), len(mg.core_nodes)))
   
    return scipy.sparse.linalg.lgmres(Am, C, atol=0.000000000001)[0]

# Initializing the coefficients, domain-size and model parameters
K_sp  = 0.000025; D = 0.005; U = 0.001; m_sp  = 1; n_sp  = 1; dx = 1.0; Nr = 101; Nc = 101

# Calculating "Xi" value 
Xi   = (K_sp*((Nc-1)*dx)**(m_sp + n_sp))/ (D**n_sp * U**(1-n_sp))
print("Xi value for this simulation run is ",Xi)

# Initializing the grid
mg = RasterModelGrid((Nr,Nc), dx)

# Initializing the elevation values
_  = mg.add_zeros('node', 'topographic__elevation')
z  = np.zeros((Nr, Nc))
mg.at_node['topographic__elevation'] = z.reshape(Nr*Nc)

# Imposing fixed elevation boundary conditions
for edge in (mg.nodes_at_left_edge, mg.nodes_at_right_edge):
    mg.status_at_node[edge] = FIXED_VALUE_BOUNDARY
for edge in (mg.nodes_at_bottom_edge,mg.nodes_at_top_edge):
    mg.status_at_node[edge] = FIXED_VALUE_BOUNDARY

# Selecting D_infinity as the flow direction method
fc = FlowAccumulator(mg, flow_director = 'FlowDirectorDINF')
fc.run_one_step()

# Minimum time-steps for simulation to run
min_try = 500

# Other important variables for the simulation
i = -1; t = 0; dt = 1.; diff_list = []; steady_state = False

while steady_state is False:
    
    # Selecting 'dt' for high accuracy (user dependent)
    i    += 1
    alpha = min((i + 1) * 0.01, 1.)
    dt    = alpha * min(dx / (K_sp * np.max((mg.at_node['drainage_area'] / dx) ** m_sp)) , 20. * dx ** 2 / (2 * D))
    dt    = min(dt, 500.)

    # Storing the previous time-step's elevation values
    ele_1 = np.copy(mg.at_node['topographic__elevation'])

    erode_done = False
    while erode_done is False:
        try:
            # Updating elevation using the erosion term (calling DinfEroder function)
            fc = FlowAccumulator(mg, flow_director = 'FlowDirectorDINF')
            fc.run_one_step()
            mg.at_node['topographic__elevation'] = DinfEroder(mg, dt, K_sp, m_sp, n_sp)
            
            # Updating elevation using the uplift and diffusion terms implicitly
            mg.at_node['topographic__elevation'][mg.core_nodes] = implicit_diffusion(mg.at_node['topographic__elevation'][mg.core_nodes] + U * dt, dt)
            erode_done = True

        except ValueError:
            print('Adjusting dt')
            mg.at_node['topographic__elevation'] = ele_1
            dt = dt / 2.

    # Copying the updated elevation values after the time-step
    ele_2 = np.copy(mg.at_node['topographic__elevation'])
    t    += dt

    # Calculating maximum change in elevation value at any node in the domain
    ele_diff          =         np.abs(ele_1 - ele_2).max()
    
    # Computing the change in the mean elevation value over the time-step
    ele_diff_mean     = np.abs(ele_1.mean() - ele_2.mean())
    
    diff_list.append([t, ele_diff, ele_diff_mean])
    
    # Printing the time and other important information to see the simulation's progress
    print(i, t, round(ele_diff/ (U * 0.001 * dt), 1), round(ele_diff_mean/ (U * 0.0001 * dt), 1))
    
    # Checking if the steady-state is reached
    if  i >= min_try and ele_diff < U * 0.001 * dt and ele_diff_mean < 0.0001 * U * dt:
        print("Steady-state is reached, relevant arrays will be stored!")
        steady_state = True

# Saving the steady-state elevation field and accumulation as numpy arrays
a_1 = './square_' + str(round(Xi)) + '_steady_ele_at_' + str(int(t)) + '.npy'
a_2 = './square_' + str(round(Xi)) + '_steady_acc_at_' + str(int(t)) + '.npy'
np.save(a_1, mg.at_node['topographic__elevation'].reshape(Nr,Nc))
np.save(a_2,         mg.at_node['drainage_area'].reshape(Nr, Nc))
