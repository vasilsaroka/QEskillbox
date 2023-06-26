# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a Differential Evolution script file.
"""
import os
import subprocess
import re
from os.path import dirname, abspath

# QE run with an arbitrary number of the parameters
def QEGenrun(params,
          int_param_num = 1,
          place_holders = ["@ecutwfc", "@celldm", "@kpt"],
          pwx_path = "pw.x",
          input_file = "template-scf.in", # must contain two compalsory place holders "@outdir" and "@pseudo_dir" that need not to be specified in the place_holders kwarg.
          output_file = "template-scf.out",
          working_dir = ".",  # working directory is a more descriptive name for this option
          outdir = "",   
          pseudo_dir = ""):  
    """
    This function runs Quantum Espresso with a given template file as a input. 
    The given parameter values are substituted into a
    template instead of the placeholders starting from '@', e.g. @ecutwfc, 
    @celldm and etc. Then the function returns the used parameter values together
    with a total energy of the system.
    
    This funciton is to be used as the cost function in the differential evolution
    optimization of the parameters of interest.
    
    """
    plen = len(params)
    
    if plen != len(place_holders):
        raise Exception("The number of optimization parameters must be equal to \
                        the number of corresponding place holders. Check place_holders kwarg.")
    if int_param_num > plen :
        raise Exception("The number of integer parameters cannot exceed the dimensionality \
                        of parameter space. Check int_param_num kwarg and the length of the first arg.")
    
    # set the path to QE installation
    if pwx_path == "":
        cmd = subprocess.run(['which', 'pw.x'], 
                                 stdout=subprocess.PIPE, 
                                 universal_newlines=True)
        path = cmd.stdout.strip()
        
        if path == "":
            raise Exception("Failed to detect path to pw.x file. Install Quantum Espresso or provide path manually using pwx_path kwarg.")
        else: 
            pwx_path = path
    
    if (working_dir == "" or working_dir == "." or working_dir == "./"):
        working_dir = dirname(abspath(__file__))
                
    if pseudo_dir == "":
        pseudo_dir = '"' + working_dir + '"'
        
    if outdir == "":
        outdir = '"' + os.path.join(working_dir,"tmp") + '"'
    
    st = open(os.path.join(working_dir,input_file),'r')    
    IN = st.read() # this list still contains special characters '\n'
    st.close()
    
    plc_holders = place_holders + ["@outdir", "@pseudo_dir"]
    paramstr = []
    for i in range(plen):
        if i < plen-int_param_num:
            paramstr.append(str(params[i]))
        else:
            paramstr.append(str(int(params[i])))
    
    paramstr.append(outdir)
    paramstr.append(pseudo_dir)
    
    for i in range(len(paramstr)):
        IN = re.sub(plc_holders[i],paramstr[i],IN)
    
    # print input file for the testing purposes
    #print(IN)
    
    pwx = subprocess.Popen([pwx_path],
                            stdin =subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            universal_newlines=True,
                            bufsize=0)
     
    # Send qe commands to stdin
    pwx.stdin.write(IN)
    pwx.stdin.close()
    
    
    # in what follow two option are possible either to save output or avoid this
    e_tot = "";
    if output_file != "":
        # the output file is produced and saved in the outdir, for testing purposes only
        # avoid files overwriting up to 100 runs
        outfname = output_file
        for i in range(100):
            if os.path.isfile(os.path.join(working_dir,outfname)):
                outfname = output_file + '_' + str(i)
            else:
                break
        
        
        f = open(os.path.join(working_dir,outfname),'a')    
        # this list still contains special characters '\n'
        for line in pwx.stdout:
            f.write(line) 
            if line.find("!") != -1:
                e_tot = float(re.split(r"\s+",line.strip())[-2])
            
        f.close()
    else:
        # the output file is not produced, should be faster
        # this list still contains special characters '\n'
        for line in pwx.stdout:
            if line.find("!") != -1:
                e_tot = float(re.split(r"\s+",line.strip())[-2])
    
    if e_tot == "":
        raise Exception("Total energy was not found in the Quantum Espresso output, probably due to the lack of self-consistent run convergence or a crash. Check the given script parent directory for the Quantum Espresso CRASH file.")
        
    return params, e_tot

# 2023-06-16 
# DE optimization function for an arbitrary number of QE parameters; 
# one, few or all parametes can be chosen to be integer, counting from the end of
# the parameter space vector

import numpy as np
import time
from datetime import datetime
import matplotlib.pyplot as plt
 
def QEDE(param_space,
          popsize = 4,       # population size cannot be less than 4
          mutation = 0.5,    # mutation parameter chosen from [0,2] interval
          crossover = 0.7,   # crossover probability chosen from [0,1] interval
          max_iter = 30,    # maximum number of iterations to break the loop in case convergence is too slow 
          DE_tol = 1e-6,     # convergence threshold
          int_param_num = 1, # number of the integer parameters
          place_holders = ["@ecutwfc", "@celldm", "@kpt"], # place holders used in the QE input file template instead of the parametes of interest
          pwx_path = "",                  # path to QE pw.x is found automatically
          input_file = "template-scf.in", # QE input file template; it must contain two compalsory place holders "@outdir" and "@pseudo_dir" that need not to be specified in the place_holders kwarg.
          output_file = "",               # QE output file; can be used for tracking errors in QE runs
          working_dir = ".",  # working directory is a more descriptive name for this option
          outdir = "",        # QE output directory
          pseudo_dir = ""):   # directory with QE pseudo potentials 
    """
    This function runs Quantum Espresso within Differential Evolution (DE) Global 
    optimization algorithm. 
    The DE algorithm is searching for the Global minimum of the total energy of the
    system within the parameter space region specified as param_space.
    The parameter space vectors are substituted into a
    template input file instead of the placeholders that by convention start
    from '@', e.g. @ecutwfc, @celldm and etc. The Quantum Espresso is run on this
    input file and the total energy is returned. This cycle is repeated within
    DE algorithm until maximum number of iterations, max_iter, is achieved or the difference
    between total energies of all parameter space vectors in the given generation
    becomes lower than the given tolerance, DE_tol.
    The function returns an array of the parameter space vectors and their corresponding total
    energies for each generation as well as same quantities for last generation. The former 
    allows one to track the convergence of the DE, while the later to choose between the best found parameters.
    
    This funciton is to be used for finding optimal parameters of interest.
    
    """
    start_time = time.time()
    
    if (working_dir == "" or working_dir == "." or working_dir == "./"):
        working_dir = dirname(abspath(__file__))
    
    # Differential Evolution optimization
    D = len(param_space) # number optimized parameters, i.e. the dimensionality of parameter space

    # random sampling of the parameters within the zeroth generation
    X = np.asarray([(a[1] - a[0]) * np.random.random(popsize) + a[0] for a in param_space])
    # test initial parameter space vectors
    print("Initial parameter space vectors (in columns):")
    print(X)
    
    # initializing costs for every parameter space vector in the 0-th generation
    costs = np.zeros(popsize)
    for i in range(popsize):
        params = X[:,i]
        costs[i] = QEGenrun(params,
                            int_param_num,
                            place_holders,
                            pwx_path,
                            input_file,
                            output_file,
                            working_dir,
                            outdir,
                            pseudo_dir)[1]
    
    print("Total energies for the initial generation of the parameter space vectors:")
    print(costs)
    best_idx = np.argmin(costs)
    best_cost = costs[best_idx]
    best_vector = X[:,best_idx]
    
    # iterating generation untill the convergence is achieved
    # initializing a trial vector and a range of indexes numbering vectors in generation
    trial_vector = np.zeros(D)
    gen_range = np.arange(popsize)

    # introducing variable to collect the parameter space vectors and total energies thoughout all the iterations
    data = np.zeros((max_iter*popsize,D+1))
    
    # convergence loop
    # set est and k parameters to enter to loop
    est = DE_tol
    k = 0

    while (DE_tol <= est) and (k < max_iter):
         kp = k*popsize; # alias for k*popsize within the while loop
         # iterating through of the generation of parameter vectors choosing each vector to be once a target vector
         for i in range(popsize):
             # Choose r1, r2, r3 randomly so that they are different from each other and the running index i
             r1, r2, r3 = np.random.permutation(np.delete(gen_range,i))[0:3]

             # trial vector is the mixture of the mutant_vector and the target_vector
             rnbr = np.random.randint(D) # for the given vector, the random integer between 0 and D - 1
             for j in range(D):
                 if (np.random.random() > crossover) and (j != rnbr):
                     # Crossover with a target vector
                     trial_vector[j] = X[j,i]
                 else:
                     # Mutation between three other vectors that are chosen randomly
                     mutant = X[j, r1] + mutation * (X[j, r2] - X[j, r3])
                     # corrections of the parameter value, it must stay within parameter limits
                     if mutant < param_space[j, 0]:
                         trial_vector[j] = param_space[j, 0]
                     elif mutant > param_space[j, 1]:
                         trial_vector[j] = param_space[j, 1]
                     else:
                         trial_vector[j] = mutant
         
    
             print("trial vector = " + str(trial_vector))
             
             # Selection
             # If the trial vector minimizes the cost function compared to the target vectors than it replaces
             # the target vector in the next generation
             data[kp+i,:D], data[kp+i,D] = QEGenrun(trial_vector,
                                                          int_param_num,
                                                          place_holders,
                                                          pwx_path,
                                                          input_file,
                                                          output_file,
                                                          working_dir,
                                                          outdir,
                                                          pseudo_dir)
             trial_cost = data[kp+i,D]

             if trial_cost < costs[i]:
                 X[:,i] = trial_vector
                 costs[i] = trial_cost
                 if trial_cost < best_cost:
                     best_cost = trial_cost
                     best_vector = trial_vector
             else:
                     data[kp+i,:D] = X[:,i]
                     data[kp+i,D] = costs[i]
         
         print(costs)
         
         k += 1
         print("iteration No. " + str(k))
         # estimating convergence
         est = 0
         for i in range(popsize-1):
             est += np.abs(costs[0] - costs[i+1])
         print("Total energies spread = " + str(est))
    # clean data array from zeros, outside the while loop
    data = data[:k*popsize]
    # print the time spent on this search
    print("--- %s seconds ---" % (time.time() - start_time))
    
    # writing last generation of parameter space vectors into an output file
    # this piece of code must be tested 19/06/2023
    date = datetime.now().strftime('%Y-%m-%d_%H_%M_%S')
    fname = "QEDE_" + date + ".out"
    # writing data to the text file
    f = open(os.path.join(working_dir,fname),'w')
    f.write('Created by Vasil Saroka 40.ovasil@gmail.com')
    # using string literals
    f.write('\nToday is ' + date + '\n')
    f.writelines('\nDifferential Evolution results:\n')
    f.write('\nThe last generation of parameter space vectors:')
    
    for i in range(popsize):
        f.write('\nparameter vector ' + str(i+1) + '\n')
        for j in range(D):
            f.write('%5.9f\n' % X[j,i])
        f.write('\nThe total energy for vector ' + str(i+1) + ':\n')
        f.write('%5.9f\n' % costs[i])
    
    f.write('\nThe best vector is\n')
    for j in best_vector:
        f.write('%5.9f\n' % j)
    f.write('\nThe total energy for the best vector:\n')
    f.write('%5.9f\n' % best_cost)
    if k == max_iter :
        f.write("\nThe maximum number of iterations has been achieved. Try to increase this number.\n")
    f.write('\nDE is done!\n')
    f.close()
    
    # track convergence of the parameter space vectors during the DE evolution
    # create figure object
    fig, ax = plt.subplots()
    ax.set_xlim(0,k)
    ax.set_ylim(min(data[:,3]),max(data[:,3]))
    ax.set_xlabel("DE iteration")
    ax.set_ylabel(r"$E_{tot}$, Ry")

    for i in range(popsize):
        v = data[i:-1:popsize,3]
        ax.plot(v, alpha=0.9)
        
    fname ="QEDE_" + date + ".png"
    fig.savefig(os.path.join(working_dir,fname),dpi=200,bbox_inches='tight',transparent=False)
    
    return data, k, X, costs, best_vector, best_cost
#%% test on bulk Si

param_limits = np.array([[30, 80], [5, 15], [4, 16]])

# QE run options
placeholders = ["@ecutwfc", "@celldm", "@kpt"]
QEpath = "pw.x"
IN = "Si.scf-template.in"
OUT = ""    # no output file shall be generated to save the space on disk
wk_dir = "." # working directory is the directory of this script
pp_dir = ""  # implies pseudo potential to be located in the working directory

data, k, vecs, costs, best_vector, best_cost = QEDE(param_limits,
                                                    place_holders = placeholders,
                                                    pwx_path = QEpath, 
                                                    input_file = IN,
                                                    output_file= OUT,
                                                    working_dir= wk_dir, 
                                                    pseudo_dir=pp_dir)
#%%
print(vecs)
print(best_vector)

#%%
print(best_vector)
print(best_cost)

#%%

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.patches import FancyArrowPatch

class Arrow3D(FancyArrowPatch):

    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)
        
    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs) 
    


def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)


setattr(Axes3D, 'arrow3D', _arrow3D)

#%% working .gif-animation for QEGenrun

from matplotlib import pyplot as plt
from celluloid import Camera
import numpy as np

gsz = int(len(data)/k)
path = dirname(abspath(__file__))

# create figure object
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlim(param_limits[0,0],param_limits[0,1])
ax.set_ylim(param_limits[1,0],param_limits[1,1])
ax.set_zlim(param_limits[2,0],param_limits[2,1])



ax.set_xlabel('$E_{cut}$')
ax.set_ylabel('$a_{lat}$')
ax.set_zlabel('$k_{pt}$')

camera = Camera(fig)

for fr in range(k):
    # dynamic title for celluloid animation
    ax.text(80, 7, 108, str(fr)+'-th generation of parameter vectors', transform=ax.transAxes)
    
    pv = data[fr*gsz:(fr+1)*gsz]
    best_idx = np.argmin(pv[:,3])

    for i in range(gsz):
        if i == best_idx:
            ax.arrow3D(0,0,0,
               pv[i,0],pv[i,1],pv[i,2],
               mutation_scale=20,
               ec ='green',
               fc='red')
        else:
            ax.arrow3D(0,0,0,
               pv[i,0],pv[i,1],pv[i,2],
               mutation_scale=20,
               arrowstyle="-|>",
               linestyle='dashed')
            
    
    
    plt.pause(0.1)
    camera.snap()

animation = camera.animate()
animation.save(os.path.join(path,'animation_Si.gif'), writer='PillowWriter', fps=2, dpi = 200)


fig.tight_layout()

#%% working .mp4-video for QEGenrun
import matplotlib.animation as animation
from IPython.display import HTML

gsz = int(len(data)/k)
path = dirname(abspath(__file__))

fig = plt.figure()

ax = fig.add_subplot(111,projection='3d')


def animate(ind):
   ax.clear()
   ax.set_xlim(param_limits[0,0],param_limits[0,1])
   ax.set_ylim(param_limits[1,0],param_limits[1,1])
   ax.set_zlim(param_limits[2,0],param_limits[2,1])
    
   pv = data[ind*gsz:(ind+1)*gsz]
   best_idx = np.argmin(pv[:,3])

   for i in range(gsz):
       if i == best_idx:
           ax.arrow3D(0,0,0,
              pv[i,0],pv[i,1],pv[i,2],
              mutation_scale=20,
              ec ='green',
              fc='red')
       else:
           ax.arrow3D(0,0,0,
              pv[i,0],pv[i,1],pv[i,2],
              mutation_scale=20,
              arrowstyle="-|>",
              linestyle='dashed')
           

   ax.set_title(str(ind)+'-th generation of parameter vectors')
   ax.set_xlabel('$E_{cut}$')
   ax.set_ylabel('$a_{lat}$')
   ax.set_zlabel('$k_{pt}$')
   fig.tight_layout()
   
anim = animation.FuncAnimation(fig, animate, frames=k, interval=400)

print(HTML(anim.to_html5_video()))

plt.show()

anim.save(os.path.join(path,"video_Si.mp4"))

