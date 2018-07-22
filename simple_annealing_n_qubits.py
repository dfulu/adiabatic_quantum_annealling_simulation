# import dependencies
import numpy as np
import matplotlib.pyplot as plt
from progressbar import printProgress
import os.path
import scipy.linalg # this line is needed as well as the one below so that this imports in Jupyter correctly
import scipy as sp
from tabulate import tabulate
import time
import sys
print 'annealing loaded - '+'%s/%s/%s - '%time.localtime()[:3][::-1]+':'.join(('0'+str(x))[-2:] for x in list(time.localtime()[3:6]))
########################################### SIGMA MATRICES ###########################################################
# define the pauli matrices we are using
sigma_x = np.matrix('0. 1.; 1. 0.')
sigma_z = np.matrix('1. 0.; 0. -1.')

def sigma_zi(i,N):
    '''This returns the equivalent sigma_z for the ith qubit in the new Hilbert space given that there are N qubits'''
    return np.kron(np.eye(2**i), np.kron(sigma_z,np.eye(2**(N-(i+1)))))
    
def sigma_xi(i,N):
    '''This returns the equivalent sigma_x for the ith qubit in the new Hilbert space given that there are N qubits'''
    return np.kron(np.eye(2**i), np.kron(sigma_x,np.eye(2**(N-(i+1)))))
    
    
########################################### SAVE AND LOAD ############################################################
def save_project_data(filename,data):
    path = r"C:\Users\User\Documents\FIZZIX\4th Year\Project\Data"
    completeName = os.path.join(path,filename)
    np.savetxt(completeName,data)

def load_project_data(filename):
    path = r"C:\Users\User\Documents\FIZZIX\4th Year\Project\Data"
    completeName = os.path.join(path,filename)
    return np.loadtxt(completeName)



####################################### DIAGONALISATION FUNCTIONS ####################################################
# the Anneal class is dependent on the next two functions
def Diagonaliser(M):
    '''This function takes a hermitian matrix and return the unitary transform matrices to change to and from the eigenbasis
    as well as the diagonalised form of this matrix'''
    eigvals, eigvecs = sp.linalg.eigh(M)
    P = np.matrix(eigvecs)
    P_in = P.H
    return eigvals,P,P_in
    
def make_diag(d):
    '''this is a function to make a diagonal matrix given a series of values to go on the diagonal'''
    return np.matrix(np.diag(d))
    


####################################### ANNEALING SCHEDULE FUNCTIONS #################################################
# use the following function as the default anealing schedule
def linear_schedule(points):
    '''A(s) is indexed with return[0] B(s) is indexed with return[1]'''
    return np.asarray([np.ones(points) - np.linspace(0,1,points),np.linspace(0,1,points)])
    
def exponential_schedule(points):
    '''A(s) is indexed with return[0] B(s) is indexed with return[1]'''
    a = 5.
    b = 5.
    return np.asarray([1/(1-np.exp(-a))*(np.exp(-a*np.linspace(0,1,points))-np.exp(-a)),1/(np.exp(b)-1)*(np.exp(b*np.linspace(0,1,points))-1)])

def exponential_schedule_ab(a,b):
    '''This functional returns the previous exponential function but allows you to set the a and b params'''
    def exponential_schedule(points):
        '''A(s) is indexed with return[0] B(s) is indexed with return[1]'''
        return np.asarray([1/(1-np.exp(-a))*(np.exp(-a*np.linspace(0,1,points))-np.exp(-a)),1/(np.exp(b)-1)*(np.exp(b*np.linspace(0,1,points))-1)])
    return exponential_schedule

def interpolate_exper(exp_ss,exp_AB,points):
    '''This function takes the experimental AB you feed in which are discrete and interpolates them onto a new set of points'''
    new_ss = np.linspace(0,1,points)
    new_AB = np.zeros((2,points))
    new_AB[0] = sp.interp(new_ss, exp_ss, exp_AB[0])
    new_AB[1] = sp.interp(new_ss, exp_ss, exp_AB[1])
    return new_ss, new_AB

def vesuvius_schedule(points):
    '''This function returns the annealing schedule which is an interpolation of that used in the ~512 qubit D-Wave 2 (Vesuvius) 
    device (specifically the one previously at the University of Southern California, although different devices from the 
    same generation have nearly identical schedules)'''
    ves_SAB = load_project_data('ISI_A&B.txt') # this is the saved text file of the schedule used in the vesuvius
    ves_ss = ves_SAB[:,0]
    ves_AB = np.transpose(ves_SAB[:,1:])
    interp_ss, interp_AB = interpolate_exper(ves_ss,ves_AB,points)
    #interp_AB[0] =  interp_AB[0]/ interp_AB[0][0]
    #interp_AB[1] =  interp_AB[1]/ interp_AB[1][-1]
    return interp_AB        

def vesuvius_schedule_normed(points):
    '''This function returns the annealing schedule which is an interpolation of that used in the ~512 qubit D-Wave 2 (Vesuvius) 
    device (specifically the one previously at the University of Southern California, although different devices from the 
    same generation have nearly identical schedules)'''
    ves_SAB = load_project_data('ISI_A&B.txt') # this is the saved text file of the schedule used in the vesuvius
    ves_ss = ves_SAB[:,0]
    ves_AB = np.transpose(ves_SAB[:,1:])
    interp_ss, interp_AB = interpolate_exper(ves_ss,ves_AB,points)
    interp_AB[0] =  interp_AB[0]/ interp_AB[0][0]
    interp_AB[1] =  interp_AB[1]/ interp_AB[1][-1]
    return interp_AB     
    
#################################### BACK AND FORTH TO BIT SOLUTIONS #################################################
def bit_values(qubits,v):
    '''this function takes the vector in the hamiltonian space and turns it into an expectation bit value'''
    bits = np.zeros(qubits)
    for i in range(qubits):
        H = sigma_zi(i,qubits)
        bits[i] = 0.5*(1-v.getH()*H*v)
    return bits

def make_eigen_vector(bits):
    '''this function takes a string of bit assignments and turns them into a vector in the hamiltonian space'''
    neg = np.matrix('0. 1.')
    pos = np.matrix('1. 0.')
    vecs = [pos,neg]
    vec = vecs[bits[0]]
    for i in range(1,len(bits)):
        vec = np.kron(vec,vecs[bits[i]])
    return vec

def bit_table(last_state,problem_x0s, qubits):
    '''given a state this fucntion prints out the table of the probability of measuring this state in it's eigenfunctions'''
    headers = ['Bit Solution','Final State Probability','Problem X0 Probability']
    bits = []
    last_state_probs = []
    problem_x0_probs = []
    table = []
    for i in range(2**qubits):
        bits.append( [(i/2**j)%2 for j in range(qubits)[::-1]])
        last_state_probs.append((abs(last_state.getH()*np.transpose(make_eigen_vector(bits[-1])))**2).item(0))
        problem_x0_probs.append(sum([abs(np.dot(make_eigen_vector(bits[-1]),prob_x0_pos)**2).item(0) for prob_x0_pos in problem_x0s]))
        if problem_x0_probs[-1]==0: table.append([bits[-1],last_state_probs[-1],''])
        else: table.append([bits[-1],last_state_probs[-1],problem_x0_probs[-1]])
    
    table.insert(0,['','Result Highlights',''])
    table.insert(1,['','Full Results',''])
    nonzeroindex = np.where(np.asarray(problem_x0_probs) !=0)[0]
    for nzi in nonzeroindex:
        table.insert(1,[bits[nzi],last_state_probs[nzi],problem_x0_probs[nzi]])
    table.append(['Total_Prob - 1: ',float(np.sum(last_state_probs)-1),''])
    print tabulate(table, headers, tablefmt="grid")
    return tabulate(table, headers, tablefmt="grid")
    

#################################### KAPPA FUNCTION AND RELIABLES #################################################

def diff(x,y):
    'RETURNS THE BEST ESTIMATE FOR THE DIFFERENTIAL OF A SIMPLE 1D DATASET'
    dx0 = np.asarray([x[1]-x[0]]); dy0 = np.asarray([y[1]-y[0]])
    dxf = np.asarray([x[-1]-x[-2]]); dyf = np.asarray([y[-1]-y[-2]])
    dxs = x[2:]-x[:-2]; dys = y[2:]-y[:-2]
    return np.concatenate((dy0/dx0, dys/dxs, dyf/dxf))

def kappa(a):
        '''This function returns the max-kappa value throughtout the course of an anneal'''
        K = []
        #lets find out the differential of the schedule
        dA = diff(a.ss,a.AB[0])
        dB = diff(a.ss,a.AB[1])
        dAB = np.asarray([dA,dB])
        
        for i in range(a.points):
            # in this loop we calulate the eigen energy gap between (degenerate) ground state and the first excited level
            
            # neeed to be able to find the ground and frist excited states
            d,P,P_in = Diagonaliser(-a.AB[0][i]*a.Hb+a.AB[1][i]*a.Hp)
            E0 = min(d)
            x0 = P[:,d==E0] #x0 = P[:,d==E0] this line should only work if the ground state is not degenerate
            if len(np.argwhere(d==E0))>1: print '[%s] This is not valid, the ground state was degerate'%i
                
            # we want to consider the kappa value for all states wso that we can consider the maximum one
            pot_k = []
            for i0 in range(len(d)):
                if d[i0]==E0: # don't consider the overlap with itself or any degenerate levels
                    pass
                else:
                    dE = abs(E0-d[i0])
                    xn = P[:,i0] 
                    # here is the real line of calculation. This should be exactly the kappa equation
                    pot_k.append(abs(np.dot(np.transpose(x0),(-dAB[0][i]*a.Hb+dAB[1][i]*a.Hp)*xn).item()/(dE**2)))
            k = max(pot_k)
            K.append(k)
        return np.asarray(K)
    
def schedule_from_const_quant(Q,points):
        '''This function returns a schedule that will keep a quantity constant given that it was generated 
        from the linear schedule'''
        global sprime
        sprime = np.cumsum(Q)
        sprime = sprime-sprime[0]
        sprime = sprime/sprime[-1]
         # note that to make general for non-linear feed in schedules replace the last entry in this function below with AB[1]
        B = np.interp(np.linspace(0,1,points), sprime, np.linspace(0,1,len(sprime)))  
        AB = np.asarray([-B+1,B])
        return AB

class Anneal():
    '''
    This object is to keep all the annealing stuff together in one neat unit for each anneal completed.
    It applies an anneal for any number of bits as long as you give it the right input. 
    
    Each function explains itself but generally: 
        -the run function does the anneal
        -show_results shows the resuls of the anneal
    
    You need to input the parameters in the form:
        params = [h,j]
        
    Furthermore J needs to be formatted correctly. If there are n bits and Jij is the interaction coefficient between the ith and jth bit,
    J can be in the following format: 
        J = [ [J12,J13,J14,...Jin],  [J23,J24,J25,...J2n],  .......,  [J(n-2)(n-1),J(n-2)n],  [J(n-1)n] ]
        Alternatively it can be in a square array where J[i,j] = Jij and all terms below the upper right triange are zeros (this includes the diagonal)
    '''
    
    def __init__(self,qubits,params, **kwargs):
        '''This function sets up the attributes of the class'''
        # management
        self.creation_date = time.gmtime()
        self.light = kwargs.get('light', False) # runs a version of this where we are only interested in the final success probability nothing else at all
        if type(self.light)==str: # should be 'DEG=N'
            self.DEG = int(self.light[4:])
            self.light=True
        self.PITNUM = kwargs.get('PROGRESS_IT_NUMBER', None)
        self.show_bar = kwargs.get('show_bar', True)
        
        # schedule functions
        self.diff_scheds = kwargs.get('diff_scheds', False) # I am assuming the form of sched_funcs will be [ [A1(s), A2(s),...,AQUBITS(S)] , [B(S)] ]
        if self.diff_scheds: schedules = kwargs.get('sched_funcs') 
        else: schedule = kwargs.get('sched_func', linear_schedule)        
        
        # parameters
        self.qubits = qubits
        self.points = kwargs.get('points', 10000)
        self.T = float(kwargs.get('T', 100))
        self.h = params[0]
        self.J = params[1]
        if type(params[1])==list: self.J = np.asarray([[0]*(i+1)+params[1][i] for i in range(qubits-1)]+[[0]*qubits])
        self.check_params() # checks validity
        self.hamiltonian_parts() # sets up Hp and Hb
        
        # set schedules for the anneal
        self.ss = np.linspace(0,1,self.points)
        if not self.diff_scheds:
            if schedule == 'inverse_energy':
                self.AB = self.inverse_energy_schedule(self.points)
            elif schedule == 'inverse_square_energy':
                self.AB = self.inverse_square_energy_schedule(self.points)
            elif schedule == 'constant_k':
                self.AB = self.constant_k_schedule(self.points)
            else:
                self.AB = schedule(self.points)
        else: 
            ABs = np.zeros((self.points,qubits+1))
            for i in range(qubits):
                ABs[:,i] = schedules[0][i](self.points)[0]
            ABs[:,-1] = schedules[1][0](self.points)[1]
            self.ABs = ABs
            
        # method Euler or RK4
        self.method = kwargs.get('method', 'Euler')
        
        
    
    def check_params(self):
      '''Checks to see of the parameters you have put in are consistant with the number of qubits you have put in'''
      cond1 = (self.J.shape==(self.qubits,self.qubits))
      cond2 = (len(self.h)==self.qubits)
      if not cond1: 
          print '\n\n\nThe length of the J parameters array does not match the number of qubits you have selected\n\n\n'
      if not cond2: 
          print '\n\n\nThe length of the h parameters array does not match the number of qubits you have selected\n\n\n'    
        
    def hamiltonian_parts(self):
        '''
        - Constructs the hamiltonian from the parameters input
        - If we are doing an anneal run where we have different annealing schedules for each qubits it will return a list 
        of hamiltonian parts for the base hamiltonian. Else it returns a single hamiltonian for this part
        '''
        #initiate the hamiltonians giving them the right dimensions
        Hb = []
        Hp = np.matrix(np.eye(2**self.qubits)*0.)
        h = self.h
        J = self.J
        # loop through adding pieces to construct the hamiltonians 
        for i in range(self.qubits):
            Hb.append(sigma_xi(i,self.qubits))
            Hp += h[i]*sigma_zi(i,self.qubits)
            # the following part adds the interaction between all qubits. There is 1 interaction coefficient for each pairing 
            for j in range(i+1,self.qubits): 
                Hp +=J[i][j]*sigma_zi(i,self.qubits)*sigma_zi(j,self.qubits)     
        if not self.diff_scheds: Hb = sum(Hb)
        #set attibutes to store data   
        self.Hb = Hb
        self.Hp = Hp
    
    def H(self,i):
        if not self.diff_scheds:
            Hi = -self.AB[0][i]*self.Hb + self.AB[1][i]*self.Hp
        else: 
            Hi = sum([-self.ABs[:,j][i]*self.Hb[j] for j in range(self.qubits)])+self.ABs[:,-1][i]*self.Hp
        return Hi
        
    
    def inverse_energy_schedule(self, points):
        '''This function returns the AB annealing schedule so that each time step is proportional to the inverse energy gap at that time'''
        sys.stdout.write('Calculating Annealing Schedule');sys.stdout.flush()
        AB = linear_schedule(int(points/10.)) # this is just used to get the eigenergies with different A(s) and B(s) values
        delta_eigenvals = [] # used to store the data
        
        for i in range(int(points/10.)):
            # in this loop we calulate the eigen energy gap between (degenerate) ground state and the first excited level
            d = sp.linalg.eigvalsh(-AB[0][i]*self.Hb+AB[1][i]*self.Hp)
            delta_eigenvals.append(min(d[d!=min(d)])-min(d))
        
        # here we do the work to set the schedule
        inverse_delta_eigenvals = np.asarray(delta_eigenvals).astype(float)**-1
        AB = schedule_from_const_quant(inverse_delta_eigenvals,points)
        return AB
    

    
        
    def constant_k_schedule(self, points):
        '''This function returns the AB annealing schedule so that the kappa parameter is kept constant.'''
        sys.stdout.write('Calculating Annealing Schedule \n');sys.stdout.flush()
        self.AB = linear_schedule(points) # this is just used to get the eigenergies with different A(s) and B(s) values
        K =  kappa(self)
        # here we do the work to set the schedule
        AB = schedule_from_const_quant(K,points)
        sys.stdout.write('');sys.stdout.flush()
        
        return AB
        
    
    def inverse_square_energy_schedule(self, points):
        '''This function returns the AB annealing schedule so that each time step is proportional to the inverse energy gap at that time'''
        sys.stdout.write('Calculating Annealing Schedule');sys.stdout.flush()
        AB = linear_schedule(int(points/10.)) # this is just used to get the eigenergies with different A(s) and B(s) values
        delta_eigenvals = [] # used to store the data
        
        for i in range(int(points/10.)):
            # in this loop we calulate the eigen energy gap between (degenerate) ground state and the first excited level
            d = sp.linalg.eigvalsh(-AB[0][i]*self.Hb+AB[1][i]*self.Hp)
            delta_eigenvals.append(min(d[d!=min(d)])-min(d))
        
        # here we do the work to set the schedule
        inverse_delta_square_eigenvals = np.asarray(delta_eigenvals).astype(float)**-2
        AB = schedule_from_const_quant(inverse_delta_square_eigenvals,points)
        return AB
    
    def euler_update(self,dt,i):
        # mathematical pieces to diagonalise
        d,P,P_in = Diagonaliser(self.H(i))
        # store energy eigenvalues
        self.eigenvals.append(d)
        # store the new state of the qubits
        self.states.append(P*make_diag(np.exp(-1.j*d*dt))*P_in*self.states[-1])
        # the following is not needed in light mode
        if not self.light:
            # find instantaneous ground states
            self.instant_x0s.append([P[:,d==min(d)][:,k] for k in range(P[:,d==min(d)].shape[1])]) #  this line can't come before the line above or things will change due to the apend. This is the current state in the diagonal basis
            # find and store the probability of measuring the qubits to be in the instantaneous ground state
            self.instant_x0s_prob.append(sum([abs(np.dot(self.states[-1].getH(),ins_x0_pos)**2).item(0) for ins_x0_pos in self.instant_x0s[-1] ])) # this line finds the probability of the state being in the instantaneous ground state
            # find and store the probability of measuring the qubits to be in the problem ground state
            self.problem_x0_prob.append(sum([abs(np.dot(self.states[-1].getH(),prob_x0_pos)**2).item(0) for prob_x0_pos in self.problem_x0s]))
            # find and store the energy eigenvalue gap
            self.delta_eigenvals.append(min(d[d!=min(d)])-min(d))
        
    def rk4_update(self,dt,i):
        
        d1,P1,P_in1 = Diagonaliser(self.H(i))
        K1 = (P1*make_diag(np.exp(-1.j*d1*dt))*P_in1)*self.states[-1]-self.states[-1]
        
        d2,P2,P_in2 = Diagonaliser(0.5*(self.H(i)+self.H(i+1)))
        K2 = (P2*make_diag(np.exp(-1.j*d2*dt))*P_in2)*(self.states[-1]+1/2.*K1)-(self.states[-1]+1/2.*K1)
        
        d3,P3,P_in3 = d2,P2,P_in2
        K3 = (P3*make_diag(np.exp(-1.j*d3*dt))*P_in3-np.matrix(np.eye(2**self.qubits)))*(self.states[-1]+1/2.*K2)
        
        d4,P4,P_in4 = Diagonaliser(self.H(i+1))
        K4 = (P4*make_diag(np.exp(-1.j*d4*dt))*P_in4-np.matrix(np.eye(2**self.qubits)))*(self.states[-1]+K3)
        if i==10:
            global KK
            KK={
                'd1':d1,'P1':P1, 'P_in1':P_in1,
                'd2':d2,'P2':P2, 'P_in2':P_in2,
                'd3':d3,'P3':P3, 'P_in3':P_in3,
                'd4':d4,'P4':P4, 'P_in4':P_in4,
                'K1':K1,'K2':K2,'K3':K3,'K4':K4,
                'state':self.states[-1]                
                }
            print K2.getH()*K2
        self.states.append(self.states[-1]+K2)#1./6.*(K1+2*K2+2*K3+K4))
        
        # store energy eigenvalues
        self.eigenvals.append(d1)
        # the following is not needed in light mode
        if not self.light:
            # find instantaneous ground states
            self.instant_x0s.append([P1[:,d1==min(d1)][:,k] for k in range(P1[:,d1==min(d1)].shape[1])]) #  this line can't come before the line above or things will change due to the apend. This is the current state in the diagonal basis
            # find and store the probability of measuring the qubits to be in the instantaneous ground state
            self.instant_x0s_prob.append(sum([abs(np.dot(self.states[-1].getH(),ins_x0_pos)**2).item(0) for ins_x0_pos in self.instant_x0s[-1] ])) # this line finds the probability of the state being in the instantaneous ground state
            # find and store the probability of measuring the qubits to be in the problem ground state
            self.problem_x0_prob.append(sum([abs(np.dot(self.states[-1].getH(),prob_x0_pos)**2).item(0) for prob_x0_pos in self.problem_x0s]))
            # find and store the energy eigenvalue gap
            self.delta_eigenvals.append(min(d1[d1!=min(d1)])-min(d1))
        return 
        
    def run(self):
        '''this function runs the annealing protocol and creates attributes to store the data within the class'''
        
        # find the ground state of the base hamiltonian
        base_eigen_vals0,base_eigen_vecs0 = np.linalg.eigh(self.H(0)) # the minus sign before the base hamiltonian is because we use H(s) = -A(s)Hb + B(s)Hp and H(0) = -Hb
        self.base_x0 = base_eigen_vecs0[:,base_eigen_vals0.argmin()]+0.j*base_eigen_vecs0[:,base_eigen_vals0.argmin()]
        
        # we need to distinguish between the ground state of the base hamiltonian and the base state of the hamiltonian at t=0 because these may be different
        start_eigen_vals0,start_eigen_vecs0 = np.linalg.eigh(self.H(0))
        self.start_x0 =start_eigen_vecs0[:,start_eigen_vals0.argmin()]+0.j*start_eigen_vecs0[:,start_eigen_vals0.argmin()]
        
        # find the ground state of the problem hamiltonian 
        d,P = np.linalg.eigh(self.Hp)
        self.problem_x0s = [P[:,d==min(d)][:,i] for i in range(P[:,d==min(d)].shape[1])]       
        
        # set the time steps
        dt = self.T/self.points
        
        # set up data stores for the anneal whilst setting the inital state to be th ground state on the base hamiltonian
        self.states = [self.base_x0] # store the vector wavestates
        self.instant_x0s = [] # store the instantaneous ground state

        self.eigenvals =[] # store the energy eigenvalues
        self.delta_eigenvals = []
        
        self.problem_x0_prob = [] # store the probabilities of being measured in the problem base state
        self.instant_x0s_prob = [] # store the probabilities of being measured in the instantaneous base state        
        
        # carry out the anneal whilst storing the state and the energy eigenvalues thoughout the process
        l = len(self.ss)-1
        
        if self.method == 'Euler':
            for i in range(len(self.ss)):
                # print loading bar but only update after each full percent has been completed
                if (((i+1)*100%(l+1)<=100) or i == l) and self.show_bar: printProgress (i, l, prefix = 'Annealing: ', suffix = 'est_time', decimals = 0, barLength = 50, PROGRESS_IT_NUMBER = self.PITNUM)
                self.euler_update(dt,i)
                if self.light and i == self.points-1: # calculate only the last probability of success
                    self.problem_x0_prob.append(sum([abs(np.dot(self.states[-1].getH(),prob_x0_pos)**2).item(0) for prob_x0_pos in self.problem_x0s])) 
                
        elif self.method == 'RK4':
            for i in range(len(self.ss)):
                # print loading bar but only update after each full percent has been completed
                if (((i+1)*100%(l+1)<=100) or i == l) and self.show_bar: printProgress (i, l, prefix = 'Annealing: ', suffix = 'est_time', decimals = 0, barLength = 50, PROGRESS_IT_NUMBER = self.PITNUM)
                
                # update and store the data with this function
                if i != self.points-1: # the Runge Kutta technique won't work on the last step since it needs to index i+1
                    self.rk4_update(dt,i)
                else:
                    self.euler_update(dt,i)
                
                if self.light and i == self.points-1: # calculate only the last probability of success
                    self.problem_x0_prob.append(sum([abs(np.dot(self.states[-1].getH(),prob_x0_pos)**2).item(0) for prob_x0_pos in self.problem_x0s]))
        # set attributes to store the data
        self.states = self.states[1:]
        # self.instant_x0s = instant_x0s I have decided I don't need these stored

        self.eigenvals = np.asarray(self.eigenvals)
        self.delta_eigenvals = np.asarray(self.delta_eigenvals)
        
        self.problem_x0_prob = np.asarray(self.problem_x0_prob)
        self.instant_x0s_prob = np.asarray(self.instant_x0s_prob)
        
        
        del self.PITNUM # Don't need this anymore. It was just to see when running
        del self.show_bar
        
    def show_results(self):
        '''This function outputs all of the important data, either printed or graphed'''
        
        # plot all of the energy eigenvalue traces
        print 'last set of eigenvalues were:',self.eigenvals[-1]
        plt.figure()
        plt.xlabel('s')
        plt.ylabel('$\lambda$')        
        for i in range(len(self.eigenvals[0,:])):
            plt.plot(self.ss,self.eigenvals[:,i])
        plt.show()
        
        # plot the probability of being found in the problem ground state/instanteaous gorund state
        plt.figure()
        ax1 = plt.subplot2grid((3,1), (0, 0), rowspan=2); ax2 = ax1.twinx()
        ax1.set_ylabel('$P(Problem\,ground\,state)$', color='r') 
        ax1.plot(self.ss, self.problem_x0_prob, 'r-')
        ax2.set_ylabel('$P(Instantaneous\,ground\,state)$', color='b')       
        ax2.plot(self.ss, self.instant_x0s_prob, 'b-')
        
        # plot the energy gap between the first 2 levels and the annealing schedule if it is uniform across bits
        ax3 = plt.subplot2grid((3,1), (2, 0), rowspan = 1,sharex = ax1)
        if not self.diff_scheds : 
            ax4 = ax3.twinx()
            ax4.set_ylabel('A(s), B(s)')
            ax4.plot(self.ss,self.AB[0], label = '$A(s)$',color = 'red')
            ax4.plot(self.ss,self.AB[1], label = '$B(s)$', color = 'blue')
        ax3.set_ylabel('$\Delta \lambda$',color = 'black')
        ax3.set_xlabel('s')     
        ax3.plot(self.ss,self.delta_eigenvals,color='black',label = '$\Delta \lambda$')    
        plt.legend(framealpha = 0, loc = 'best')
        plt.show()
        
        # if there has been different schedules we will plot these schedules seperately
        if self.diff_scheds:
            ax5 = plt.subplots()[1]
            ax5.set_ylabel('A(s), B(s)')
            for i in range(self.qubits):
                ax5.plot(self.ss,self.ABs[:,i], label = str(i+1),color = 'red', alpha = 0.5)
            ax5.plot(self.ss,self.ABs[:,-1], label = '$B(s)$', color = 'blue')
        plt.show()
        
        # print the bit table and the final probability of being in the ground state
        bit_table(self.states[-1],self.problem_x0s, self.qubits)
        print '\nProbability of being measured in lowest energy state: %s' %self.problem_x0_prob[-1]

if __name__ == "__main__": 
    #x = np.pi/6
    #qubits = 8
    #h = np.random.random(qubits)
    #J = np.random.random((qubits,qubits))
    #sfB = lambda points: exponential_schedule_ab(2,2)(points)
    #sfl = lambda points: np.cos(x)*exponential_schedule_ab(2,2)(points)
    #sfa = lambda points: np.sin(x)*exponential_schedule_ab(2,2)(points)
    
    # this following section defines the neutral hamiltonian 
    #define number of logical and ancilla qubits
    log_qub = 3
    anc_qub = log_qub
    qubits = log_qub*2
    
    # define the coupling parameters that comply with the conditions above
    JN = 0.1
    q0 = 50.
    Ja = 100.
    Jl = Ja
    hl = -Ja+q0
    
    
    # construct h
    hl = np.asarray([-Ja+q0]*log_qub) # h for the logical qubits
    ha = np.asarray([-Ja*(2*i-log_qub)+q0 for i in range(1,log_qub+1)]) + np.asarray([JN if (log_qub-i)%2 == 1 else -JN for i in range(1,log_qub+1)])# h for the ancilla qubits
    h = np.concatenate((hl,ha)) # the full h
    
    # constuct J
    J_log_log = Jl*np.triu(np.ones((log_qub,log_qub)), k=1) # the J component for the logical qubit interactions
    J_log_anc = Ja*np.ones((log_qub,log_qub)) # interactions between logical and ancilla bits
    Jt = np.concatenate((J_log_log,J_log_anc),axis=1) # make top half of J array
    Jb = np.zeros(Jt.shape) # make bottom half of J array
    J = np.concatenate((Jt,Jb),axis=0) # the full J

   # sched_funcs = [[sfl]*log_qub+[sfa]*anc_qub,[sfB]]
    a = Anneal(qubits,[h,J],T=10.,points = 1000, diff_scheds = False, show_bar = True, light = False) 
    a.run()
    a.show_results()
    #test1 = Anneal(qubits,[h,J],T=100,points = 1000, light = False,sched_func = sf1)
    #test1.run()
    #test2 = Anneal(qubits,[h,J],T=100,points = 250, light = False,sched_func = sf2)
    #test2.run()
    #raw_input('type to continue')
    #test1.show_results()
    #test2.show_results()

    #[plt.close(i) for i in [1,3]]
