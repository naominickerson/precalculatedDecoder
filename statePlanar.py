import numpy as np
#import matplotlib.pyplot as plt


class PlanarLattice:
    """ State is reponsible for holding the configuration
    of a planar lattice - both errors and stabiliser states.
    
    The lattice has side length L, which must be odd. 

    The qubits lie in positions (i,j), where i + j is even.
    Their state is encoded as:
        - 0: no error
        - 1: a Z flip
        - 2: an X flip
        - 3: a Y flip (X and Z)

    The X plaquettes lie at (i, j), where i is even and j is odd.
    They fire if surrounded by an odd number of {Z,Y} flips

    The Z plaquettes lie at (i, j), where i is odd and j is even. 
    The fire if surrounded by and odd numper of {X, Y} flips

    Measure_vert_z detects any logical vertical z errors, by measuring for zs along a vertical z row.
    Apply_vert_z applies a vertical z error, by flipping along a horizontal x row
    
    X 1 X 1 X 1     <-- this is a vertical Z error
    . Z . Z . Z         it can't be seen by the Xs
    X . X . X .
    . Z . Z . Z
    X . X . X .
    . Z . Z . Z
    
    """
    def __init__(self, L):
        self.L = L
        self.array = np.zeros((L,L), dtype='uint8')
        self._matching = None
        self.error_types = ['None', 'X']
        self.x_i_indices = range(0, L, 2)
        self.x_j_indices = range(1, L, 2)
        self.z_i_indices = range(1, L, 2)
        self.z_j_indices = range(0, L, 2)
        self.x_stabiliser_indices = [(i,j) for i in self.x_i_indices for j in self.x_j_indices]
        self.z_stabiliser_indices = [(i,j) for i in self.z_i_indices for j in self.z_j_indices]
        self.qubit_indices = [(i,j) for j in range(0, self.L) for i in range(j%2, self.L, 2)]
        self._n_errors = None

    def add_errors(self,p):
        self._n_errors = None # reset error_count
        n_qubits = len(self.qubit_indices)
        errors = np.random.rand(n_qubits) < p
        self.array[zip(*self.qubit_indices)] = errors

    def apply_y_errors(self,locations):
        self._n_errors = None # reset error_count
        n_qubits = len(self.qubit_indices)
        errors = np.random.rand(n_qubits) < p
        self.array[zip(*self.qubit_indices)] = errors

    
    @property
    def matching(self):
        if self._matching is None:
            self.generatematching()
        return self._matching

   

    def neighbours(self, i, j):
        if i==0:
            return [(i, j+1), (i+1,j), (i,j-1)] 
        elif j==0:
            return [(i-1, j), (i, j+1), (i+1,j)] 
        elif i==self.L-1:
            return [(i-1, j), (i, j+1), (i,j-1)] 
        elif j==self.L-1:
            return [(i-1, j), (i+1,j), (i,j-1)] 
        else:
            return [(i-1, j), (i, j+1), (i+1,j), (i,j-1)] 

   
    def site_type(self, i, j):
        """ Returns:
                0 - for a qubit
                1 - for an X plaquette (as they detect z errors)
                2 - for a Z plaquette
        """
        i_mod2, j_mod2 = i%2, j%2
        if i_mod2==0 and j_mod2==1:
            return 1
        elif i_mod2==1 and j_mod2==0:
            return 2
        else:
            return 0





    # Querying
    # ========
    
    def n_errors(self):
        # cast as an int, otherwise it returns a uint8, which
        # leads to all types of problems if you try to do 
        # arithmetic (e.g. 113-115 = 383498534198239348)
        if self._n_errors is None:
            self._n_errors = int(np.sum(self.array[zip(*self.qubit_indices)]))
        return self._n_errors

    def logical_error(self):
        return 'X' if self.has_logical_x_error() else 'None'
    """ 
    Measure_vert_z detects any logical vertical z errors, by measuring for zs along a vertical z row.
    Apply_vert_z applies a vertical z error, by flipping along a horizontal x row
    
    X 1 X 1 X 1     <-- this is a vertical Z error
    . Z . Z . Z         it can't be seen by the Xs
    X . X . X .
    . Z . Z . Z
    X . X . X .
    . Z . Z . Z
    
    """



    def measure_z(s,m=None):
        m=s.matching if m==None else m
        error_sum=0
        #Z logical error is horizontal (eg i fixed)
        # measure Zs along a column
        j=s.z_j_indices[0]
        for i in s.x_i_indices:
            n=s.array[i,j]^m[i,j]
            error_sum=error_sum ^ n
        return error_sum & 1 == 1

    def measure_x(s,m=None):
        m=s.matching if m==None else m
        error_sum = 0
        # X logical error is vertical (e.g. j fixed)
        # measure Xs along a row
        i=s.x_i_indices[0]
        for j in s.z_j_indices:
            n = s.array[i, j]^m[i,j]
            error_sum = error_sum ^ n
        return error_sum & 2 == 2            
        

    def add_x(s):
        s._n_errors = None # reset error_count
        # X-flip the qubits on a Z-column
        j = s.z_j_indices[0]
        for i in s.x_i_indices:
            s.array[i, j] = s.array[i,j] ^ 2
        return s
    
    def add_z(s):
        s._n_errors = None # reset error_count
        # Z-flip the qubits on an X-row
        i = s.x_i_indices[0]
        for j in s.z_j_indices:
            s.array[i, j] = s.array[i,j] ^ 1
        return s



    
    # Actions
    # =======

    def measure_stabiliser(s, x, y):

        #print (x,y),s.site_type(x,y) ,[s.array[c]  for c in s.neighbours(x,y)],reduce(np.bitwise_xor,[s.array[c]  for c in s.neighbours(x,y)])
        return bool(s.site_type(x,y) & reduce(np.bitwise_xor,[s.array[c]  for c in s.neighbours(x,y)]) ) 


    def get_syndrome(self):

        stabilisers=self.x_stabiliser_indices+self.z_stabiliser_indices
        coords=[]    
        for i,j in stabilisers:
            
            if self.measure_stabiliser(i,j):
                self.array[i,j]=1
                #syndrome+=[1]
                coords.append((i,j))
            else:
                self.array[i,j]=0
                #syndrome+=[0]
            
        return coords
                
                          

            

    def generate_syndrome(s):
        coords = []
        for i, j in s.x_stabiliser_indices + s.z_stabiliser_indices:
            if s.measure_stabiliser(i,j):
                s.array[i, j] = 1
                coords.append((i,j))
            else:
                s.array[i, j] = 0
        return coords


    def matching_from_syndrome(self,coords):
        self._n_errors = None
        m= np.zeros(np.shape(self.array), dtype='uint8')
        for I, J in coords:
            e=self.site_type(I,J)
            if e==2:
                for i in range(0, I, 2): #know first qubit is at 1
                    m[i, J] = m[i, J] ^ e
                    

            if e==1: 
                for j in range(0, J, 2):
                    m[I, j] = m[I, j] ^ e
        return m

    
    def generate_matching(self):
        self._n_errors = None
        coords = self.generate_syndrome()
        m = self._matching = np.zeros(np.shape(self.array), dtype='bool')
        # find coords of x anyons
        # for each one Z flip the qubits required to 
        # connect it to (0,0)
        for I, J in coords:
            # Z flip first row up to I
           
            for i in range((I%2)+1, I+1, 2): #know first qubit is at 1
                m[i, (I+1)%2] = m[i, (I+1)%2] ^ 1
            # Z flip Ith column up to J
            for j in range((I+1)%2+1, J+1, 2):
                m[I, j] = m[I, j] ^ 1
        return m
    

    def apply_stabiliser(s, x, y):
        s._n_errors = None # reset error_count
        e=1 if s.site_type(x,y) is 2 else 2
        
        for i,j in s.neighbours(x,y):
            s.array[i,j] = s.array[i,j] ^ e
        return s

    def count_errors(s):
        # x,y and z errors each count as 1
        return sum([1 for q in s.array[zip(*s.qubit_indices)] if q!=0])
    

    def copy(self):
        s = self.__class__(self.L)
        self.copy_onto(s)
        return s

    def copy_onto(self, other_state):
        other_state.array = self.array.copy()
        # note that matching is NOT copied - it's a reference
        # this is fine - each time we call generate_matching()
        #                a new one is created
        other_state.matching = self.matching
        return other_state
### done to here

    

    # Displaying
    # ==========
    def dump_s(self):
        return __class__.a_to_s(self.array)

    def to_n(self):
        # WARNING - only represents x errors
        ans = 0
        for i,j in self.qubit_indices:
            ans = ans*2 + self.array[i,j]
        return ans

    def show(self):
        fig = plt.gcf()# or plt.figure()
        fig.clf()
        ax = fig.add_subplot(111)
        # Decide how to show array:
        #  - 0 is a non-firing X stabiliser
        #  - 1 is a non-firing Z stabiliser
        #  - 2 is a non-error qubit
        #  - 6 is an error qubit
        #  - 4 is a firing stabiliser
        show_array = np.zeros((self.L, self.L), dtype='int')
        a = self.array
        for (i,j) in self.x_stabiliser_indices:
            show_array[i,j] = 4 if a[i,j]==1 else 0
        for (i,j) in self.z_stabiliser_indices:
            show_array[i,j] = 4 if a[i,j]==1 else 1
        for (i,j) in self.qubit_indices:
            show_array[i,j] = 6 if a[i,j]==1 else 2
        cax = ax.imshow(show_array, interpolation='nearest', vmax=6)
        cbar = plt.colorbar(cax, ticks=[0, 1, 2, 4, 6])
        cbar.ax.set_yticklabels(['X stab', 'Z stab','qubit', 'firing X stab', 'error qubit'])
        #for x, y in np.argwhere(self.matching > 0):
         #   circ = plt.Circle((y, x), radius = 0.3)
            # it seems that matplotlib transposes the coords
            # of an imshow such that the coords where 
            # a[i,j] are shown are actually (j, i)
            # this makes sense, as for arrays j corresponds
            # to horizontal movement, whereas the convention
            # for graph coords is (x=horiz, y=vert)
          #  ax.add_patch(circ)
        plt.show() # in case not already drawn
        plt.draw() # refresh if drawn
        return show_array

    @staticmethod
    def a_to_s(array):
        def translate(elt, i, j):
            if (i+j)%2 == 0: # X or Z stabiliser
                return "." if elt==0 else '1'
            else: # qubit
                if elt == 0:
                    return "."
                elif elt == 1:
                    return "Z"
                elif elt == 2:
                    return "X"
                else:
                    return "Y"
        return '\n'.join([" ".join([translate(elt, i, j) for i, elt in enumerate(row)]) for j, row in enumerate(array)])

    @staticmethod
    def s_to_a(string):
        def t(elt):
            if elt == ".":
                return 0
            elif elt == '1':
                return 1
            elif elt == 'Z':
                return 1
            elif elt == 'X':
                return 2
            else : #elt == 'Y'
                return 3
        return np.array([[t(elt) for elt in s.split(" ")] for s in string.strip().split('\n')])



### done from here
class UniformPlanarState(PlanarLattice):
    def __init__(self, L, p):
        PlanarLattice.__init__(self, L)
        self.p = p

    def generate_errors(self):
        self._n_errors = None # reset error_count
        n_qubits = len(self.qubit_indices)
        errors = np.random.rand(n_qubits) < self.p
        self.array[zip(*self.qubit_indices)] = errors



    def likelihood(self):
        n = self.n_errors()
        #N = len(self.qubit_indices)
        p = self.p/(1-self.p)
        return p**n #* (1-p)**(N-n)


    def relative_prob(self, s2):
        # returns p(s2)/p(self)
        n = self.n_errors()
        n2 = s2.n_errors()
        x = self.p/(1-self.p)
        diff = n2 - n
        #print(n2, n, self.p,x, diff, x**diff )

        return x**diff



    def copy(self):
        s = self.__class__(self.L, self.p)
        self.copy_onto(s)
        return s

    def generate_next(s):
        s._n_errors = None # reset error_count
        # pick a random z site
        n = len(s.z_stabiliser_indices)
        r = np.random.random_integers(0, n-1)
        x, y = s.z_stabiliser_indices[r]
        # apply the stabilizer at that point
        s.apply_stabiliser(x,y)


