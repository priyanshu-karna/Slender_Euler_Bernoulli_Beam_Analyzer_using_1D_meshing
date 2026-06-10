import numpy as np 
import matplotlib.pyplot as plt

class Beam:
    def __init__(self, L, E, I, P, load_node, n_elem, bc_type):
        self.L = L
        self.E = E
        self.I = I
        self.P = P                  
        self.load_node = load_node  
        self.n_elem = n_elem
        self.l_elem = L / n_elem
        self.ndof = 2 * (n_elem + 1)
        self.bc_type = bc_type
        
    def stiffness_matrix(self):
        l = self.l_elem
        k = (self.E * self.I) / l**3 * np.array([
            [12,      6*l,    -12,     6*l],
            [6*l,     4*l**2, -6*l,    2*l**2],
            [-12,    -6*l,     12,    -6*l],
            [6*l,     2*l**2, -6*l,    4*l**2]
        ])
        return k
    
    def assemble(self):
        
        K = np.zeros((self.ndof, self.ndof))
        f = np.zeros(self.ndof)
        
        k_elem = self.stiffness_matrix()
        
        
        for i in range(self.n_elem):
            idx = slice(2*i, 2*i + 4)
            K[idx, idx] += k_elem
            
        
        target_dof = 2 * self.load_node
        f[target_dof] = -self.P
            
        return K, f
    
    def get_bc_dofs(self):
        if self.bc_type == "cantilever":
            fixed_dof = [0, 1] 
        elif self.bc_type == "simply_supported":
            fixed_dof = [0, self.ndof - 2] 
        elif self.bc_type == "fixed_fixed":
            fixed_dof = [0, 1, self.ndof - 2, self.ndof - 1] 
        else:
            raise ValueError(f"Unknown boundary condition: {self.bc_type}")
            
        free_dof = [i for i in range(self.ndof) if i not in fixed_dof]
        return {"fixed_dof": fixed_dof, "free_dof": free_dof}
    
    def solve(self):
        K, f = self.assemble()
        d = np.zeros(self.ndof)
        bc = self.get_bc_dofs()
        
        free = bc["free_dof"]
        
        K_ff = K[np.ix_(free, free)]
        f_ff = f[free]
        
        d_free = np.linalg.solve(K_ff, f_ff)
        d[free] = d_free
        return d
    def plotting(self, d):
        x = np.linspace(0, self.L, self.n_elem + 1)
        deflection = d[0::2] 

        plt.figure(figsize=(8, 4))
        plt.plot(x, deflection, '-ro', label='FEM Point Load Deflection')
        plt.axvline(x=x[self.load_node], color='blue', linestyle='--', label='Point Load Location')
        plt.title(f'Beam Deflection Profile under Point Load ({self.bc_type.title()})')
        plt.xlabel('Beam Length (m)')
        plt.ylabel('Deflection (m)')
        plt.grid(True)
        plt.legend()
        plt.show()


# --- Example Usage ---
L = 10   
E = 210e9
I = 1e-6
n_elem = 10

P = 5000            # 5 kN point load
load_node = 5       # Must be between 0 and n_elem 
bc_type = "fixed_fixed"  # Change as described in the class for different boundary conditions

beam1 = Beam(L, E, I, P, load_node, n_elem, bc_type)
d = beam1.solve()

beam1.plotting(d)