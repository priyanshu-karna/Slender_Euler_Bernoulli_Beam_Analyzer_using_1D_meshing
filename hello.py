import numpy as np 
import matplotlib.pyplot as plt

class Beam:
    def __init__(self, L, E, I, P, load_pos, n_elem, bc_type):
        self.L = L
        self.E = E
        self.I = I
        self.P = P                  
        self.load_pos = load_pos  
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
            
        
        if np.isclose(self.load_pos % self.l_elem, 0) or np.isclose(self.load_pos % self.l_elem, self.l_elem):
            node_idx = int(round(self.load_pos / self.l_elem))
            f[2 * node_idx] = -self.P

        else:
            
            left_node = int(self.load_pos // self.l_elem)
            a = float(self.load_pos % self.l_elem) 
            b = self.l_elem - a
            l = self.l_elem
         
            dof_indices = [2*left_node, 2*left_node+1, 2*left_node+2, 2*left_node+3]
            
            
            f_equivalent = np.array([
                -self.P * b**2 * (3*a + b) / l**3,      
                -self.P * a * b**2 / l**2,              
                -self.P * a**2 * (3*b + a) / l**3,      
                self.P * a**2 * b / l**2                
            ])
            
            for idx, dof in enumerate(dof_indices):
                f[dof] += f_equivalent[idx]
                
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
        deflection = d[0::2]  # Extracts only vertical displacement DOFs

        plt.figure(figsize=(8, 4))
        plt.plot(x, deflection, 'r._',markersize=4,linewidth=1, label='FEM Point Load Deflection')
        
        plt.axvline(x=self.load_pos, color='blue', linestyle='--', label='Point Load Location')
        plt.title(f'Beam Deflection Profile under Point Load ({self.bc_type.title()})')
        plt.xlabel('Beam Length (m)')
        plt.ylabel('Deflection (m)')
        plt.grid(True)
        plt.legend()
        plt.show()

# --- Corrected Example Usage ---
L = 10   
E = 210e9
I = 1e-6

n_elem = 50

P = 5000           
load_pos = 5.59      

beam1 = Beam(L, E, I, P, load_pos, n_elem, bc_type="fixed_fixed")
d = beam1.solve()

beam1.plotting(d)