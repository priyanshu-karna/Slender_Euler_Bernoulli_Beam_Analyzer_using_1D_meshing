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
    def def_plotting(self, d):
        points = 20
        l = self.l_elem
        x_smooth, defl = [], []
        
        for i in range(self.n_elem):
            x0, x1 = i * l, (i + 1) * l
            v1, t1, v2, t2 = d[2*i : 2*i + 4]
            
            x_local = np.linspace(x0, x1, points)
            for xl in x_local:
                xi = (xl - x0) / l
                v = (1 - 3*xi**2 + 2*xi**3)*v1 + (xi - 2*xi**2 + xi**3)*l*t1 + (3*xi**2 - 2*xi**3)*v2 + (-xi**2 + xi**3)*l*t2
                x_smooth.append(xl)
                defl.append(v)
                
        plt.figure(figsize=(8, 4))
        plt.plot(x_smooth, defl, 'r-', linewidth=2, label='Interpolated deflection using cubic Hermite shape functions')
        plt.axvline(self.load_pos, color='black', linestyle='--', alpha=0.5)
        plt.axhline(0, color='black', linewidth=1.8, linestyle='--')
        plt.ylabel('Deflection (m)')
        plt.xlabel('Length (m)')
        plt.title(f'Interpolated Deflection ({self.bc_type.title()})')
        plt.grid(True)
        plt.legend()
        plt.show()

    def get_reactions(self, d):
        
        K, f = self.assemble()
        
        
        R = np.dot(K, d) - f
        
        print("\n" + "="*35)
        print(f"  REACTIONS FOR {self.bc_type.upper()}  ")
        print("="*35)
        
        # Extract left support (Node 0)
        R_left_force = R[0]
        R_left_moment = R[1]
        print(f"Left Support (x = 0 m):")
        print(f"  Vertical Reaction Force: {R_left_force:12.2f} N")
        print(f"  Bending Moment Reaction: {R_left_moment:12.2f} N*m")
        
        # Extract right support (Last Node)
        R_right_force = R[self.ndof - 2]
        R_right_moment = R[self.ndof - 1]
        print(f"\nRight Support (x = {self.L} m):")
        print(f"  Vertical Reaction Force: {R_right_force:12.2f} N")
        print(f"  Bending Moment Reaction: {R_right_moment:12.2f} N*m")
        print("="*35)
        
        return R
    
    def plot_sfd_bmd(self, d):
        # Arrays to store plotting coordinates
        x_coords = []
        shear_forces = []
        bending_moments = []
        
        # Calculate stiffness factor multiplier
        EI = self.E * self.I
        l = self.l_elem
        
        # Loop through each element to find local internal forces
        for i in range(self.n_elem):
            x0 = i * l
            x1 = (i + 1) * l
            
            # Extract the 4 local DOFs for this element
            v1, theta1, v2, theta2 = d[2*i : 2*i + 4]
            local_d = np.array([v1, theta1, v2, theta2])
            
            # 1. BENDING MOMENT CALCULATION (M = EI * d2v/dx2)
            # Evaluated at left node (xi = 0) and right node (xi = 1) of the element
            M_left = (EI / l**2) * (-6 * v1 - 4 * l * theta1 + 6 * v2 - 2 * l * theta2)
            M_right = (EI / l**2) * (6 * v1 + 2 * l * theta1 - 6 * v2 + 4 * l * theta2)
            
            # 2. SHEAR FORCE CALCULATION (V = EI * d3v/dx3)
            # Constant across this element
            V_element = (EI / l**3) * (-12 * v1 - 6 * l * theta1 + 12 * v2 - 6 * l * theta2)
            
            # Append values for plotting
            # For Shear: We use two points per element to create clean vertical "jumps" at point loads
            x_coords.extend([x0, x1])
            shear_forces.extend([V_element, V_element])
            
            # For Moment: Linear variation between nodes
            bending_moments.extend([M_left, M_right])
            
        # --- Plotting Subplots ---
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        
        # 1. Shear Force Diagram
        ax1.plot(x_coords, np.array(shear_forces) / 1000, 'b-', linewidth=2) # Convert to kN
        ax1.fill_between(x_coords, np.array(shear_forces) / 1000, color='blue', alpha=0.1)
        ax1.axhline(0, color='black', linewidth=1)
        ax1.axvline(self.load_pos, color='red', linestyle='--', alpha=0.7, label='Load Position')
        ax1.set_title('Shear Force Diagram (SFD)')
        ax1.set_ylabel('Shear Force (kN)')
        ax1.grid(True)
        
        # 2. Bending Moment Diagram
        ax2.plot(x_coords, np.array(bending_moments) / 1000, 'g-', linewidth=2) # Convert to kN*m
        ax2.fill_between(x_coords, np.array(bending_moments) / 1000, color='green', alpha=0.1)
        ax2.axhline(0, color='black', linewidth=1)
        ax2.axvline(self.load_pos, color='red', linestyle='--', alpha=0.7)
        ax2.set_title('Bending Moment Diagram (BMD)')
        ax2.set_xlabel('Beam Length (m)')
        ax2.set_ylabel('Bending Moment (kN*m)')
        ax2.grid(True)
        
        plt.tight_layout()
        plt.show()

# --- Corrected Example Usage ---
L = 10   
E = 210e9
I = 1e-6

n_elem = 250

P = 5000           
load_pos = 4.77    

beam1 = Beam(L, E, I, P, load_pos, n_elem, bc_type="simply_supported")
d = beam1.solve()
beam1.def_plotting(d)
R = beam1.get_reactions(d)
beam1.plot_sfd_bmd(d)