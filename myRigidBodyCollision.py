from numpy import sqrt,reciprocal, zeros

class myRigidBodyCollision(Equation):
    def __init__(self,dest,sources,Cn=1e-5):
        super(myRigidBodyCollision,self).__init__(dest,sources)
        self.Cn = Cn
        
    def _get_Rstar(self, d_idx, s_idx, d_x, d_y, d_z, s_x, s_y, s_z):
        r_d = sqrt(d_x[d_idx]**2 + d_y[d_idx]**2 + d_z[d_idx]**2)
        r_s = sqrt(s_x[s_idx]**2 + s_y[s_idx]**2 + s_z[s_idx]**2)
        return (r_s*r_d / r_s+r_d)
    
    def loop(self,d_idx,s_idx,d_E,s_E,d_nu,s_nu,s_mu,d_cm,s_cm,HIJ,RIJ,VIJ,
             d_x,d_y,d_z,s_x,s_y,s_z,d_fx,d_fy,d_fz,d_total_mass,s_total_mass,
             d_body_id, s_body_id):
        '''
        '''
        # Figure out which among the (possibly) many bodies we are dealing with
        d_body , s_body = d_body_id[d_idx], s_body_id[s_idx]

        ##### Calculate Normal Forces using the Modified, Non-Linear, #####
        #####                    Hertzian Model                       #####

        # Calculate Parameters
        M_star = reciprocal( 1.0/d_total_mass[d_body] + 1.0/s_total_mass[s_body])
        E_star = reciprocal( (1-d_nu**2)/d_E + (1-s_nu**2)/s_E )
        R_star = self._get_Rstar(d_idx, s_idx, d_x, d_y, d_z, s_x, s_y, s_z)
        # Calculate Normal Stiffness and Damping Constants
        k_n_ij = 4./3. * E_star * sqrt(R_star)
        gamma_n_ij = self.Cn*sqrt(M_star * E_star * sqrt(R_star))    
        # Calculate Particle Overlap
        DELTAIJ = max(0.,HIJ - RIJ)
        # Calculate the Unit Vector between the two Center of Masses
        EIJ = zeros_like(d_cm)
        for component in range(3):
            EIJ[component] = d_cm[component] - s_cm[component]
        # Normalize to obtain unit vector
        EIJ =  EIJ / sqrt(EIJ[0]**2+EIJ[1]**2+EIJ[2]**2)
        # Calculate Rate of Normal Deformation, DELTADOTIJ
        DELTADOTIJ = VIJ[0]*EIJ[0] + VIJ[1]*EIJ[1] + VIJ[2]*EIJ[2]
        # Calculate Normal Contact Force 
        FnIJ = numpy.zeros(3)
        FnrIJ = k_n_ij*DELTAIJ**(3./2.)                # Normal Repulsive force
        FndIJ = gamma_n_ij*DELTAIJ**(1./4.)*DELTADOTIJ # Normal Damping force
        for component in range(3):
            FnIJ[component] = Eij[component] * ( FnrIJ - FndIJ )
        
        ### Calculate Contact Forces using the Modified Coulomb's Law ###
        
        
