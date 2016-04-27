from pysph.sph.equation import Equation

class myRigidBodyCollision(Equation):
    def __init__(self,dest,sources,Cn=1e-5):
        super(myRigidBodyCollision,self).__init__(dest,sources)
        self.Cn = Cn

    def loop(self,d_idx,s_idx,d_E,s_E,d_nu,s_nu,d_mu,s_mu,d_cm,s_cm,HIJ,RIJ,
             VIJ,d_x,d_y,d_z,s_x,s_y,s_z,d_Fx,d_Fy,d_Fz,d_body_id,s_body_id,
             d_total_mass,s_total_mass):
        '''
        '''
        FnIJ = declare('matrix((3,))')
        FtIJ = declare('matrix((3,))')
        EIJ = declare('matrix((3,))')
        ETIJ = declare('matrix((3,))')
        
        if d_body_id[d_idx] != s_body_id[s_idx] :
            # Calculate Material Constants
            EI, EJ = d_E[d_idx], s_E[s_idx]
            nuI, nuJ = d_nu[d_idx], s_nu[s_idx]
            muI, muJ = d_mu[d_idx], s_mu[s_idx]
            MI = d_total_mass[d_body_id[d_idx]]
            MJ = s_total_mass[d_body_id[s_idx]]
            Ri = sqrt(d_x[d_idx]**2 + d_y[d_idx]**2 + d_z[d_idx]**2)
            Rj = sqrt(s_x[s_idx]**2 + s_y[s_idx]**2 + s_z[s_idx]**2)

            ##### Calculate Normal Forces using the Modified, Non-Linear, #####
            #####                    Hertzian Model                       #####

            # Calculate Parameters
            M_star = (MI*MJ)/(MI+MJ)
            E_star = (EI*EJ)/(EJ*(1.0-nuI**2) + EI*(1.0-nuJ**2))
            R_star = (Ri*Rj)/(Ri+Rj)
            # Calculate Normal Stiffness and Damping Constants
            k_n_ij = 4.0/3.0 * E_star * sqrt(R_star)
            gamma_n_ij = self.Cn*sqrt(M_star * E_star * sqrt(R_star))
            # Calculate Particle Overlap
            DELTAIJ = DELTATIJ = max(0.0,HIJ - RIJ)
            # Calculate the Unit Vector between the two Center of Masses
            E2IJ = 0.0
            for i in range(3):
                EIJ[i] = d_cm[3*d_idx + i] - s_cm[3*s_idx + i]
                E2IJ += EIJ[i]**2
            # Normalize to obtain unit vector
            for i in range(3): 
                EIJ[i] =  EIJ[i] / sqrt(E2IJ)
            # Calculate Rate of Normal Deformation, DELTADOTIJ
            DELTADOTIJ = VIJ[0]*EIJ[0] + VIJ[1]*EIJ[1] + VIJ[2]*EIJ[2]
            # Calculate Normal Contact Force, FnIJ
            FnrIJ = k_n_ij*DELTAIJ**(1.5)                 # mag of Repulsive force
            FndIJ = gamma_n_ij*DELTAIJ**(0.25)*DELTADOTIJ # mag of Damping force
            for i in range(3):
                FnIJ[i] = EIJ[i] * (FnrIJ - FndIJ)

            ####               Calculate Tangential Forces                ####

            # Calculate Tangential Stiffness and Damping Constants
            k_t_ij = 2./7. * k_n_ij
            gamma_t_ij = 2./7. * gamma_n_ij
            # Calculate Coefficient of friction, mu_ij
            mu_f_ij = 0.5*(muI+muJ)
            # Calculate Unit Vector ETIJ
            for i in range(3):
                ETIJ[i] = VIJ[i] - (DELTADOTIJ * EIJ[i])
            # Calculate Rate of Tangential Deformation, DELTATDOTIJ
            DELTATDOTIJ = VIJ[0]*ETIJ[0] + VIJ[1]*ETIJ[1] + VIJ[2]*ETIJ[2]
            # Calculate Tangential Contact Force
            FtrIJ = k_t_ij*DELTATIJ                 # mag of repulsive force
            FtdIJ = gamma_t_ij*DELTATIJ*DELTATDOTIJ # mag of damping force
            # Modified Coulomb Friction Force using the Sigmoid Function, Fc
            Fc = mu_f_ij*(FnrIJ-FndIJ)*tanh(8*DELTATDOTIJ)
            # Calculate Tangential Contact Force
            for i in range(3):
                FtIJ[i] = min(Fc,(FtrIJ-FtdIJ)) * ETIJ[i]

            ####               Total Contact Forces                       ####
            d_Fx[d_body_id[d_idx]] += FnIJ[0] + FtIJ[0]
            d_Fy[d_body_id[d_idx]] += FnIJ[1] + FtIJ[1]
            d_Fz[d_body_id[d_idx]] += FnIJ[2] + FtIJ[2]

    def post_loop(self,d_idx,d_body_id,d_Fx,d_Fy,d_Fz,d_fx,d_fy,d_fz):
       d_fx[d_idx] += d_Fx[d_body_id[d_idx]] 
       d_fy[d_idx] += d_Fy[d_body_id[d_idx]]
       d_fz[d_idx] += d_Fz[d_body_id[d_idx]] 
