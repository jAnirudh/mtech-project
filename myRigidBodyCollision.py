from pysph.sph.equation import Equation

class myRigidBodyCollision(Equation):
    def __init__(self,dest,sources,Cn=1e-5):
        super(myRigidBodyCollision,self).__init__(dest,sources)
        self.Cn = Cn

    def loop(self,d_idx,s_idx,d_E,s_E,d_nu,s_nu,d_mu,s_mu,d_cm,s_cm,HIJ,RIJ,
             VIJ,d_x,d_y,d_z,s_x,s_y,s_z,d_Fx,d_Fy,d_Fz,d_body_id,s_body_id,
             d_total_mass,s_total_mass,d_check_bit,d_Fmag):
        '''
        '''
        FnIJ = declare('matrix((3,))')
        FtIJ = declare('matrix((3,))')
        EIJ = declare('matrix((3,))')
        ETIJ = declare('matrix((3,))')

        d_bID = declare('int')
        s_bID = declare('int')

        d_bID = d_body_id[d_idx]
        s_bID = s_body_id[s_idx]
   
        if d_bID != s_bID :
            # Calculate Material Constants
            EI, EJ = d_E[d_idx], s_E[s_idx]
            nuI, nuJ = d_nu[d_idx], s_nu[s_idx]
            muI, muJ = d_mu[d_idx], s_mu[s_idx]
            MI = d_total_mass[d_bID]
            MJ = s_total_mass[s_bID]
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
            EIJ[0] = d_cm[3*d_idx + 0] - s_cm[3*s_idx + 0]
            EIJ[1] = d_cm[3*d_idx + 1] - s_cm[3*s_idx + 1]
            EIJ[2] = d_cm[3*d_idx + 2] - s_cm[3*s_idx + 2]
            E2IJ = EIJ[0]**2 + EIJ[1]**2 + EIJ[2]**2 
            # Normalize to obtain unit vector
            EIJ[0] =  EIJ[0] / sqrt(E2IJ)
            EIJ[1] =  EIJ[1] / sqrt(E2IJ)
            EIJ[2] =  EIJ[2] / sqrt(E2IJ)
            # Calculate Rate of Normal Deformation, DELTADOTIJ
            DELTADOTIJ = VIJ[0]*EIJ[0] + VIJ[1]*EIJ[1] + VIJ[2]*EIJ[2]
            # Calculate Normal Contact Force, FnIJ
            FnrIJ = k_n_ij*DELTAIJ**(1.5)                 # mag of Repulsive force
            FndIJ = gamma_n_ij*DELTAIJ**(0.25)*DELTADOTIJ # mag of Damping force
            FnIJ[0] = EIJ[0] * (FnrIJ - FndIJ)
            FnIJ[1] = EIJ[1] * (FnrIJ - FndIJ)
            FnIJ[2] = EIJ[2] * (FnrIJ - FndIJ)
            #print "%s %s" %(HIJ,RIJ)#DELTAIJ#k_n_ij#FnrIJ# - FndIJ

            ####               Calculate Tangential Forces                ####

            # Calculate Tangential Stiffness and Damping Constants
            k_t_ij = 2./7. * k_n_ij
            gamma_t_ij = 2./7. * gamma_n_ij
            # Calculate Coefficient of friction, mu_ij
            mu_f_ij = 0.5*(muI+muJ)
            # Calculate Unit Vector ETIJ
            ETIJ[0] = VIJ[0] - (DELTADOTIJ * EIJ[0])
            ETIJ[1] = VIJ[1] - (DELTADOTIJ * EIJ[1])
            ETIJ[2] = VIJ[2] - (DELTADOTIJ * EIJ[2])
            # Calculate Rate of Tangential Deformation, DELTATDOTIJ
            DELTATDOTIJ = VIJ[0]*ETIJ[0] + VIJ[1]*ETIJ[1] + VIJ[2]*ETIJ[2]
            # Calculate Tangential Contact Force
            FtrIJ = k_t_ij*DELTATIJ                 # mag of repulsive force
            FtdIJ = gamma_t_ij*DELTATIJ*DELTATDOTIJ # mag of damping force
            # Modified Coulomb Friction Force using the Sigmoid Function, Fc
            Fc = mu_f_ij*(FnrIJ-FndIJ)*tanh(8*DELTATDOTIJ)
            # Calculate Tangential Contact Force
            FtIJ[0] = min(Fc,(FtrIJ-FtdIJ)) * ETIJ[0]
            FtIJ[1] = min(Fc,(FtrIJ-FtdIJ)) * ETIJ[1]
            FtIJ[2] = min(Fc,(FtrIJ-FtdIJ)) * ETIJ[2]

            ####               Total Contact Forces                       ####
            Fmag=sqrt((FnIJ[0]+FtIJ[0])**2+(FnIJ[1]+FtIJ[1])**2+\
                      (FnIJ[2]+FtIJ[2])**2)

            #print "======================================================"
            #print "Before checkbit"
            #print d_check_bit[d_bID]
            #print "Fmag = %s" %Fmag
            #print "======================================================"
            if d_check_bit[d_bID] ==  0 :
                #print "Checkbit!"
                #print Fmag
                d_Fmag[d_bID] = Fmag
                d_Fx[d_bID] = FnIJ[0] + FtIJ[0]
                d_Fy[d_bID] = FnIJ[1] + FtIJ[1]
                d_Fz[d_bID] = FnIJ[2] + FtIJ[2]
                d_check_bit[d_bID] = 1
            elif Fmag - d_Fmag[d_bID] > 1e-9:
                #print "Else checkbit"
                d_Fmag[d_bID] = Fmag
                d_Fx[d_bID] = FnIJ[0] + FtIJ[0]
                d_Fy[d_bID] = FnIJ[1] + FtIJ[1]
                d_Fz[d_bID] = FnIJ[2] + FtIJ[2]

            #print "Fmag = %s, Fmag[%s]=%s" %(Fmag,d_bID,d_Fmag[d_bID])
            #print "d_Fx[%s]=%s d_Fy[%s]=%s d_Fz[%s]=%s" %(d_bID,d_Fx[d_bID],d_bID,d_Fy[d_bID],d_bID,d_Fz[d_bID])
            #print "After checkbit"
            #print d_check_bit[d_bID]
            #print "======================================================"
 
    def post_loop(self,d_idx,d_body_id,d_Fx,d_Fy,d_Fz,d_fx,d_fy,d_fz,
                  d_Fmag,d_check_bit):
        d_fx[d_idx] += d_Fx[d_body_id[d_idx]] 
        d_fy[d_idx] += d_Fy[d_body_id[d_idx]]
        d_fz[d_idx] += d_Fz[d_body_id[d_idx]]
        # reset check_bits and force magnitude
        d_check_bit[d_body_id[d_idx]] = 0
        #print "======================================================"
        #print "id=%s fx=%s fy=%s" %(d_body_id[d_idx],d_fx[d_idx],d_fy[d_idx])
        #print "======================================================"
