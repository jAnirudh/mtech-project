from numpy import mgrid,ones_like,sin,cos,sqrt,pi,concatenate

from pysph.base.kernels import CubicSpline
from pysph.base.utils import get_particle_array_rigid_body
from pysph.sph.equation import Group

from pysph.sph.integrator import EPECIntegrator

from pysph.solver.solver import Solver

from pysph.sph.rigid_body import (BodyForce, RigidBodyMotion, 
                                RK2StepRigidBody, RigidBodyMoments)
from pysph.solver.application import Application

class SimplePendulum(Application):
    def initialize(self, rho0=10.0, d=4.0, hdx=1.0, theta=30.0, l = 10.0):
        self.rho0 = rho0
        self.dia = d
        self.hdx = hdx
        self.tr = ( pi * theta ) / 180.0
        self.l = l
        
    def _modify_bob(self,bob):
 #       '''
 #       Recalculates the bob's Position (based on the initial angular 
 #       displacement theta) of the bob particles
 #       
 #       Parameters:
 #       -----------
 #           bob: pySPH Rigid Body ParticleArray
 #               The bob is a standard 2D rigid sphere centered about the origin
 #       Returns:
 #       --------
 #           bob: pySPH Rigid Body ParticleArray
 #               Returns the bob with corrected spatial postion
 #       '''
        x = bob.x*cos(self.tr) - bob.y*sin(self.tr)
        y = bob.x*sin(self.tr) + bob.y*cos(self.tr)        
        bob.x = x
        bob.y = y
        return bob

    def create_particles(self):
        # create the pendulum as a composite rigid body of the string and bob
        nx , ny = 100, 100
        dx = 1.0 / (nx-1)
        r = self.dia / 2
        x_bob,y_bob = mgrid[-r:r:nx*1j,-(self.l-r):-(self.l+r):ny*1j]
        x_string,y_string = mgrid[0:0:nx*1j,r-self.l:0:ny*1j]
        x = concatenate((x_bob,x_string),axis=0)
        y = concatenate((y_bob,y_string),axis=0)         
        x = x.flat
        y = y.flat
        m = ones_like(x) * dx * dx * self.rho0
        h = ones_like(x) * self.hdx * dx
        bob = get_particle_array_rigid_body(name='bob',x=x,y=y,m=m,h=h)
        # indices of particles not comprising of the pendulum
        indices = []
        for i in range(len(x)):
            if  sqrt(x[i]**2 + (y[i]+self.l)**2) - r > 1e-10 and x[i] != 0: 
               indices.append(i)     
        bob.remove_particles(indices)
        bob = self._modify_bob(bob)
        #bob.omega[0] = 0.5
        #bob.omega[1] = 5.0
        bob.vc[0] = 5.0
        bob.vc[1] = 5.0             
        return [bob]
        
    def create_solver(self):
        kernel = CubicSpline(dim=2)
        integrator = EPECIntegrator(bob=RK2StepRigidBody())
        solver = Solver(kernel=kernel, dim=2, integrator=integrator,
                    dt=5e-3, tf=5, adaptive_timestep=False)
        solver.set_print_freq(10)
        return solver

    def create_equations(self):
        equations = []
            #Group(equations=[BodyForce(dest='bob', sources=None, gz=-9.81)]),]
            #Group(equations=[RigidBodyMoments(dest='bob', sources=None)]),
            #Group(equations=[RigidBodyMotion(dest='bob', sources=None)]),]
        return equations

if __name__ == '__main__':
    app = SimplePendulum()
    app.run()
    #bob = app.create_particles()
    #plot_particles(bob)