import numpy

from pysph.base.kernels import CubicSpline
from pysph.base.utils import get_particle_array_rigid_body
from pysph.sph.equation import Group

from pysph.solver.solver import Solver
from pysph.sph.integrator import EPECIntegrator

from pysph.solver.application import Application

from pysph.sph.rigid_body import RigidBodyMoments, RigidBodyMotion, RK2StepRigidBody


# constants

l = 10.
r = 0.5
rho0 = 10.
hdx = 1.

class SimplePendulum(Application):
    def create_particles(self):
        nx, ny = 7,6
        dx = 1.0 / (nx-1)
        x_bob,y_bob,z_bob= numpy.mgrid[-r:r:nx*1j,-(l-r):-(l+r):ny*1j,0:1:1j]
#        x_string,y_string,z_string= numpy.mgrid[0:0:1j*nx,r-l:0:1j*ny,0:1:1j]
#        x = numpy.concatenate((x_bob,x_string),axis=0)
#        y = numpy.concatenate((y_bob,y_string),axis=0)
#        z = numpy.concatenate((z_bob,z_string),axis=0)
        x = x_bob.flat
        y = y_bob.flat
        z = z_bob.flat
        m = numpy.ones_like(x) * dx * dx * rho0
        h = numpy.ones_like(x) * dx *hdx

        pendulum_bob = get_particle_array_rigid_body(name='bob',x=x,y=y,z=z,m=m,h=h)

        # remove particle outside the bob
        indices = []
        for i in range(len(x)):
            if  numpy.sqrt(x[i]**2 + (y[i]+l)**2) - r > 1e-10 and x[i] != 0:
               indices.append(i)
        pendulum_bob.remove_particles(indices)

        pendulum_bob.omega[2] = 1

        return [pendulum_bob]

    def create_solver(self):
        kernel = CubicSpline(dim=2)
        integrator = EPECIntegrator(bob = RK2StepRigidBody())
        solver = Solver(kernel=kernel,dim=2,integrator=integrator,
                        tf=5,dt=5e-3,adaptive_timestep=False)
        solver.set_print_freq(10)
        return solver

    def create_equations(self):
        equations = [
                    Group(equations=[RigidBodyMotion(dest='bob',sources=None)]),
                    Group(equations=[RigidBodyMoments(dest='bob',sources=None)])
                    ]
        return equations

if __name__=='__main__':
    app = SimplePendulum()
    bob = app.create_particles()
    #app.run()
