'''_create_geometry.py'''

import numpy

from pysph.base.utils import get_particle_array_rigid_body

def make_cylinder_2d(dx):
    '''
    Make a 2D cylinder of Unit diameter Centered about the origin.

    Paramters:
    ----------
        dx: float
            spacing between two consecutive points

    Returns:
    --------
        x, y, z: [array-like]
            arrays of the x, y, z coordinates
    '''
    x,y,z = numpy.mgrid[-0.5:0.5+dx:dx, -0.5:0.5+dx:dx, 0:1:1j]
    x = x.ravel()
    y = y.ravel()
    z = z.ravel()

    indices = []

    for i in range(len(x)):
        if numpy.sqrt(x[i]**2 + y[i]**2) - 0.5 > 1e-9:
            indices.append(i)

    x = numpy.delete(x, indices)
    y = numpy.delete(y, indices)
    z = numpy.delete(z, indices)
   
    return x,y,z


class CollapsingCylinderGeometry():
    def __init__(self, nCylinder_layers = 6, hdx = 0.5):
        self.container_rho = 15.0 # arbitrarily chosen
        self.hdx = hdx
        self.nCylinder_layers = nCylinder_layers

    def _cylinder_mass(self):
        '''Mass of the cylinder
        '''
        rho = 2.7e-3 # kg/cm**3
        r = 0.5      # cm
        l = 9.9      # cm

        return numpy.pi * r**2 * rho  * l 

    def create_particles(self):
        # create Container
        nx, ny = 500, 500
        dx = 1.0 / (nx -1)

        x, y, z = numpy.mgrid[-1:27:nx*1j, -1:27:ny*1j, 0:1:1j]

        interior = ((x > 0) & (x < 26)) & (y > 0) 
        container = numpy.logical_not(interior)
        x = x[container].flat
        y = y[container].flat
        z = z[container].flat
        
        container_m = numpy.ones_like(x) * self.container_rho * dx * dx 
        container_h = numpy.ones_like(x) * self.hdx * dx
        
        container = get_particle_array_rigid_body(name = 'container' , x = x, 
                      y = y , z = z , m = container_m, h = container_h )

        container.total_mass[0] = numpy.sum(container_m)
        
        # Create Cylinder Arrays
        
        r = 0.5
        nx , ny = 15 , 15
        dx = 1.0 / (nx - 1)

        _x, _y, _z = make_sphere_2d(dx)
        _id = numpy.ones_like(_x,dtype=int)
        n_sphere_particles = len(_x)
        
        disp = []
        for layer in range(self.nCylinder_layers):
            yc = layer + r
            if layer % 2 == 0:
                for i in range(1):
                    xc = i + r
                    disp.append((xc, yc, 0.0))
            else:
                for i in range(5):
                    xc = i + 2*r
                    disp.append((xc, yc, 0.0))

        x, y, z, body_id = [], [], [], []

        for i, d in enumerate(disp):
            x.append(_x + d[0])
            y.append(_y + d[1])
            z.append(_z + d[2])
            body_id.append(_id * i )

        x = numpy.concatenate(x)
        y = numpy.concatenate(y)
        z = numpy.concatenate(z)
        body_id = numpy.concatenate(body_id)
        m = numpy.ones_like(x) * self._cylinder_mass()/n_sphere_particles
        h = numpy.ones_like(x) * self.hdx * 1.0 / (n_sphere_particles - 1)

        cylinder = get_particle_array_rigid_body(name='cylinder', x=x, y=y,
                        z=z, m=m, h=h, body_id=body_id )

        return [cylinder, container]

##############################################################################




'''rigid_body_validation.py'''

from pysph.base.kernels import CubicSpline

from pysph.sph.equation import Group
from pysph.sph.integrator import EPECIntegrator

from pysph.sph.rigid_body import (RigidBodyCollision, RigidBodyMoments,
    RigidBodyMotion, BodyForce, RK2StepRigidBody)

from pysph.solver.application import Application
from pysph.solver.solver import Solver
from collapsing_cylinders_2d import CollapsingCylinderGeometry

class CollapsingCylinders(Application):
    def initialize(self,layers = 6):
        self.layers = layers
 
    def create_particles(self):
        geometry = CollapsingCylinderGeometry(nCylinder_layers=self.layers)
        return geometry.create_particles()

    def create_solver(self):
        kernel = CubicSpline(dim=2)
        integrator = EPECIntegrator(cylinder = RK2StepRigidBody())
        solver = Solver(kernel=kernel, dim=2, integrator=integrator, tf=0.52, 
                        dt=5e-4, adaptive_timestep=False )
        solver.set_print_freq(10)
        return solver

    def create_equations(self):
        equations = [ 
            Group(equations=[
                       BodyForce(dest='cylinder', sources=None, gy=-981),
              RigidBodyCollision(dest='cylinder', 
                                 sources=['cylinder', 'container'], 
                                 k=8, d=2.5, eta=0.01, kt=0.1)]),
            Group(equations=[RigidBodyMoments(dest='cylinder',sources=None)]),
            Group(equations=[RigidBodyMotion(dest='cylinder',sources=None)]),]
        return equations

if __name__ == '__main__':
    app = CollapsingCylinders()
    app.run()
