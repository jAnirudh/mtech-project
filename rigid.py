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
                        dt=1e-3, adaptive_timestep=False )
        solver.set_print_freq(10)
        return solver

    def create_equations(self):
        equations = [ 
            Group(equations=[
                       BodyForce(dest='cylinder', sources=None, gy=-981),
              RigidBodyCollision(dest='cylinder', 
                                 sources=['cylinder', 'container'], 
                                 k=1e-2, d=2, eta=0.1, kt=0.1)]),
            Group(equations=[RigidBodyMoments(dest='cylinder',sources=None)]),
            Group(equations=[RigidBodyMotion(dest='cylinder',sources=None)]),]
        return equations

if __name__ == '__main__':
    app = CollapsingCylinders()
    app.run()
