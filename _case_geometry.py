'''Contains the create_particles method needed for setting up the Collapsing 
Cylinder Test Case'''

import numpy
from display_particles import plot_particles
from pysph.base.utils import get_particle_array_rigid_body

def make_2dCyl(dx):
    '''Makes a 2D Cylinder of radius 1cm centered about the orgin
    
    Parameters
    ----------
        dx: float
            Grid Spacing
    
    Returns
    -------
        x,y,z: array-like
            Co-ordinates of the 2D cylinder
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
    '''
Parameters
----------
    1.) nCylinder_layers: int [default = 6]
            Number of layers in the cylinder stack. 
            Experimental results are available for 6 layer stacks for 
            validation of the Collision model and for the 6, 12 layer stack for
            the Coupled Fluid-Structure interaction model.
    2.) hdx: float [default = 1.0]
            Smoothing Length associated with the generated particle arrays
    3.) container_rho: float [default = 15.0 kg/cm**3]
            Density of the container. 
            (Value is chosen arbitrarily and may need tuning) 
    '''
    def __init__(self,nCylinder_layers = 6,hdx = 1.0,container_rho=15.0):
        msg = 'Allowed number of layers in the Cylinder stack = 6/12.\n'
        assert(nCylinder_layers in [6,12]), msg+self.__doc__
        self.container_rho = container_rho
        self.hdx = hdx
        self.nCylinder_layers = nCylinder_layers

    def _cylinder_mass(self):
        '''Mass of the cylinder
        '''
        rho = 2.7e-3 # kg/cm**3
        r = 0.5      # cm
        l = 9.9

        return numpy.pi * r**2 * rho # * l 

    def create_particles(self):
        '''Creates the Container and the Cylinder stack particle arrays.

        The Dimensions
        '''
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
        container_h = numpy.ones_like(x) * self.hdx * 5 * dx
        E = numpy.ones_like(x)*30e4
        nu = numpy.ones_like(x)*0.3
        mu = numpy.ones_like(x)*0.45
        cm = numpy.asarray([13.0,26.0/3.0,0.0]*len(x))
        body_id = numpy.ones_like(x,dtype='int32')*33
        constants = {'E':E,'nu':nu,'mu':mu,'cm':cm}

        container = get_particle_array_rigid_body(name='container',x=x,y=y,z=z,
             m=container_m,h=container_h,constants=constants,body_id=body_id)

        container.total_mass[0] = 10 #numpy.sum(container_m)
        
        # Create Cylinder Arrays
        
        r = 0.5
        nx , ny = 25, 25
        dx = 1.0 / (nx - 1)

        _x, _y, _z = make_2dCyl(dx)
        _id = numpy.ones_like(_x,dtype=int)
        n_sphere_particles = len(_x)
        
        disp = []
        for layer in range(self.nCylinder_layers):
            yc = layer + r
            if layer % 2 == 0:
                for i in range(6):
                    xc = i + r
                    disp.append((xc, yc, 0.0))
            else:
                for i in range(5):
                    xc = i + 2*r
                    disp.append((xc, yc, 0.0))
        self.disp = disp
        x, y, z, body_id, cm = [], [], [], [], []

        for i, d in enumerate(disp):
            x.append(_x + d[0])
            y.append(_y + d[1])
            z.append(_z + d[2])
            body_id.append(_id * i )
            cm += list(len(_x)*d)

        x = numpy.concatenate(x)
        y = numpy.concatenate(y)
        z = numpy.concatenate(z)
        body_id = numpy.concatenate(body_id)
        m = numpy.ones_like(x) * self._cylinder_mass()/n_sphere_particles
        h = numpy.ones_like(x) * self.hdx * 1.0 / (n_sphere_particles - 1)
        E = numpy.ones_like(x) * 69e5
        nu = numpy.ones_like(x) * 0.3
        mu = numpy.ones_like(x) * 0.45
        Fx = numpy.zeros(33)
        Fy = numpy.zeros(33)
        Fz = numpy.zeros(33)
        constants = {'E':E, 'nu':nu, 'mu':mu,'cm':cm,'Fx':Fx,'Fy':Fy,'Fz':Fz}
        cylinder = get_particle_array_rigid_body(name='cylinder',x=x,y=y,z=z,
                        m=m,h=h,body_id=body_id,constants=constants)
        cylinder.total_mass = numpy.asarray(len(disp)*[self._cylinder_mass()])
        return [cylinder, container]

if __name__== '__main__':
    app = CollapsingCylinderGeometry()
#    app = CollapsingCylinderGeometry(nCylinder_layers = 8)
    cyl,box = app.create_particles()
#    plot_particles(app.create_particles())
