import numpy

from pysph.base.utils import get_particle_array_rigid_body
from display_particles import plot_particles

def make_sphere_2d(dx):
    '''
    Make a unit sphere Centered about the origin. 
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
    def __init__(self, nCylinder_layers = 6, hdx = 1.0):
        self.container_mass = 15.0 # arbitrarily chosen
        self.hdx = hdx
        self.nCylinder_layers = nCylinder_layers

    def sphere_mass(self,n_particles):
        '''Mass of Each Particle constituting the sphere
        '''
        rho = 2.7e-3 # kg/cm**3
        r = 0.5      # cm
        return 1.0/n_particles * 4.0/3.0 * numpy.pi * r**3 * rho

    def create_particles(self):
        # create Container
        n = 100
        # left face
        lx, ly, lz = numpy.mgrid[0:0:1j*n,0:26:1j*n,0:1:1j]
        # bottom face
        bx, by, bz = numpy.mgrid[0:26:1j*n,0:0:1j*n,0:1:1j]
        # right face
        rx, ry, rz = numpy.mgrid[26:26:1j*n,0:26:1j*n,0:1:1j]
        
        x = numpy.concatenate((lx,bx,rx))
        y = numpy.concatenate((ly,by,ry))
        z = numpy.concatenate((lz,bz,rz))
        
        x = x.flat
        y = y.flat
        z = z.flat
        
        container_m = numpy.ones_like(x) * self.container_mass / (3*n)
        container_hdx = numpy.ones_like(x) * self.hdx / (n-1)
        
        container = get_particle_array_rigid_body(name = 'container' , x = x, 
                      y = y , z = z , m = container_m, hdx = container_hdx )

        container.total_mass[0] = numpy.sum(container_m)
        
        # Create Cylinder Arrays
        
        r = 0.5
        nx , ny = 10 , 10
        dx = 1.0 / (nx - 1)

        _x, _y, _z = make_sphere_2d(dx)
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
        m = numpy.ones_like(x) * self.sphere_mass(n_sphere_particles)
        hdx = numpy.ones_like(x) * self.hdx * 1.0 / (n_sphere_particles - 1)

        cylinder = get_particle_array_rigid_body(name='cylinder', x=x, y=y,
                        z=z, m=m, hdx=hdx, body_id=body_id )

        return [cylinder, container]

if __name__== '__main__':
    app = CollapsingCylinderGeometry()
#    app = CollapsingCylinderGeometry(nCylinder_layers = 8)
    plot_particles(app.create_particles())
