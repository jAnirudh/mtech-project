import matplotlib.pyplot as plot

def plot_particles(args):
    cylinder, container = args
    xcyl = cylinder.x
    ycyl = cylinder.y
    xcon = container.x
    ycon = container.y

    plot.xlim(-1,27)
    plot.ylim(-1,27)
    plot.plot(xcyl,ycyl,'b.')
    plot.plot(xcon,ycon,'r.')
#    fig = plot.figure()
#    ax = fig.add_subplot(111,projection='3d')
#    ax.scatter(x,y,z)
    return plot.show()

    
