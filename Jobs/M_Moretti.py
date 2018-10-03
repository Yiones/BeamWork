from monwes.Beam import Beam
from monwes.OpticalElement import Optical_element
from monwes.CompoundOpticalElement import CompoundOpticalElement
from monwes.Shape import BoundaryRectangle
from monwes.SurfaceConic import SurfaceConic
from monwes.Vector import Vector
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.stats import norm
import time
#import h5py
#from srxraylib.plot.gol import plot_scatter

main = "__Moretti_1__"



def theta():

    p = 0.23e-3
    f = 0.23e-3/2
    xp = 0.0127046
    yp = 0.350885
    m = xp/p


    v = Vector(xp, yp-f, 0.)
    v.normalization()
    t = Vector(xp/p, -1, 0.)
    t.normalization()

    print((np.arccos(v.dot(t)))*180./np.pi)

    return np.arccos(v.dot(t))







if main == "__Moretti_1__":

    beam0 = Beam()
    beam0.set_gaussian_divergence(25. / (2 * np.sqrt(2 * np.log(2))) * 1e-6, 25. / (2 * np.sqrt(2 * np.log(2))) * 1e-6)
    beam0.set_rectangular_spot(xmax=1e-6,xmin=-1e-6,zmax=1e-6,zmin=-1e-6)

    xmax = 0.0
    xmin = -0.1
    ymax = 0.150
    ymin = -0.150
    zmax = 0.1
    zmin = -0.0

    bound = BoundaryRectangle(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)
    b = [0]*7
    salpha = [0]*7

    bins_list = [1, 2, 3, 4, 5]
    y_list = [1, 2, 3, 4, 5]


    Nn = 7


    for i in range (Nn):

        print("iteration %d" %i)

        beam = beam0.duplicate()
        alpha =(-0.03 + 0.01 * i) * np.pi / 180
        print(alpha)
        salpha[i] = str(round(-0.03 + 0.01 * i, 3))



        montel =  CompoundOpticalElement.initialize_as_montel_parabolic(p=0.351, q=1., theta_z=theta(), bound1=bound, bound2=bound, angle_of_mismatch=alpha, infinity_location='q')
        beam1, beam2, beam3 = montel.trace_montel(beam,print_footprint=0,mode=1)

        b[i] = beam3[2]


    plt.figure()
    plt.hist(b[0].vx * 1e6, bins=900, normed=1, histtype='step', color='b')
    plt.hist(b[1].vx * 1e6, bins=900, normed=1, histtype='step', color='orange')
    plt.hist(b[2].vx * 1e6, bins=900, normed=1, histtype='step', color='g')
    plt.hist(b[3].vx * 1e6, bins=900, normed=1, histtype='step', color='r')
    plt.hist(b[4].vx * 1e6, bins=900, normed=1, histtype='step', color='darkviolet')
    plt.hist(b[5].vx * 1e6, bins=900, normed=1, histtype='step', color='saddlebrown')
    plt.hist(b[6].vx * 1e6, bins=900, normed=1, histtype='step', color='lightpink')
    plt.legend(salpha)
    plt.xlabel("x'[urad]")
    plt.ylabel("frequency [a.u.]")

    plt.show()