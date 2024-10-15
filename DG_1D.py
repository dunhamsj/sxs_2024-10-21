#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )
from scipy.interpolate import interp1d

from gaussxw import gaussxw
from dgUtilities import *

xL = 0.6
xH = 0.8
x = np.linspace( xL, xH, 100 )
xC = 0.5 * ( xL + xH )
dx = xH - xL

def createElement():

    fig, ax = plt.subplots( 1, 1 )

    ax.set_xticks( [ xL, xC, xH ], \
                   [ r'$0.6$', \
                     r'$0.7$', \
                     r'$0.8$' ] )
    ax.axvline( xL )
    ax.axvline( xH )

    return fig, ax

def plotExact( fig, ax ):

    rho = np.empty( (x.shape[0]), np.float64 )
    for i in range( rho.shape[0] ):
        rho[i] = rhoExact( x[i] )

    vmin = rho.min()
    vmax = rho.max()

    ax.plot( x, rho, 'k-', label = 'Exact' )

    return vmin, vmax

def plotDensity( N, NN, fig, ax, vmin, vmax, c ):

    rhoE = interp1d( x, rhoExact( x ) )

    etaqN , wqN  = gaussxw( N  )
    etaqNN, wqNN = gaussxw( NN )

    M = ComputeMassMatrix( N, NN )

    intU = np.empty( N, np.float64 )

    xqN  = xC + dx * etaqN
    xqNN = xC + dx * etaqNN

    for i in range( N ):
        intU[i] \
          = np.sum( wqNN * rhoExact( xqNN ) * Lagrange( etaqNN, etaqN, i ) )

    rho_q = np.dot( np.linalg.inv( M ), intU )

    rho = np.empty( (x.shape[0]), np.float64 )
    for i in range( rho.shape[0] ):
        rho[i] = rhoh( x[i], rho_q, xqN )

    ax.plot( x, rho, c + '--', label = r'$k={:d}$'.format( N-1 ) )
    ax.plot( xqN, rho_q, c + 'o' )

    return

if __name__ == '__main__':

    fig, ax = createElement()
    vmin, vmax = plotExact( fig, ax )

    NN = 10

    N = 1
    c = 'r'
    plotDensity( N, NN, fig, ax, vmin, vmax, c )

    N = 2
    c = 'm'
    plotDensity( N, NN, fig, ax, vmin, vmax, c )

    N = 3
    c = 'b'
    plotDensity( N, NN, fig, ax, vmin, vmax, c )

    ax.legend( loc = 'upper center' )

    #plt.show()

    figName = 'fig.DG_1D.png'
    plt.savefig( figName, dpi = 300 )
    print( '\n  Saved {:}'.format( figName ) )

    plt.close()
