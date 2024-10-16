#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.style.use( 'publication.sty' )

from gaussxw import gaussxw
from dgUtilities import *

def uExact( x ):
    return np.tanh((x-0.05)/0.1)

def plotBase( fig, ax ):

    ax.set_xlim( xL - 0.1, xH + 0.1 )

    for i in range( x.shape[0]+1 ):
        ax.axvline( xL+i*dx, c = 'k' )

    ax.axis( 'off' )

def nodalToModal( N, xl ):

    NN = 10
    etaqN , wqN  = gaussxw( N  )
    etaqNN, wqNN = gaussxw( NN )

    M = ComputeMassMatrix( N, NN )

    intU = np.empty( N, np.float64 )
    Cn   = np.empty( N, np.float64 )

    xh = xl + dx
    xC = xl + 0.5 * dx

    xqN  = xC + dx * etaqN
    xqNN = xC + dx * etaqNN

    for i in range( N ):
        intU[i] \
          = np.sum( wqNN * uExact( xqNN ) * Lagrange( etaqNN, etaqN, i ) )

    u_q = np.dot( np.linalg.inv( M ), intU )

    for n in range( N ):
        Cn[n] \
          = np.sum( wqNN * uh( etaqNN, u_q, etaqN ) \
                      * Legendre( etaqNN, n ) )

    return Cn

def plotExactSolution( fig, ax, ls1, c1, ls2, c2, ls3, c3 ):

    xl = xL
    xh = xl + dx
    xx = np.linspace( xl, xh, 100 )
    l1 = ax.plot( xx, uExact( xx ), ls1, color = c1 )

    xl += dx
    xh += dx
    xx = np.linspace( xl, xh, 100 )
    l2 = ax.plot( xx, uExact( xx ), ls2, color = c2 )

    xl += dx
    xh += dx
    xx = np.linspace( xl, xh, 100 )
    l3 = ax.plot( xx, uExact( xx ), ls3, color = c3 )

    return# l1, l2, l3

def plotDensity( fig, ax, c, N, iX1 ):

    NN = 10
    etaqN , wqN  = gaussxw( N  )
    etaqNN, wqNN = gaussxw( NN )

    M = ComputeMassMatrix( N, NN )

    intU = np.empty( N, np.float64 )
    Cn   = np.empty( N, np.float64 )

    xl = x[0] + iX1 * dx
    xh = xl + dx
    xC = xl + 0.5 * dx

    xqN  = xC + dx * etaqN
    xqNN = xC + dx * etaqNN

    for i in range( N ):
        intU[i] \
          = np.sum( wqNN * uExact( xqNN ) * Lagrange( etaqNN, etaqN, i ) )

    u_q = np.dot( np.linalg.inv( M ), intU )

    for i in range( N ):
        Cn[i] \
          = np.sum( wqNN * uh( xqNN, u_q, xqN ) * Legendre( etaqNN, i ) )

    xx = np.linspace( xl, xh, 100 )
    u = uh( xx, u_q, xqN )

    ax.plot( xx , u * np.ones( xx.shape[0], dtype = np.float64 ), \
             '-', color = c )

    xl += dx
    xh += dx

    return u_q

if ( __name__ == '__main__' ) :

    xL = -0.3
    xH = 0.3
    x  = np.linspace( xL, xH, 3 )
    dx = ( x[-1] - x[0] ) / np.float64( x.shape[0] )

    fig, ax = plt.subplots()
    plotBase( fig, ax )
    plotExactSolution( fig, ax, '-', 'k', '-', 'k', '-', 'k' )
    plt.show()
    plt.close()

    fig, ax = plt.subplots()
    plotBase( fig, ax )
    plotExactSolution( fig, ax, '-', 'r', '-', 'm', '-', 'b' )
    plt.show()
    plt.close()

    #C1m = plotDensity( fig, ax, 'r', 1, 0 )
    #C1  = plotDensity( fig, ax, 'm', 1, 1 )
    #C1p = plotDensity( fig, ax, 'b', 1, 2 )
    #ax.plot( [xL+0.5*dx,xL+1.5*dx], [C1m,C1], 'r--' )
    #ax.plot( [xL+1.5*dx,xL+2.5*dx], [C1,C1p], 'b--' )
    #
    #N = 3
    #Cn = nodalToModal( N, xL+dx )
    #eta = np.linspace( -0.5, +0.5, 10 )
    #y = sum( [ Cn[i] * Legendre( eta, i ) for i in range( N ) ] )
    #xx = np.linspace( xL+dx, xL+2*dx, 10 )
    #ax.plot( xx, Cn[0] + Cn[1] * Legendre( eta, 1 ), 'm' )

    plt.show()
    plt.close()
