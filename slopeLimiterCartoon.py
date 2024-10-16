#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.style.use( 'publication.sty' )

from gaussxw import gaussxw
from dgUtilities import *

xL = 0.4
xH = 1.0
xExact = np.linspace( xL, xH, 100 )
x      = np.linspace( xL, xH, 3 )
dx     = ( x[-1] - x[0] ) / np.float64( x.shape[0] )

def rhoExact( x ):
    return 1.0 + 0.1 * np.sin( 3.1 * np.pi * x )

def plotBase( fig, ax, vmin, vmax ):

#    ax.set_xlim( xL - 0.1, xH + 0.1 )
#    ax.set_ylim( vmin, vmax )

    ax.axvline( 0.4, c = 'k' )
    ax.axvline( 0.6, c = 'k' )
    ax.axvline( 0.8, c = 'k' )
    ax.axvline( 1.0, c = 'k' )

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
          = np.sum( wqNN * rhoExact( xqNN ) * Lagrange( etaqNN, etaqN, i ) )

    rho_q = np.dot( np.linalg.inv( M ), intU )

    for n in range( N ):
        Cn[n] \
          = np.sum( wqNN * rhoh( etaqNN, rho_q, etaqN ) \
                      * Legendre( etaqNN, n ) )

    return Cn

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
          = np.sum( wqNN * rhoExact( xqNN ) * Lagrange( etaqNN, etaqN, i ) )

    rho_q = np.dot( np.linalg.inv( M ), intU )

    for i in range( N ):
        Cn[i] \
          = np.sum( wqNN * rhoh( xqNN, rho_q, xqN ) * Legendre( etaqNN, i ) )

    xx = np.linspace( xl, xh, 100 )
    rho = rhoh( xx, rho_q, xqN )

    ax.plot( xx , rho * np.ones( xx.shape[0], dtype = np.float64 ), \
             '-', color = c )

    if ( iX1 == 1 ) :
        ax.plot( xx , rhoExact( xx ), '-', color = 'k' )
    else :
        ax.plot( xx , rhoExact( xx ), '-', color = 'k' )

    xl += dx
    xh += dx

    return rho_q

vmin = 0.95 * rhoExact( x ).min()
vmax = 1.04 * rhoExact( x ).max()

fig, ax = plt.subplots()
plotBase( fig, ax, vmin, vmax )
C1m = plotDensity( fig, ax, 'r', 1, 0 )
C1  = plotDensity( fig, ax, 'm', 1, 1 )
C1p = plotDensity( fig, ax, 'b', 1, 2 )
ax.plot( [xL+0.5*dx,xL+1.5*dx], [C1m,C1], 'r--' )
ax.plot( [xL+1.5*dx,xL+2.5*dx], [C1,C1p], 'b--' )

N = 3
Cn = nodalToModal( N, xL+dx )
eta = np.linspace( -0.5, +0.5, 10 )
y = sum( [ Cn[i] * Legendre( eta, i ) for i in range( N ) ] )
xx = np.linspace( xL+dx, xL+2*dx, 10 )
ax.plot( xx, Cn[0] + Cn[1] * Legendre( eta, 1 ), 'm' )

plt.show()
plt.close()
