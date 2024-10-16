#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )
plt.rcParams.update({
'text.latex.preamble': r'\usepackage{amsfonts}'})
import matplotlib.patches as patches

from gaussxw import gaussxw
from dgUtilities import *

xL = 0.0
xH = 1.0
xExact = np.linspace( xL, xH, 100 )
x      = np.linspace( xL, xH, 5 )
dx     = ( x[-1] - x[0] ) / np.float64( x.shape[0] )

def rhoExact( x ):
    return 1.0 + 0.1 * np.sin( 2.0 * np.pi * x )

def plotBase( fig, ax ):

    ax.plot( xExact, rhoExact( xExact ), 'k-' )
    ax.set_xlabel( r'$x$' )
    ax.set_ylabel( r'$u$' )

    xx = x[0]
    for i in range( x.shape[0]+1 ):
        ax.axvline( xx )
        xx += dx

    ax.set_xlim( -0.1, 1.1 )
    ax.set_ylim( 0.89, 1.11 )

fig, ax = plt.subplots()
plotBase( fig, ax )
plt.title( r'$u\left(x\right) = 1 + 0.1 \, \sin\left(2\,\pi\,x\right)$', \
           fontsize = 18 )
plt.savefig( 'fig.sine.png', dpi = 300 )
#plt.show()
plt.close()

def plotDensity( N, x, fig, ax, vmin, vmax, vmid, c ):

    NN = 10
    etaqN , wqN  = gaussxw( N  )
    etaqNN, wqNN = gaussxw( NN )

    M = ComputeMassMatrix( N, NN )

    intU = np.empty( N, np.float64 )

    xl = xL
    xh = xl + dx

    for iX1 in range( x.shape[0] ):

        xC = xl + 0.5 * dx

        xqN  = xC + dx * etaqN
        xqNN = xC + dx * etaqNN

        for i in range( N ):
            intU[i] \
              = np.sum( wqNN * rhoExact( xqNN ) * Lagrange( etaqNN, etaqN, i ) )

        rho_q = np.dot( np.linalg.inv( M ), intU )

        xx = np.linspace( xl, xh, 10 )
        rho = rhoh( xx, rho_q, xqN )

        ax.plot( xx , rho * np.ones( xx.shape[0], dtype = np.float64 ), \
                 '-', color = c )
        ax.plot( xqN, rho_q, 'o', color = c )

        xl += dx
        xh += dx

    return

vmin = rhoExact( x.min() )
vmax = rhoExact( x.max() )
vmid = 0.5 * ( vmin + vmax )

N = [ 1, 2, 3 ]
c = [ 'r', 'm', 'b' ]

for i in range( len( N ) ):

    fig, ax = plt.subplots()
    plotBase( fig, ax )

    plotDensity( N[i], x, fig, ax, vmin, vmax, vmid, c[i] )
    plt.title( r'$\forall K \in \mathcal{{T}}_{{h}} : u^{{K}}_{{h}} \in \mathbb{{P}}^{:}\left(K\right)$'.format( N[i]-1 ), \
               fontsize = 18, color = c[i] )
    plt.savefig( 'fig.sine_k{:d}.png'.format( N[i]-1 ), dpi = 300 )
    #plt.show()
    plt.close()

    if ( N[i] == 3 ) :

        fig, ax = plt.subplots()
        plotBase( fig, ax )

        rect = patches.Rectangle \
                 ( (0.59,0.895), 0.22, 0.05, \
                   linewidth = 2, edgecolor = 'r', facecolor = 'none' )
        ax.add_patch(rect)

        plotDensity( N[i], x, fig, ax, vmin, vmax, vmid, c[i] )
        plt.title( r'$\forall K \in \mathcal{{T}}_{{h}} : u^{{K}}_{{h}} \in \mathbb{{P}}^{:}\left(K\right)$'.format( N[i]-1 ), \
                   fontsize = 18, color = c[i] )
        plt.savefig( 'fig.sine_k{:d}_box.png'.format( N[i]-1 ), dpi = 300 )
        #plt.show()

