#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.style.use( 'publication.sty' )

from gaussxw import gaussxw
from dgUtilities import *

def plotBase( fig, ax ):

    ax.set_xlim( etaL - 0.1, etaH + 0.1 )
    ax.set_ylim( etaL - 0.1, etaH + 0.1 )

    ax.plot( [ etaL, etaH ], [ etaL, etaL ], 'k', lw = 2 )
    ax.plot( [ etaL, etaH ], [ etaH, etaH ], 'k', lw = 2 )
    ax.plot( [ etaL, etaL ], [ etaL, etaH ], 'k', lw = 2 )
    ax.plot( [ etaH, etaH ], [ etaL, etaH ], 'k', lw = 2 )

    ax.set_aspect( 'equal' )
    ax.axis( 'off' )

def plotRefinedEdges( fig, ax ):

    ax.plot( [ 0.0, 0.0 ], [ etaL, etaH ], 'b', lw = 2 )
    ax.plot( [ etaL, etaH ], [ 0.0, 0.0 ], 'b', lw = 2 )

def plotVolumePoints( fig, ax, N, ref = True, ms = '.' ):

    for i in range( N ):
        ax.plot( etaq, [ etaq[i] for g in range( N ) ], 'k{:}'.format(ms) )

    if ref:
        for j2 in range( 2 ):
            for j1 in range( 2 ):
                for i2 in range( N ):
                    for i1 in range( N ):
                        ax.plot( etaOfXi( etaq[i1], j1 ), \
                                 etaOfXi( etaq[i2], j2 ), 'b.' )

def plotFacePoints( fig, ax, N ):

#    for i1 in range( N ):
#        ax.plot( etaq[i1], etaL, 'ks' )
#        ax.plot( etaq[i1], etaH, 'ks' )
#    for i2 in range( N ):
#        ax.plot( etaL, etaq[i2], 'ks' )
#        ax.plot( etaH, etaq[i2], 'ks' )
#
#    for j2 in range( 2 ):
#        for j1 in range( 2 ):
#            for i1 in range( N ):
#                ax.plot( etaOfXi( etaq[i1], j1 ), etaOfXi( etaL, j2 ), 'bs' )
#                ax.plot( etaOfXi( etaq[i1], j1 ), etaOfXi( etaH, j2 ), 'bs' )
#            for i2 in range( N ):
#                ax.plot( etaOfXi( etaL, j2 ), etaOfXi( etaq[i2], j1 ), 'bs' )
#                ax.plot( etaOfXi( etaH, j2 ), etaOfXi( etaq[i2], j1 ), 'bs' )

    k = 0
    for j2 in range( 2 ):
        for j1 in range( 2 ):
            k += 1
            ax.text( -0.3 + 0.5 * j1, -0.06 + 0.5 * j2, \
                     r'$j={:}$'.format(k), fontsize = 14 )

    for i2 in range( N ):
        ax.plot( etaH, etaq[i2], 'ks' )
    for j2 in range( 2 ):
        for i2 in range( N ):
            ax.plot( etaOfXi( etaH, 1 ), etaOfXi( etaq[i2], j2 ), 'bs' )

def etaOfXi( xi, j ):
    return 0.5 * ( xi + 0.5 * ( -1 )**(j+1) )

def xiOfEta( eta, j ):
    return 2.0 * eta - 0.5 * ( -1 )**(j+1)

if ( __name__ == '__main__' ) :

    etaL = -0.5
    etaH = +0.5

    N = 2

    etaq, wq  = gaussxw( N )

    fig, ax = plt.subplots()
    plotBase( fig, ax )
    plotRefinedEdges( fig, ax )
    plotVolumePoints( fig, ax, N )
    plt.savefig( 'fig.amrElement.png', dpi = 300 )
    plt.show()
    plt.close()

    fig, ax = plt.subplots()
    plotBase( fig, ax )
    plotRefinedEdges( fig, ax )
    plotVolumePoints( fig, ax, N )
    plotFacePoints( fig, ax, N )
    plt.savefig( 'fig.amrElement_facePoints.png', dpi = 300 )
    plt.show()
    plt.close()

    N = 1002

    for n in range( 1, N ):

        etaq, wq  = gaussxw( n )

        fig, ax = plt.subplots()
        ax.set_title( r'$N = {:}$'.format( n ) )
        plotBase( fig, ax )
        plotVolumePoints( fig, ax, n, False, ',' )
        fn = 'fig.dg{:s}.png'.format(str(n-1).zfill(3))
        plt.savefig( fn )
        print( 'Saved {:}'.format( fn ) )
        #plt.show()
        plt.close()
        del etaq, wq, fig, ax

    exit()
