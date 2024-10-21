#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.style.use( 'publication.sty' )

from gaussxw import gaussxw
from dgUtilities import *

def uExact( x ):
    y = []
    for i in range( x.shape[0] ):
        if ( x[i] > 0.0 ) :
            y.append( np.tanh( ( x[i] - 0.01 ) / 0.01 ) \
                        * np.exp( -( x[i] + 0.01 )**2 / 0.1 ) )
        else:
            y.append( np.tanh( ( x[i] - 0.01 ) / 0.01 ) \
                        * ( 0 + np.exp( -( x[i] - 0.01 )**2 / 4.0 ) ) )

    return np.array( y )

def plotBase( fig, ax ):

    ax.set_xlim( xL - 0.1, xH + 0.1 )

    ax.set_ylim( -1.5, 1.1 )

    for i in range( x.shape[0]+1 ):
        ax.axvline( xL+i*dx, c = 'k', alpha = 0.2 )

    ax.axis( 'off' )

def plotExactSolution( fig, ax, ls1, c1, ls2, c2, ls3, c3, alpha = 1.0 ):

    xl = xL
    xh = xl + dx
    xx = np.linspace( xl, xh, 100 )
    l1 = ax.plot( xx, uExact( xx ), ls1, color = c1, alpha = alpha )

    xl += dx
    xh += dx
    xx = np.linspace( xl, xh, 100 )
    l2 = ax.plot( xx, uExact( xx ), ls2, color = c2, alpha = alpha )

    xl += dx
    xh += dx
    xx = np.linspace( xl, xh, 100 )
    l3 = ax.plot( xx, uExact( xx ), ls3, color = c3, alpha = alpha )

    return# l1, l2, l3

def computeNodalValues( xl ):

    xh = xl + dx
    xC = xl + 0.5 * dx

    xqN  = xC + dx * etaqN
    xqNN = xC + dx * etaqNN

    intU = np.empty( N, np.float64 )

    for i in range( N ):
        intU[i] \
          = np.sum( wqNN * uExact( xqNN ) * Lagrange( etaqNN, etaqN, i ) )

    uq = np.dot( np.linalg.inv( M ), intU )

    return uq

def computeCellAverage( uq, offset ):

    u = uh( etaqNN, uq, etaqN + offset )

    uK = np.sum( wqNN * u )

    return uK

def plotApproximateSolution( fig, ax, ls, c, xl, uq, offset, alpha = 1.0 ):

    xC = xl + 0.5 * dx

    u = uh( eta, uq, etaqN + offset )

    xx = xC + dx * eta
    ax.plot( xx , u * np.ones( xx.shape[0], dtype = np.float64 ), \
             ls, color = c, alpha = alpha )

    return

def plotCellAverage( fig, ax, ls, c, xl, uK, text ):

    xh = xl + dx

    xx = np.linspace( xl, xh, 100 )
    ax.plot( xx , uK * np.ones( xx.shape[0], dtype = np.float64 ), \
             ls, color = c )
    ax.text( xl+0.24, uK-0.15, text, color = c, fontsize = 14 )

    return uK

if ( __name__ == '__main__' ) :

    xL = -1.5
    xH = +1.5
    nX = 3
    x  = np.linspace( xL, xH, nX )
    dx = ( x[-1] - x[0] ) / np.float64( nX )

    x1l = xL
    x1h = x1l + dx
    x2l = x1h
    x2h = x2l + dx
    x3l = x2h
    x3h = x3l + dx

    N = 3
    NN = 100
    etaqN , wqN  = gaussxw( N  )
    etaqNN, wqNN = gaussxw( NN )
    M = ComputeMassMatrix( N, NN )
    eta = np.linspace( -0.5, +0.5, 100 )
    x2 = np.linspace( x2l, x2h, eta.shape[0] )

    u1q = computeNodalValues( x1l )
    u2q = computeNodalValues( x2l )
    u3q = computeNodalValues( x3l )
    u1K = computeCellAverage( u1q, 0.0 )
    u2K = computeCellAverage( u2q, 0.0 )
    u3K = computeCellAverage( u3q, 0.0 )

    u1KK = computeCellAverage( u1q, -1.0 )
    u3KK = computeCellAverage( u3q, +1.0 )

    fnn = 0
    fig, ax = plt.subplots()
    plotBase( fig, ax )
    plotExactSolution( fig, ax, '-', 'k', '-', 'k', '-', 'k' )
    plt.savefig( 'fig.tci_{:s}.png'.format( str( fnn ).zfill(2) ), dpi = 300 )
    plt.show()
    plt.close()

    fnn += 1
    fig, ax = plt.subplots()
    plotBase( fig, ax )
    plotExactSolution( fig, ax, '-', 'k', '-', 'k', '-', 'k', 0.5 )
    plotApproximateSolution( fig, ax, '-' , 'r', x1l, u1q,  0.0 )
    plotApproximateSolution( fig, ax, '-' , 'k', x2l, u2q,  0.0 )
    plotApproximateSolution( fig, ax, '-' , 'b', x3l, u3q,  0.0 )
    plt.savefig( 'fig.tci_{:s}.png'.format( str( fnn ).zfill(2) ), dpi = 300 )
    plt.show()
    plt.close()

    fnn += 1
    fig, ax = plt.subplots()
    plotBase( fig, ax )
    plotApproximateSolution( fig, ax, '-' , 'r', x1l, u1q,  0.0 )
    plotApproximateSolution( fig, ax, '-' , 'k', x2l, u2q,  0.0 )
    plotApproximateSolution( fig, ax, '-' , 'b', x3l, u3q,  0.0 )
    plotCellAverage( fig, ax, '-.', 'r', x1l, u1K, r'$\bar{u}_{1}$' )
    plotCellAverage( fig, ax, '-.', 'k', x2l, u2K, r'$\bar{u}_{0}$' )
    plotCellAverage( fig, ax, '-.', 'b', x3l, u3K, r'$\bar{u}_{2}$' )
    plt.savefig( 'fig.tci_{:s}.png'.format( str( fnn ).zfill(2) ), dpi = 300 )
    plt.show()
    plt.close()

    fnn += 1
    fig, ax = plt.subplots()
    plotBase( fig, ax )
    plotApproximateSolution( fig, ax, '-' , 'r', x1l, u1q,  0.0 )
    plotApproximateSolution( fig, ax, '-' , 'k', x2l, u2q,  0.0 )
    plotApproximateSolution( fig, ax, '-' , 'b', x3l, u3q,  0.0 )
    plotCellAverage( fig, ax, '-.', 'r', x1l, u1K, r'$\bar{u}_{1}$' )
    plotCellAverage( fig, ax, '-.', 'k', x2l, u2K, r'$\bar{u}_{0}$' )
    plotCellAverage( fig, ax, '-.', 'b', x3l, u3K, r'$\bar{u}_{2}$' )
    plotApproximateSolution( fig, ax, '--', 'r', x2l, u1q,  -1.0 )
    plotApproximateSolution( fig, ax, '--', 'b', x2l, u3q,  +1.0 )
    plt.savefig( 'fig.tci_{:s}.png'.format( str( fnn ).zfill(2) ), dpi = 300 )
    plt.show()
    plt.close()

    fnn += 1
    fig, ax = plt.subplots()
    plotBase( fig, ax )
    plotApproximateSolution( fig, ax, '-' , 'r', x1l, u1q,  0.0 )
    plotApproximateSolution( fig, ax, '-' , 'k', x2l, u2q,  0.0, 0.5 )
    plotApproximateSolution( fig, ax, '-' , 'b', x3l, u3q,  0.0 )
    plotApproximateSolution( fig, ax, '--', 'r', x2l, u1q,  -1.0 )
    plotApproximateSolution( fig, ax, '--', 'b', x2l, u3q,  +1.0 )
    plotCellAverage( fig, ax, '-.', 'r', x1l, u1K, r'$\bar{u}_{1}$' )
    plotCellAverage( fig, ax, '-.', 'k', x2l, u2K, r'$\bar{u}_{0}$' )
    plotCellAverage( fig, ax, '-.', 'b', x3l, u3K, r'$\bar{u}_{2}$' )
    plotCellAverage( fig, ax, ':', 'r', x2l, u1KK, r'$\bar{\bar{u}}_{1}$' )
    plotCellAverage( fig, ax, ':', 'b', x2l, u3KK, r'$\bar{\bar{u}}_{2}$' )
    plt.savefig( 'fig.tci_{:s}.png'.format( str( fnn ).zfill(2) ), dpi = 300 )
    plt.show()
    plt.close()

    fnn += 1
    fig, ax = plt.subplots()
    plotBase( fig, ax )
    plotApproximateSolution( fig, ax, '-' , 'r', x1l, u1q,  0.0 )
    plotApproximateSolution( fig, ax, '-' , 'k', x2l, u2q,  0.0, 0.5 )
    plotApproximateSolution( fig, ax, '-' , 'b', x3l, u3q,  0.0 )
    plotApproximateSolution( fig, ax, '--', 'r', x2l, u1q,  -1.0 )
    plotApproximateSolution( fig, ax, '--', 'b', x2l, u3q,  +1.0 )
    plotCellAverage( fig, ax, '-.', 'r', x1l, u1K, r'$\bar{u}_{1}$' )
    plotCellAverage( fig, ax, '-.', 'k', x2l, u2K, r'$\bar{u}_{0}$' )
    plotCellAverage( fig, ax, '-.', 'b', x3l, u3K, r'$\bar{u}_{2}$' )
    plotCellAverage( fig, ax, ':', 'r', x2l, u1KK, r'$\bar{\bar{u}}_{1}$' )
    plotCellAverage( fig, ax, ':', 'b', x2l, u3KK, r'$\bar{\bar{u}}_{2}$' )
    d = 1.0e-2 * abs( u2K - u1KK )
    ax.annotate( text = '', xy = (-0.4,u1KK-d), xytext = (-0.4,u2K+d), \
                 arrowprops = dict( arrowstyle = '<->', color = 'r' ) )
    d = 1.0e-2 * abs( u3KK - u2K )
    ax.annotate( text = '', xy = (+0.4,u2K-d), xytext = (+0.4,u3KK+d), \
                 arrowprops = dict( arrowstyle = '<->', color = 'b' ) )
    plt.savefig( 'fig.tci_{:s}.png'.format( str( fnn ).zfill(2) ), dpi = 300 )
    plt.show()
    plt.close()
