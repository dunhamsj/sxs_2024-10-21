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

def computeCellAverage( uq ):

    u = uh( etaqNN, uq, etaqN )

    uK = np.sum( wqNN * u )

    return uK

def plotApproximateSolution( fig, ax, ls, c, xl, uq ):

    xC = xl + 0.5 * dx

    u = uh( eta, uq, etaqN )

    xx = xC + dx * eta
    ax.plot( xx , u * np.ones( xx.shape[0], dtype = np.float64 ), \
             ls, color = c )

    return

def plotCellAverage( fig, ax, ls, c, xl, uK ):

    xh = xl + dx

    xx = np.linspace( xl, xh, 100 )
    ax.plot( xx , uK * np.ones( xx.shape[0], dtype = np.float64 ), \
             ls, color = c )

    return uK

def nodalToModal( N, uq ):

    Cn = np.empty( N, np.float64 )

    for n in range( N ):
        Cn[n] \
          = np.sum( wqNN * uh( etaqNN, uq, etaqN ) * Legendre( etaqNN, n ) )

    return Cn

def modalToNodal( N, Cn ):

    uq = np.empty( N, np.float64 )

    for i in range( N ):
        uq[i] = sum( [ Cn[n] * Legendre( etaqN[i], n ) for n in range( 2 ) ] )

    return uq

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

    N = 10
    NN = 100
    etaqN , wqN  = gaussxw( N  )
    etaqNN, wqNN = gaussxw( NN )
    M = ComputeMassMatrix( N, NN )
    eta = np.linspace( -0.5, +0.5, 100 )
    x2 = np.linspace( x2l, x2h, eta.shape[0] )

    u1q = computeNodalValues( x1l )
    u2q = computeNodalValues( x2l )
    u3q = computeNodalValues( x3l )
    u1K = computeCellAverage( u1q )
    u2K = computeCellAverage( u2q )
    u3K = computeCellAverage( u3q )

    Cn = nodalToModal( N, u2q )
    uM = sum( [ Cn[i] * Legendre( eta, i ) for i in range( 2 ) ] )

    C1t = ( u3K - u2K ) / np.sqrt( 12.0 )
    u2M = Cn[0] + C1t * Legendre( x2, 1 )

    fig, ax = plt.subplots()
    plotBase( fig, ax )
    plotExactSolution( fig, ax, '-', 'k', '-', 'k', '-', 'k' )
    plt.savefig( 'fig.sl_00.png', dpi = 300 )
    plt.show()
    plt.close()

    fig, ax = plt.subplots()
    plotBase( fig, ax)
    plotExactSolution( fig, ax, '-', 'k', '-', 'k', '-', 'k' )
    plotApproximateSolution( fig, ax, '-.', 'r', x1l, u1q )
    plotApproximateSolution( fig, ax, '-.', 'r', x2l, u2q )
    plotApproximateSolution( fig, ax, '-.', 'r', x3l, u3q )
    plt.savefig( 'fig.sl_01.png', dpi = 300 )
    plt.show()
    plt.close()

    fig, ax = plt.subplots()
    plotBase( fig, ax)
    plotApproximateSolution( fig, ax, '-.', 'r', x1l, u1q )
    plotApproximateSolution( fig, ax, '-.', 'm', x2l, u2q )
    plotApproximateSolution( fig, ax, '-.', 'b', x3l, u3q )
    plt.savefig( 'fig.sl_02.png', dpi = 300 )
    plt.show()
    plt.close()

    fig, ax = plt.subplots()
    plotBase( fig, ax)
    plotApproximateSolution( fig, ax, '-.', 'r', x1l, u1q )
    plotApproximateSolution( fig, ax, '-.', 'm', x2l, u2q )
    plotApproximateSolution( fig, ax, '-.', 'b', x3l, u3q )
    plotCellAverage        ( fig, ax, '-' , 'r', x1l, u1K )
    plotCellAverage        ( fig, ax, '-' , 'm', x2l, u2K )
    plotCellAverage        ( fig, ax, '-' , 'b', x3l, u3K )
    plt.savefig( 'fig.sl_03.png', dpi = 300 )
    plt.show()
    plt.close()

    fig, ax = plt.subplots()
    plotBase( fig, ax)
    plotApproximateSolution( fig, ax, '-.', 'm', x2l, u2q )
    plotCellAverage( fig, ax, '-' , 'r', x1l, u1K )
    plotCellAverage( fig, ax, '-' , 'm', x2l, u2K )
    plotCellAverage( fig, ax, '-' , 'b', x3l, u3K )
    ax.plot( x2, uM, 'm-' )
    ax.text( -0.4, u2K+0.05, r'$C_{0}$', color = 'm', fontsize = 14 )
    ax.text( 0.15, u2K+0.2, r'$C_{1}$', color = 'm', fontsize = 14 )
    plt.savefig( 'fig.sl_04.png', dpi = 300 )
    plt.show()
    plt.close()

    fig, ax = plt.subplots()
    plotBase( fig, ax)
    plotCellAverage( fig, ax, '-' , 'r', x1l, u1K )
    plotCellAverage( fig, ax, '-' , 'm', x2l, u2K )
    plotCellAverage( fig, ax, '-' , 'b', x3l, u3K )
    ax.plot( x2, uM, 'm-' )
    ax.plot( [x1l+0.5*dx,x1h+0.5*dx], [u1K,u2K], 'r--' )
    ax.plot( [x2l+0.5*dx,x2h+0.5*dx], [u2K,u3K], 'b--' )
    d = 5.0e-1 * abs( u2K - u1K )
    ax.text( x1h-0.1, u1K+d, r'$a$', color = 'r', fontsize = 14 )
    d = 9.0e-1 * abs( u3K - u2K )
    ax.text( x3l+0.1, u3K-d, r'$b$', color = 'b', fontsize = 14 )
    ax.text( 0.15, u2K+0.2, r'$C_{1}$', color = 'm', fontsize = 14 )
    plt.savefig( 'fig.sl_05.png', dpi = 300 )
    plt.show()
    plt.close()

    fig, ax = plt.subplots()
    plotBase( fig, ax)
    plotCellAverage( fig, ax, '-' , 'r', x1l, u1K )
    plotCellAverage( fig, ax, '-' , 'm', x2l, u2K )
    plotCellAverage( fig, ax, '-' , 'b', x3l, u3K )
    ax.plot( x2, u2M, 'm-' )
    ax.plot( [x1l+0.5*dx,x1h+0.5*dx], [u1K,u2K], 'r--' )
    ax.plot( [x2l+0.5*dx,x2h+0.5*dx], [u2K,u3K], 'b--' )
    d = 5.0e-1 * abs( u2K - u1K )
    ax.text( x1h-0.1, u1K+d, r'$a$', color = 'r', fontsize = 14 )
    d = 9.0e-1 * abs( u3K - u2K )
    ax.text( x3l+0.1, u3K-d, r'$b$', color = 'b', fontsize = 14 )
    ax.text( 0.1, u2K+0.15, r'$\tilde{C}_{1}=b$', color = 'm', fontsize = 14 )
    plt.savefig( 'fig.sl_06.png', dpi = 300 )
    plt.show()
    plt.close()

    Cn[1] = C1t
    Cn[2:] = 0.0
    u2qt = modalToNodal( N, Cn )
    fig, ax = plt.subplots()
    plotBase( fig, ax)
    plotExactSolution( fig, ax, '-', 'k', '-', 'k', '-', 'k', 0.5 )
    plotApproximateSolution( fig, ax, '-.', 'r', x1l, u1q )
    plotApproximateSolution( fig, ax, '-.', 'm', x2l, u2q )
    plotApproximateSolution( fig, ax, '-' , 'm', x2l, u2qt )
    plotApproximateSolution( fig, ax, '-.', 'b', x3l, u3q )
    plt.savefig( 'fig.sl_07.png', dpi = 300 )
    plt.show()
    plt.close()
