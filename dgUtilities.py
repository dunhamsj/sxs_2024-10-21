#!/usr/bin/env python3

import numpy as np

from gaussxw import gaussxw

def Lagrange( x, xq, i ):

    L = 1.0

    for j in range( xq.shape[0] ):

        if j != i:

            L *= ( x - xq[j] ) / ( xq[i] - xq[j] )

    return L

def Legendre( x, n ):

    L = 1.0

    y = 2.0 * x

    if   ( n == 0 ) :
        L = np.ones( y.shape[0] )
    elif ( n == 1 ) :
        L = np.sqrt( 3.0 ) * y
    elif ( n == 2 ) :
        L = np.sqrt( 5.0 ) * ( 3.0 * y**2 - 1.0 ) / 2.0

    return L

def rhoh( x, rho_q, xq ):

    N = xq.shape[0]

    rho = 0.0
    for i in range( N ):
        rho += rho_q[i] * Lagrange( x, xq, i )

    return rho

def ComputeMassMatrix( N, NN ):

    xqN , wqN  = gaussxw( N  )
    xqNN, wqNN = gaussxw( NN )

    M = np.zeros( (N,N), np.float64 )

    for i in range( N ):
        for j in range( N ):

            Li \
              = np.array( [ Lagrange( xqNN[q], xqN, i ) for q in range( NN ) ] )
            Lj \
              = np.array( [ Lagrange( xqNN[q], xqN, j ) for q in range( NN ) ] )

            M[i,j] = np.sum( wqNN * Li * Lj )

    return M
