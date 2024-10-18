#!/usr/bin/env python3

import numpy as np

from gaussxw import gaussxw

def Lagrange( x, xq, i ):

    L = 1.0

    for j in range( xq.shape[0] ):

        if j != i:

            L *= ( x - xq[j] ) / ( xq[i] - xq[j] )

    return L

def Legendre( y, n ):

    L = 0.0

    x = 2.0 * y

    if   ( n == 0 ) :
        try:
            L = np.ones( x.shape[0] )
        except:
            L = 1.0
    elif ( n == 1 ) :
        L = np.sqrt( 3.0 ) * x
    elif ( n == 2 ) :
        L = np.sqrt( 5.0 ) * ( 3.0 * x**2 - 1.0 ) / 2.0
    elif ( n == 3 ) :
        L = np.sqrt( 7.0 ) * ( 5.0 * x**3 - 3.0 * x ) / 2.0

    return L

def uh( x, u_q, xq ):

    N = xq.shape[0]

    u = 0.0
    for i in range( N ):
        u += u_q[i] * Lagrange( x, xq, i )

    return u

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
