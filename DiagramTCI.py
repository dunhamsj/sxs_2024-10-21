#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family']      = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

N   = 1000
nOS = 10

xL  = -1.5
xLU = -0.5
xUL = +0.5
xU  = +1.5

yL = -1.0
yU = +1.0

cp = 'red'
ck = 'black'
cn = 'blue'

fs = 13 #  Font size

a = 0.6 # Transparency of extrapolated polynomials

fig, ax =  plt.subplots( 1, 1 )
plt.axis('off')

# Define extent of element
ax.set_xlim( xL, xU )
ax.set_ylim( yL, yU )

# Interfaces
yAxis = ax.axvline( xL , c = 'black', lw = 1   )
yAxis = ax.axvline( xLU, c = 'black', lw = 0.5 )
#yAxis = ax.axvline(  0.0, c = 'black', ls = '--' )
yAxis = ax.axvline( xUL, c = 'black', lw = 0.5 )
yAxis = ax.axvline( xU , c = 'black', lw = 1   )

# x
xp = np.linspace( xL , xLU, N )
xk = np.linspace( xLU, xUL, N )
xn = np.linspace( xUL, xU , N )

#ax.text( xL  - 0.1 ,yL - 0.1, r'$x_{K-3/2}$', fontsize = fs )
ax.text( xLU - 0.068, yL - 0.1, r'$x_{L}$', fontsize = fs )
ax.text( xUL - 0.075, yL - 0.1, r'$x_{H}$', fontsize = fs )
#ax.text( xU  - 0.1 ,yL - 0.1, r'$x_{K+3/2}$', fontsize = fs )

### Previous element ###

yp  = 0.25*(xp+0.3)**2+0.5            # Polynomial
ypE = 0.25*(xk+0.3)**2+0.5            # Extrapolated polynomial
Ap  = 0.5 * ( ypE.min() + ypE.max() ) # Cell-average

ax.text( xp[N//2-10*nOS], 0.9*yU, r'$K^{\left(1\right)}$', color = cp )

# Polynomial
ax.plot( xp, yp , color = cp, ls = '-' )
ax.text( xp[N//2+1*nOS], 0.35*(yp.min()+yp.max()), \
         r'$G_{h}^{\left(1\right)}(x)$', color = cp )

# Cell average
ax.plot( xp, 0.5*(yp.min()+yp.max())*np.ones(xp.shape[0]), \
         color = cp, ls = '-', alpha = a )
ax.text( xp[N//2], 0.53*(yp.min()+yp.max()), \
         r'$G_{K^{\left(1\right)}}^{\left(1\right)}$', \
         color = cp )

# Extrapolated polynomial
ax.plot( xk, ypE, color = cp, ls = '--' )

# Average of extrapolated polynomial
ax.plot( xk, Ap*np.ones(xk.shape[0]), color = cp, ls = '--', alpha = a )
ax.text( xk[N-18*nOS], Ap-0.1, r'$G_{K}^{\left(1\right)}$', \
         color = cp )

### Target element ###

# Polynomial in target element
yk = 0.3*(xk-1.2)**2 - 0.5                      # Polynomial
Ak = 0.5 * ( yk.min() + yk.max() )              # Cell-average

ax.text( xk[N//2-10*nOS], 0.9*yU, r'$K$', color = ck )

# Polynomials
ax.plot( xk, yk, color = ck, ls = '-' )
ax.text( xk[N//2-30*nOS], 15.0*(yk.min()+yk.max()), r'$G_{h}(x)$', color = ck )

# Cell average
ax.plot( xk, Ak*np.ones(xk.shape[0]), color = ck, ls = '-', alpha = a )
ax.text( xk[N-18*nOS], Ak+0.05, r'$G_{K}$', \
         color = ck )

### Next element ###

yn  = 0.5*(xn-1.0)**2 - 0.9                     # Polynomial
ynE = 0.5*(xk-1.0)**2 - 0.9                 # Extrapolated polynomial
An = 0.5 * ( ynE[0] + ynE[-1] )                 # Cell-average

ax.text( xn[N//2-10*nOS], 0.9*yU, r'$K^{\left(2\right)}$', color = cn )

# Polynomial
ax.plot( xn, yn , color = cn, ls = '-' )
ax.text( xn[N//2-5*nOS], 0.65*(yn[0]+yn[-1]), \
         r'$G_{h}^{\left(2\right)}(x)$', color = cn )

# Cell average
ax.plot( xn, 0.5*(yn.min()+yn.max())*np.ones(xn.shape[0]), \
         color = cn, ls = '-', alpha = a )
ax.text( xn[N//2-10*nOS], 0.48*(yn.min()+yn.max()), \
         r'$G_{K^{\left(1\right)}}^{\left(1\right)}$', \
         color = cn )

# Extrapolated polynomial
ax.plot( xk, ynE, color = cn, ls = '--' )

# Average of extrapolated polynomial
ax.plot( xk, An*np.ones(xk.shape[0]), color = cn, ls = '--', alpha = a )
ax.text( xk[N-18*nOS], An+0.05, r'$G_{K}^{\left(1\right)}$', \
         color = cn )

# TCI
ax.annotate( text = '', xy = (xk[0]+0.5,Ak), xytext = (xk[0]+0.5,Ap), \
             arrowprops = dict( arrowstyle = '<->', color = cp ) )
#ax.text( xk[0]+0.5, 0.5*(Ak+Ap), \
#         r'$\left|G^{\left(j-1\right)}_{K}-G^{\left(j\right)}_{K}\right|$', \
#         color = cp )
ax.annotate( text = '', xy = (xk[0]+0.5,An), xytext = (xk[0]+0.5,Ak), \
             arrowprops = dict( arrowstyle = '<->', color = cn ) )
#ax.text( xk[0]+0.5, 0.8*(Ak+An), \
#         r'$\left|G^{\left(j\right)}_{K}-G^{\left(j+1\right)}_{K}\right|$', \
#         color = cn )

## Equation
#ax.text( xn[0]+0.05, Ak, \
#r'$I_{K}=\frac{\left|G^{\left(j\right)}_{K}-G^{\left(j-1\right)}_{K}\right|+\left|G^{\left(j\right)}_{K}-G^{\left(j+1\right)}_{K}\right|}{\mathrm{max}_{j}\left(\left|G_{K}^{\left(j\right)}\right|,\left|G_{K^{\left(j-1\right)}}^{\left(j-1\right)}\right|,\left|G_{K^{\left(j+1\right)}}^{\left(j+1\right)}\right|\right)}$' )

#plt.show()
plt.savefig( 'fig.TCI.png', dpi = 300 )
