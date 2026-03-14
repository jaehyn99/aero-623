#include <stdio.h>
#include <stdlib.h>

/******************************************************************/
//   FUNCTION Definition: shapeL1D
int
shapeL1D(double sig, int q, double **pphi)
{
  /* Returns 1D Lagrange shape functions on [0,1]
     for order 1 <= q <= 3 at reference coordinate sig.
     pphi is reallocated to size (q+1). */
    int nshape;
    double *phi;
    
    nshape = q + 1;

    if ((*pphi) == NULL)
        (*pphi) = (double *) malloc(nshape*sizeof(double));
    else
        (*pphi) = (double *) realloc( (void *) (*pphi), nshape*sizeof(double));
    phi = *pphi;
    switch (q){

    // Shape functions for 1D Lagrange polynomials
    case 1:
        phi[0] = 1 - sig;
        phi[1] = sig;
        break;
        
    case 2: 
        phi[0] = 2*sig*sig - 3*sig + 1;
        phi[1] = -4*sig*sig + 4*sig;
        phi[2] = 2*sig*sig - sig;
        break;

    case 3:
        phi[0] = -9/2*sig*sig*sig + 9*sig*sig -11/2*sig +1;
        phi[1] = 27/2*sig*sig*sig - 45/2*sig*sig + 9*sig;
        phi[2] = -27/2*sig*sig*sig + 18*sig*sig - 9/2*sig;
        phi[3] = 9/2*sig*sig*sig - 9/2*sig + sig;
        break;

    default:
        return -1;
    }
    return 0;
}

/******************************************************************/
//   FUNCTION Definition: gradientL1D
int 
gradientL1D(double sig, int q, double **qgphi)
{ 
  /* Returns dphi/dsig for 1D Lagrange shape functions on [0,1]
     for order 1 <= q <= 3 at reference coordinate sig.
     pgphi is reallocated to size (q+1). */

    int nshape;
    double *gphi;

    nshape = q+1;

    if ((*qgphi) == NULL)
        (*qgphi) = (double *) malloc(2*nshape*sizeof(double));
    else
        (*qgphi) = (double *) realloc( (void *) (*qgphi), 2*nshape*sizeof(double));

    gphi = *qgphi;
    
    switch (q){
        case 1:
            gphi[0] = -1;
            gphi[1] = 1;
            break;
                
        case 2:
            gphi[0] = 4*sig - 3;
            gphi[1] = -8*sig + 4;
            gphi[2] = 4*sig - 1;
            break;

        case 3:
            gphi[0] = -27/2*sig*sig + 18*sig - 11/2;
            gphi[1] = 81/2*sig*sig - 45*sig + 9;
            gphi[2] = -81/2*sig*sig + 36*sig - 9/2;
            gphi[3] = 27/2*sig*sig - 9*sig + 1;
            break;

        default:
            return -1;
    }
    return 0;
} // gradientL1D