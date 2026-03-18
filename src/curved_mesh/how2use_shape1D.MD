<!-- 1. Evaluate shape functions:
   shape1D(sig, q, &phi)

2. Evaluate derivatives:
   gradient1D(sig, q, &gphi)

3. Compute edge derivatives:
   dxdsig = sum_j x[j] * gphi[j]
   dydsig = sum_j y[j] * gphi[j]

4. Compute Jacobian:
   J = sqrt(dxdsig^2 + dydsig^2)

5. Compute unit tangent:
   tx = dxdsig / J
   ty = dydsig / J

6. Compute normal:
   n = [ty, -tx]   or   n = [-ty, tx]
   choose the sign so the normal is outward -->