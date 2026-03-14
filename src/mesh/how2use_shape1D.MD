1. Evaluate shape functions:
shape1D(sig, q, &phi)

2. Evaluate derivatives
gradient1D(sig, q, &gphi)

3. Compute edge derivative
dxdsig = sum(lagrange points[j], gphi[j])

4. Compute Jacobian
J = sqrt(dxdsig^2 + dydsig^2)

5. Compute unit tangent
t = dxdssig / J

6. Compute normal (should be outward)
n = [ty, -tx] 