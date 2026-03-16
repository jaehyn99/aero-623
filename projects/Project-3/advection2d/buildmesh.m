function [V, E2N, I2E, h, normals] = buildmesh(L,N);
% builds a mesh of NxN square elements over [0,L]^2

% h = element side length (all elements are squares)
h = L/N;

% V = vertex coordinates
V = zeros((N+1)^2,2);
s = linspace(0,L,N+1);
k = 0;
for j=1:N+1, for i=1:N+1,
    k = k+1;
    V(k,:) = [s(i), s(j)];
end; end;

% E2N = element-to-node array
E2N = zeros(N*N,4);
k = 0;
for j=1:N, for i=1:N,
    k = k+1;
    ii = (j-1)*(N+1) + i;
    E2N(k,:) = [ii, ii+1, ii+N+2, ii+N+1];
end; end;

% I2E = internal edge-to-elem array (all boundaries are periodic)
%       each row is [elemL, elemR, edgeL, edgeR]
% normals = normal vector on each edge, same ordering as I2E
k=0;
ie = 0;
for j=1:N, for i=1:N,
    k = k+1;
    kR = k+1; if (i==N), kR = kR-N; end;
    kU = k+N; if (j==N), kU = k-N*(N-1); end;
    % edge between elements k and kR
    ie = ie+1;
    I2E(ie,:) = [k, kR, 2, 4]; normals(ie,:) = [1,0];
    % edge between elements k and kU
    ie = ie+1;
    I2E(ie,:) = [k, kU, 3, 1]; normals(ie,:) = [0,1];
  end
end
      

