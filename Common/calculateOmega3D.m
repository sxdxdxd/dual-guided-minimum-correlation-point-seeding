function Omega = calculateOmega3D(x, y, z, u, v, w)
dx = x(1, 2, 1) - x(1, 1, 1);
dy = y(2, 1, 1) - y(1, 1, 1);
dz = z(1, 1, 2) - z(1, 1, 1);
[Ux, Uy, Uz] = gradient(u, dx, dy, dz);
[Vx, Vy, Vz] = gradient(v, dx, dy, dz);
[Wx, Wy, Wz] = gradient(w, dx, dy, dz);
Q = -1/2*(Ux.^2 + Vy.^2 + Wz.^2) - Uy.*Vx - Uz.*Wx - Vz.*Wy;
A = Ux.^2 + Vy.^2 + Wz.^2 + 1/2.*(Uy + Vx).^2 + 1/2.*(Uz + Wx).^2 + 1/2.*(Vz + Wy).^2;
B = 1/2.*((Uy - Vx).^2 + (Uz - Wx).^2 + (Vz - Wy).^2);
Omega = B ./ (A + B + 1/500*max(Q, [], 'all'));
end