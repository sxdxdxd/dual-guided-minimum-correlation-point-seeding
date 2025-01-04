function Omega = calculateOmega2D(x, y, u, v)
dx = x(1, 2) - x(1, 1);
dy = y(2, 1) - y(1, 1);
[Ux, Uy] = gradient(u, dx, dy);
[Vx, Vy] = gradient(v, dx, dy);
Q = -1/2.*(Ux.^2 + Vy.^2 + 2.*Uy.*Vx);
A = Ux.^2 + 1/2.*(Uy + Vx).^2 + Vy.^2;
B = 1/2.*(Uy - Vx).^2;
Omega = B ./ (A + B + 1/500 * max(Q, [], 'all'));
end