function ThetaDifferenceField = calculateThetaDifference(u, u_r)
assert(ndims(u) == ndims(u_r) && (ndims(u) ==3 || ndims(u) == 4));
if ndims(u) == 4
    if size(u, 4) == 4
        u(:, :, :, 1) = u(:, :, :, 1) ./ mean(abs(u(:, :, :, 1)), 'all');
        u(:, :, :, 2) = u(:, :, :, 2) ./ mean(abs(u(:, :, :, 2)), 'all');
        u(:, :, :, 3) = u(:, :, :, 3) ./ mean(abs(u(:, :, :, 3)), 'all');
        u(:, :, :, 4) = u(:, :, :, 4) ./ mean(abs(u(:, :, :, 4)), 'all');
        u_r(:, :, :, 1) = u_r(:, :, :, 1) ./ mean(abs(u_r(:, :, :, 1)), 'all');
        u_r(:, :, :, 2) = u_r(:, :, :, 2) ./ mean(abs(u_r(:, :, :, 2)), 'all');
        u_r(:, :, :, 3) = u_r(:, :, :, 3) ./ mean(abs(u_r(:, :, :, 3)), 'all');
        u_r(:, :, :, 4) = u_r(:, :, :, 4) ./ mean(abs(u_r(:, :, :, 4)), 'all');
    end
    VecLength = sqrt(sum(u.^2, 4));
    DotField = sum(u.*u_r, 4) ./ VecLength ./ sqrt(sum(u_r.^2, 4));
else
    if size(u, 3) == 3
        u(:, :, 1) = u(:, :, 1) ./ mean(abs(u(:, :, 1)), 'all');
        u(:, :, 2) = u(:, :, 2) ./ mean(abs(u(:, :, 2)), 'all');
        u(:, :, 3) = u(:, :, 3) ./ mean(abs(u(:, :, 3)), 'all');
        u_r(:, :, 1) = u_r(:, :, 1) ./ mean(abs(u_r(:, :, 1)), 'all');
        u_r(:, :, 2) = u_r(:, :, 2) ./ mean(abs(u_r(:, :, 2)), 'all');
        u_r(:, :, 3) = u_r(:, :, 3) ./ mean(abs(u_r(:, :, 3)), 'all');
    end
    VecLength = sqrt(sum(u.^2, 3));
    DotField = sum(u.*u_r, 3) ./ VecLength ./ sqrt(sum(u_r.^2, 3));
end
DotField(DotField < -1) = -1;
DotField(DotField > 1) = 1;
TehtaField = acos(DotField);
TehtaField(isnan(DotField)) = pi;
TehtaField(VecLength == 0) = 0;
ThetaDifferenceField = TehtaField;
end