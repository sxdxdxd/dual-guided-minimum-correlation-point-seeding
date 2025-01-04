function Mask = generateMask3D(N_x, N_y, N_z, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, PointsPassedByStreamLine)
[Index_X, Index_Y, Index_Z] = convertAbsoluteCoordinatesToGrid3D(PointsPassedByStreamLine, ...
    N_x, Xmax, Xmin, N_y, Ymax, Ymin, N_z, Zmax, Zmin);
Index = unique(sub2ind([N_y, N_x, N_z], Index_Y, Index_X, Index_Z));
Field = zeros(N_y, N_x, N_z);
Field(Index) = 1;
Mask = logical(Field);
end
