function Mask = generateMask2D(N_x, N_y, Xmin, Xmax, Ymin, Ymax, PointsPassedByStreamLine)
[Index_X, Index_Y] = convertAbsoluteCoordinatesToGrid2D(PointsPassedByStreamLine, N_x, Xmax, Xmin, N_y, Ymax, Ymin);
Index = unique(sub2ind([N_y, N_x], Index_Y, Index_X));
Field = zeros(N_y, N_x);
Field(Index) = 1;
Mask = logical(Field);
end
