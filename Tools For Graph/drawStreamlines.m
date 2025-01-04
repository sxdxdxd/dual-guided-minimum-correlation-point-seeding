function f = drawStreamlines(StreamLines, X, ScalarField, MinAlpha, Seed, CLabel, OpenColorBar, AdjustIntegralConstant)
if ndims(X) == 4
    x = X(:, :, :, 1); y = X(:, :, :, 2); z = X(:, :, :, 3);
    Xmin = x(1, 1, 1); Xmax = x(1, end, 1);
    Ymin = y(1, 1, 1); Ymax = y(end, 1, 1);
    Zmin = z(1, 1, 1); Zmax = z(1, 1, end);
    f = draw3Dstreamlines(StreamLines, Xmax, Xmin, Ymax, Ymin, Zmax, Zmin, ScalarField, MinAlpha, Seed, CLabel, OpenColorBar, AdjustIntegralConstant);
else
    x = X(:, :, 1); y = X(:, :, 2);
    f = draw2DstreamlinesWithContour(StreamLines, x, y, ScalarField, Seed, CLabel, OpenColorBar, AdjustIntegralConstant);
end
end