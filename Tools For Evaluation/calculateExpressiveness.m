function ExpressionLevel = calculateExpressiveness(alpha, ThetaDifferenceField, OmegaField, ScalarField_R, ScalarField, handlingNegative)
if min(ScalarField, [], 'all') < 0
    if strcmp(handlingNegative, 'AddConst')
        MinS = abs(min(ScalarField, [], 'all'));
        ScalarField = ScalarField + MinS;
        ScalarField_R = ScalarField_R + MinS;
    else %strcmp(Options.handlingNegative,'Abs')
        ScalarField = abs(ScalarField);
    end
end
if alpha ~= 1
    ExpressionLevelField_V = (1 - ThetaDifferenceField ./ pi) .* OmegaField;
    ExpressionLevel_V = sum(ExpressionLevelField_V, 'all') / sum(OmegaField, 'all');
end
if alpha ~= 0

    D_S = abs(ScalarField - ScalarField_R);
    D_S(isnan(D_S)) = 1;
    D_S = reshape(normalize(D_S(:), 1, 'range'), size(ScalarField)).* ScalarField;
    ExpressionLevel_Phi = 1 - sum(D_S, 'all') / sum(ScalarField, 'all');
end
if alpha == 0
    ExpressionLevel = ExpressionLevel_V;
elseif alpha == 1
    ExpressionLevel = ExpressionLevel_Phi;
else
    ExpressionLevel = (1 - alpha)*ExpressionLevel_V + alpha*ExpressionLevel_Phi;
end
end
