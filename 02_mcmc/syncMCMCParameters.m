function parameters = syncMCMCParameters(parameters, nP)

switch nP
    case {105,106}
        omega = 54:58; omega2 = 85:89;
        slopes = 59:64; slopes2 = 90:95;
        knots = 96:100; knots2 = 101:105;
    case 113
        omega = 57:61; omega2 = 92:96;
        slopes = 62:68; slopes2 = 97:103;
        knots = 104:108; knots2 = 109:113;
    case 99
        omega = 48:52; omega2 = 79:83;
        slopes = 53:58; slopes2 = 84:89;
        knots = 90:94; knots2 = 95:99;
    case 64
        parameters(:, nP+(54:58)) = parameters(:, 54:58);
        parameters(:, 60:64) = repmat(parameters(:, 59), 1, 5);
        parameters(:, nP+(60:64)) = repmat(parameters(:, nP+59), 1, 5);
        return;
    otherwise
        error('MeanageSync:ParameterCount','Unsupported meanage layout.');
end
parameters(:, omega2) = parameters(:, omega);
parameters(:, nP+omega) = parameters(:, omega);
parameters(:, nP+omega2) = parameters(:, omega);
parameters(:, slopes) = repmat(parameters(:, slopes(1)), 1, numel(slopes));
parameters(:, slopes2) = repmat(parameters(:, slopes(1)), 1, numel(slopes2));
parameters(:, nP+slopes) = repmat(parameters(:, nP+slopes(1)), 1, numel(slopes));
parameters(:, nP+slopes2) = repmat(parameters(:, nP+slopes(1)), 1, numel(slopes2));
parameters(:, knots2) = parameters(:, knots);
parameters(:, nP+knots) = parameters(:, knots);
parameters(:, nP+knots2) = parameters(:, knots);
if nP == 106, parameters(:, 2*nP) = parameters(:, nP); end
end
