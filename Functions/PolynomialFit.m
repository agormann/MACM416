function polyCoeff = PolynomialFit(xData, yData, maxDegree)
%POLYFIT Returns best-fit polynomial coefficients of maximal degree.
%   Uses the adjusted coefficient of determination (R^2) and a sanity check
%   (to avoid too much oscillation) to determine the polynomial of
%   best-fit.

% CLEANING INPUT
if isrow(xData)
    xData = xData';
end
if isrow(yData)
    yData = yData';
end

% GLOBALS
n = length(xData);
N = 0 : maxDegree;
xx = linspace(min(xData), max(xData));
minY = min(yData);  % required to detect oscillation
maxY = max(yData);
TSS = sum( (yData - mean(yData)).^2 );  % total sum of squares

% BUILDING MODELS
A = xData .^ N; % constructing all vandermonde matrices
models = cell(maxDegree-1, 2);

for degree = 1 : maxDegree-1
    V = A(:, 1:degree+1); % current matrix
    [Q, R] = qr(V);     % QR decomposition for least-squares
    polyCoeff = flipud( R \ (Q'*yData) ); % least-squares solve
    
    predictedY = polyval(polyCoeff, xData);
    SSE = sum( (yData - predictedY).^2 ); % sum of squared errors
    squaredR = 1 - SSE/TSS; % coeff of determination
    adjustedR = 1 - (1-squaredR) * (n-1) / (n-degree-1);

    models{degree, 1} = polyCoeff; % store coeffs
    models{degree, 2} = adjustedR; % store parameter
end

% OPTIMIZATION
bestR = models{1, 2};
index = 1;

for degree = 2 : maxDegree-1
    polyCoeff = models{degree, 1};
    currentR = models{degree, 2};
    
    yy = polyval(polyCoeff, xx);
    minYY = min(yy);
    maxYY = max(yy);
    
    % detecting oscillation
    if minYY <= 0.75*minY || maxYY >= 1.25*maxY
        sanityCheck = 0;
    else
        sanityCheck = 1;
    end
    % verifying adj R^2 improves
    if currentR > bestR && sanityCheck
        bestR = currentR;
        index = degree;
    end
end

% RESULT
polyCoeff = models{index, 1};

end

