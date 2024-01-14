% ------- LIME low-light image enhancement -------

% main LIME function
function [I, Ti, Tf] = LIME(L, vals_lambda, vals_sigma, vals_gamma)
    % max values in R,G,B channels
    Ti = max(L,[],3) + 0.02; % Compute maximum values in RGB channels and add a small constant
    [wx, wy] = computeTextureWeights(Ti, vals_sigma); % Compute texture weights
    % solve linear equation AX=B
    Tf = solveLinearEquation(Ti, wx, wy, vals_lambda); % Solve linear equation for texture adjustment
    % adjust image - colormap intensity (values)
    Tf = imadjust(Tf, [], [], 1/vals_gamma); % Adjust colormap intensity
    I(:,:,1) = L(:,:,1) ./ Tf; % Enhance Red channel
    I(:,:,2) = L(:,:,2) ./ Tf; % Enhance Green channel
    I(:,:,3) = L(:,:,3) ./ Tf; % Enhance Blue channel
end

function [retx, rety] = computeTextureWeights(fin, sigma)
    fx = diff(fin,1,2); % Compute partial derivative along x-axis
    fx = padarray(fx, [0 1 0], 'post');
    fy = diff(fin,1,1); % Compute partial derivative along y-axis
    fy = padarray(fy, [1 0 0], 'post');

    vareps_s = 0.02; % Small constant for stability
    vareps = 0.001; % Small constant for stability
    wto = max(sum(sqrt(fx.^2 + fy.^2),3) / size(fin,3), vareps_s).^(-1); % Compute texture weights for the entire image

    % matrix resulted from convolution
    FBImg = fin;
    for ic = 1 : size(fin,3)
       ksize = bitor(round(5*sigma), 1);
       % create spatial filters for image processing
       g = fspecial('gaussian', [1, ksize], sigma);
       % 2D-convolution operation
       ret = conv2(fin(:,:,ic), g, 'same');
       ret = conv2(ret, g', 'same');
       FBImg(:,:,ic) = ret;
    end

    fbin = FBImg(:,:,ic);

    gfx = diff(fbin,1,2);
    gfx = padarray(gfx, [0 1], 'post');
    gfy = diff(fbin,1,1);
    gfy = padarray(gfy, [1 0], 'post');

    wtbx = max(sum(abs(gfx),3) / size(fin,3), vareps).^(-1);
    wtby = max(sum(abs(gfy),3) / size(fin,3), vareps).^(-1);
    retx = wtbx.*wto;
    rety = wtby.*wto;
    retx(:,end) = 0;
    rety(end,:) = 0;
end

% solution to linear equation AX=B
function OUT = solveLinearEquation(IN, wx, wy, lambda)
    [r, c, ch] = size(IN);
    dx = -lambda * wx(:);
    dy = -lambda * wy(:);
    B(:, 1) = dx;
    B(:, 2) = dy;
    d = [-r,-1];
    % matrix A of dimensions (r*c, r*c) is created
    % with the diagonals given by the columns of B
    A = spdiags(B, d, r*c, r*c);
    % creating matrix w from array dx
    w = padarray(dx, r, 'pre');
    w = w(1 : end-r);
    % creating matrix n from array dy
    n = padarray(dy, 1, 'pre');
    n = n(1 : end-1);
    D = 1 - (dx + w + dy + n);
    % matrix A_new of dimensions (r*c, r*c) is created
    % with the diagonals given by the columns of D
    A_new = spdiags(D, 0, r*c, r*c);
    A = A + A' + A_new;
    % ichol - incomplete Cholesky factorization
    if exist('ichol', 'builtin')
        % michol - incomplete Cholesky decomposition
        L = ichol(A, struct('michol', 'on'));
        OUT = IN;
        for ii = 1 : ch
            tin = IN(:, :, ii);
            % Solve the linear system of equations A * x = b
            % by means of the Preconditioned Conjugate Gradient iterative method
            [tout, flag] = pcg(A, tin(:), 0.1, 100, L, L');
            OUT(:, :, ii) = reshape(tout, r, c);
        end
    else
        OUT = IN;
        for ii = 1 : ch
            tin = IN(:, :, ii);
            tout = A \ tin(:);
            OUT(:, :, ii) = reshape(tout, r, c);
        end
    end
end

