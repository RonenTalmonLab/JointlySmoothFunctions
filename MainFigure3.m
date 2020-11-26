close all
clear

rng(2);

%%
N    = 4000;
vZ   = rand(1, N);
vEps = rand(1, N);
vEta = rand(1, N);

vRx = 1.5 * vEps + (vZ - 1/2) / 3 + .5;
vPx = 2*pi * 2 * vEps;
mX  = [vRx .* cos(vPx);
       vRx .* sin(vPx)];

R     = 1;
r     = .3;
vCosr = cos(2 * pi * vZ + pi);
vSinr = sin(2 * pi * vZ + pi);
vCosR = cos(2 * pi * vEta);
vSinR = sin(2 * pi * vEta);
mY    = [(R + r * vCosr) .* vCosR;
         (R + r * vCosr) .* vSinR;
         r * vSinr];
     
%%
mW1  = squareform( pdist(mX') );
eps1 = .3 * median(mW1(:));
mK1  = exp(-mW1.^2 / eps1^2);

mW2  = squareform( pdist(mY') );
eps2 = .3 * median(mW2(:));
mK2  = exp(-mW2.^2 / eps2^2);

%%
% [mV1, ~] = eig(mK1, d);
% [mV2, ~] = eig(mK2, d);
[mV1, vEig1] = eig(mK1, 'vector');
[mV2, vEig2] = eig(mK2, 'vector');
[~, vIdx1]   = sort(vEig1, 'descend');
[~, vIdx2]   = sort(vEig2, 'descend');

%%
d          = round(N / 4);
d          = 1000;
mW1        = mV1(:,vIdx1(1:d));
mW2        = mV2(:,vIdx2(1:d));

% mW         = [mW1, mW2];
% [mF, S, V] = svd(mW, 'econ');
% [mF, S, V] = svds(mW, 15);
mF         = CommonSVD(mW1, mW2);

%%
vE = sum( (mW1' * mF).^2 );
figure; stem(vE(1:100));

%%
vPerm = randperm(N);
mYY   = mY(:,vPerm);

mW2  = squareform( pdist(mYY') );
eps2 = .3 * median(mW2(:));
mK2  = exp(-mW2.^2 / eps2^2);

[mVV2, vEig2] = eig(mK2, 'vector');
[~, vIdx2]   = sort(vEig2, 'descend');

mWW2        = mVV2(:,vIdx2(1:d));

vG = svds(mW1' * mWW2, 2);
E0 = 1/2 * (1 + vG(2)) % 0.955840065181394


%%
E0B = 1/2 + sqrt(d - 1/2) * sqrt(N - d - 1/2) / (N - 1)

%%
figure('Position', [300, 300, 1450, 800]);
for ii = 1 : 8
    vF = mF(:,ii);
    e  = norm(mW1' * vF)^2;
    subplot(2,4,ii);
    scatter(vZ, vF, 30, 'b', 'Fill', 'MarkerEdgeColor', 'None');
    xlabel('$z$',                          'Interpreter', 'latex');
    legend(['$f_{', num2str(ii-1), '}$'], 'Interpreter', 'latex', 'Location', 'best');
    axis tight;
    str = ['\left\Vert W_{x}^{T}f_{', num2str(ii-1), '}\right\Vert _{2}^{2} = ', num2str(e)];
    if e > E0
        title(['$', str, '>E_0$'], 'Interpreter', 'Latex', 'color', [0, .7, 0]);
    else
        title(['$', str, '<E_0$'], 'Interpreter', 'Latex', 'color', 'r');
    end
    set(gca, 'FontSize', 16);
    if ii == 1
        ylim([-0.1 0.1] + ylim);
    end
end


%%
function U = CommonSVD(W1, W2)
    [Q, L, R] = svd(W1' * W2);
    I         = eye(size(L));
    L         = min(L, 1);
    Sigma     = sqrt([I + L,         zeros(size(L));
                      zeros(size(L)) I - L]);
    U         = 1 / sqrt(2) * [W1, W2] * [Q, Q; R -R] * pinv(Sigma);
    
%     figure; stem(diag(L));
%     figure; stem(diag(Sigma));
end