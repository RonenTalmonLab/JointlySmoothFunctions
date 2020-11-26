close all
clear

rng(1);

addpath('./NCCA/');

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
r     = 1/3;
vCosr = cos(2 * pi * vZ + pi);
vSinr = sin(2 * pi * vZ + pi);
vCosR = cos(2 * pi * vEta);
vSinR = sin(2 * pi * vEta);
mY    = [(R + r * vCosr) .* vCosR;
         (R + r * vCosr) .* vSinR;
         r * vSinr];
     
%%
vi   = [-18.1428 40.3022];
vBox = [680 598 623 500];

%%
figure; hold on; grid on; set(gca, 'FontSize', 16);
scatter(mX(1,:), mX(2,:), 30, vZ, 'Fill', 'MarkerEdgeColor', 'none');
xlabel('$x_{1}$',      'Interpreter', 'latex');
ylabel('$x_{2}$',      'Interpreter', 'latex');
axis equal; axis tight

figure('Position', vBox); hold on; grid on; set(gca, 'FontSize', 16);
scatter3(mY(1,:), mY(2,:), mY(3,:), 30, vZ, 'Fill', 'MarkerEdgeColor', 'none');
xlabel('$y_{1}$',      'Interpreter', 'latex');
ylabel('$y_{2}$',      'Interpreter', 'latex');
zlabel('$y_{3}$',      'Interpreter', 'latex');
axis equal; axis tight
view(vi);

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
[mV1,   vEig1] = eig(mK1, 'vector');
[mV2,   vEig2] = eig(mK2, 'vector');
[vEig1, vIdx1]   = sort(vEig1, 'descend');
[~, vIdx2]   = sort(vEig2, 'descend');

%%
d          = round(N / 4);
% d          = 1000;
mW1        = mV1(:,vIdx1(1:d));
mW2        = mV2(:,vIdx2(1:d));

% mW         = [mW1, mW2];
% [mF, S, V] = svd(mW, 'econ');
% [mF, S, V] = svds(mW, 15);
mF         = CommonSVD(mW1, mW2);

%%
N2    = 100;
vZ2   = rand(1, N2);
vEps2 = rand(1, N2);

vRx2 = 1.5 * vEps2 + (vZ2 - 1/2) / 3 + .5;
vPx2 = 2*pi * 2 * vEps2;
mX2  = [vRx2 .* cos(vPx2);
        vRx2 .* sin(vPx2)];
mDx2 = pdist2(mX', mX2')';
mKx2 = exp(-mDx2.^2 / eps1.^2);
mWx2 = (mKx2 * mW1) ./ vEig1(1:d)';

mCx  = mW1' * mF(:,1:d);
mF2  = mWx2 * mCx;

%%
vBox = [300, 300, 1000, 300];
figure('Position', vBox);
for ii = 1 : 3
    subplot(1,3,ii); hold on; 
    scatter(vZ,  mF(:,ii+1),  10, 'b', 'Fill', 'MarkerEdgeColor', 'None');
    scatter(vZ2, mF2(:,ii+1), 10, 'g', 'Fill', 'MarkerEdgeColor', 'None');
    xlabel('$z$',                       'Interpreter', 'latex');
    legend(['$f_{', num2str(ii), '}$'], ['$\tilde{f}_{', num2str(ii), '}$'], 'Interpreter', 'latex', 'Location', 'best');
    axis tight;
    set(gca, 'FontSize', 16);
end


%%
% figure; scatter3(mF(:,1), mF(:,2), mF(:,3),  30,  vZ, 'Fill', 'MarkerEdgeColor', 'none');

%% NCCA
[mU, mV] = NCCA(mX', mY', 16);

%%
figure('Position', vBox);
for ii = 1 : 3
    subplot(1,3,ii); hold on;
    scatter(vZ, mU(:,ii), 10, 'b', 'Fill', 'MarkerEdgeColor', 'None');
    scatter(vZ, mV(:,ii), 10, 'r', 'Fill', 'MarkerEdgeColor', 'None');
    xlabel('$z$',                          'Interpreter', 'latex');
    legend(['$\phi_{', num2str(ii), '}$'], ['$\psi_{', num2str(ii), '}$'], 'Interpreter', 'latex', 'Location', 'best');
    axis tight;
    set(gca, 'FontSize', 16);
end

%%
function U = CommonSVD(W1, W2)
    [Q, L, R] = svd(W1' * W2);
    I         = eye(size(L));
    Sigma     = sqrt([I + L,         zeros(size(L));
                      zeros(size(L)) I - L]);
    U         = 1 / sqrt(2) * [W1, W2] * [Q, Q; R -R] * pinv(Sigma);
    
%     figure; stem(diag(L));
%     figure; stem(diag(Sigma));
end