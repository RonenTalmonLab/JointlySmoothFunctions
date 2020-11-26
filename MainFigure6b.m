close all
clear

rng(8)

%%
M = [-2, 1
     -1, -1];

p1      = .1;
p2      = -.2;
p3      = .2;

g    = @(mX) [mX(1,:) - mX(2,:).^2; ...
              -mX(1,:).^2 + mX(2,:) + 2 * mX(1,:) .* mX(2,:).^2 - mX(2,:).^4];

%%
N  = 2000;
mP = 2 * (rand(3, N) - 1/2);

mX1 = mP(1:2,:);
v1  = mP(1,:) + mP(2,:).^3;
mX2 = g([v1;...
        mP(3,:)]);
     

%%
mW1  = squareform( pdist(mX1') );
eps1 = .3 * median(mW1(:));
mK1  = exp(-mW1.^2 / eps1^2);

mW2  = squareform( pdist(mX2') );
eps2 = .3 * median(mW2(:));
mK2  = exp(-mW2.^2 / eps2^2);

%%
[mV1, vEig1]   = eig(mK1, 'vector');
[mV2, vEig2]   = eig(mK2, 'vector');
[vEig1, vIdx1] = sort(vEig1, 'descend');
[vEig2, vIdx2] = sort(vEig2, 'descend');
mV1            = mV1(:,vIdx1);
mV2            = mV2(:,vIdx2);

%%
% d          = round(N / 4);
d          = 400;
mW1        = mV1(:,1:d);
mW2        = mV2(:,1:d);

% mW         = [mW1, mW2];
% [mF, S, V] = svd(mW, 'econ');
% [mF, S, V] = svds(mW, 15);
mF         = CommonSVD(mW1, mW2);

%%
vBox = [300, 300, 1000, 300];
figure('Position', vBox);
for ii = 1 : 3
    subplot(1,3,ii); scatter(v1, mF(:,ii), 10, 'b', 'Fill', 'MarkerEdgeColor', 'None');
    xlabel('$p_1 + p_2^3$',                          'Interpreter', 'latex');
    legend(['$\psi_{', num2str(ii), '}$'], 'Interpreter', 'latex', 'Location', 'best');
    axis tight;
    set(gca, 'FontSize', 16);
end

%%
v          = linspace(-1.2, 1.2, 45);
[PP1, PP2] = meshgrid(v, v);
mPP        = [PP1(:)'; PP2(:)'; randn(1, 45^2)];


R  = @(mX) reshape(mX, size(PP1));
figure('Position', vBox);
for ii = 1 : 3
    F = scatteredInterpolant(mP(1,:)', mP(2,:)', mF(:,ii));
%     F.Method = 'natural';
    C = F(mPP(1,:)', mPP(2,:)');
    maxVal = max(mF(:,ii));
    minVal = min(mF(:,ii));
    C(C > maxVal) = maxVal;
    C(C < minVal) = minVal;
        
    subplot(1,3,ii); hold on; grid on; set(gca, 'FontSize', 16);
    scatter(mX1(1,:), mX1(2,:), 30, mF(:,ii), 'Fill', 'MarkerEdgeColor', 'none');
    axis tight;
    vAxis = axis();
    contour(R(mPP(1,:)), R(mPP(2,:)), R(C), 'r', 'LineWidth', 2, 'LineStyle', '-');
    % title('$\mathcal{X}$', 'Interpreter', 'latex');
    xlabel('$p_{1}$',      'Interpreter', 'latex');
    axis equal;
    axis(vAxis);
    if ii > 1
        ax = gca;
        ax.YTick = [];
    else
        ylabel('$p_{2}$',      'Interpreter', 'latex');
    end
%     figure; imagesc(R(C)); colorbar;
end

%%
figure('Position', vBox);
for ii = 1 : 3
        
    subplot(1,3,ii); hold on; grid on; set(gca, 'FontSize', 16);
    scatter(mX2(1,:), mX2(2,:), 30, mF(:,ii), 'Fill', 'MarkerEdgeColor', 'none');
%     scatter(mX2(1,:), mX2(2,:), 30, mP(3,:),S 'Fill', 'MarkerEdgeColor', 'none');
    % title('$\mathcal{X}$', 'Interpreter', 'latex');
    xlabel('$x_{1}$',      'Interpreter', 'latex');
%     axis equal;
    axis tight;

    if ii > 1
        ax = gca;
        ax.YTickLabel = [];
    else
        ylabel('$x_{2}$',      'Interpreter', 'latex');
    end
%     figure; imagesc(R(C)); colorbar;
end

%%
% figure; imagesc(R(C)); colorbar;

%%
% R = @(mX) reshape(mX, size(PP1));
% figure;
% for ii = 1 : 15
%     subplot(3,5,ii); hold on; grid on; set(gca, 'FontSize', 16);
%     scatter(mX1(1,:), mX1(2,:), 30, mF(:,ii), 'Fill', 'MarkerEdgeColor', 'none');
%     contour(R(mX1(1,:)), R(mX1(2,:)), R(mF(:,ii)), 'k');
%     % title('$\mathcal{X}$', 'Interpreter', 'latex');
%     xlabel('$x_{1}$',      'Interpreter', 'latex');
%     ylabel('$x_{2}$',      'Interpreter', 'latex');
%     axis equal; axis tight
% end

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