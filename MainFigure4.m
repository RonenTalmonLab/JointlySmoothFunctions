close all
clear

load Sub_n2

%%
Ns       = length(Kernels);
mPhi{Ns} = [];
for ss = 1 : Ns
    mK               = Kernels{ss};
    [mPhi{ss}, mLam] = eigs(mK, 1000);
end

%%
mV      = cat(2, mPhi{:});
[mF, S] = svd(mV, 'econ');

%%
figure('Position', [300, 300, 1000, 400]); hold on; grid on; set(gca, 'FontSize', 16);
colormap(lines(6));
scatter3(mF(:,2), mF(:,3), mF(:,4), 50, vY, 'Fill', 'MarkerEdgeColor', 'k');
vTicks = 5/12 + (0:5) * 5/6 ;
colorbar('Ticks', vTicks, 'TickLabels', {'Wake', 'REM', 'S1', 'S2', 'S3', 'S4'});
% colorbar('Ticks', vTicks);
% axis equal;
axis tight;
xlabel('$f_1$', 'Interpreter', 'latex');
ylabel('$f_2$', 'Interpreter', 'latex');
zlabel('$f_3$', 'Interpreter', 'latex');

%%
vSmooth = sum((mPhi{1}' * mF).^2);
figure; hold on; grid on; stem(vSmooth(1:100));

%%
mX           = mF(:,1:20);
SvmTemplate  = templateSVM('KernelFunction', 'gaussian', 'kernelScale', 0.02);
SvmMdl       = fitcecoc(mX, vY, 'Learners', SvmTemplate, 'KFold', 10);
SvmMdl.kfoldLoss
1 - SvmMdl.kfoldLoss

%%
vStr = {'Wake', 'REM', 'S1', 'S2', 'S3', 'S4'};
figure('Position', [300, 300, 500, 400]);
cm = confusionchart(vStr(vY+1), vStr(SvmMdl.kfoldPredict+1), ...
                    'ColumnSummary','column-normalized', ...
                    'RowSummary','row-normalized');
cm.Title = ['Total Accuracy = ', num2str(100 * (1 - SvmMdl.kfoldLoss)), '%'];

%% KCCA:
L   = size(Kernels{1}, 1);
Z   = zeros(L);
kI  = 1e-3 * eye(L);
K1  = Kernels{1};
K2  = Kernels{2};
K3  = Kernels{3};
K4  = Kernels{4};
K5  = Kernels{5};
K6  = Kernels{6};
mG1 = [Z,     K1*K2, K1*K3, K1*K4, K1*K5, K1*K6;
       K2*K1, Z,     K2*K3, K2*K4, K2*K5, K2*K6;
       K3*K1, K3*K2, Z,     K3*K4, K3*K5, K3*K6;
       K4*K1, K4*K2, K4*K3, Z,     K4*K5, K4*K6;
       K5*K1, K5*K2, K5*K3, K5*K4, Z,     K5*K6;
       K6*K1, K6*K2, K6*K3, K6*K4, K6*K5, Z];
mG1 = (mG1 + mG1') / 2;
   
mG2 = blkdiag(K1+kI, K2+kI, K3+kI, K4+kI, K5+kI, K6+kI);

%%
[VV, LL]   = eigs(mG1, mG2, 20, 'largestabs');
[vL, vIdx] = sort(diag(LL), 'descend');
VV         = VV(:,vIdx);

%%
VV  = real(VV);
L   = length(vIdx);

%%
VV1 = VV(0*L+1:1*L,:);
VV2 = VV(1*L+1:2*L,:);
VV3 = VV(2*L+1:3*L,:);
VV4 = VV(3*L+1:4*L,:);
VV5 = VV(4*L+1:5*L,:);
VV6 = VV(5*L+1:6*L,:);

V   = VV1 + VV2 + VV3 + VV4 + VV5 + VV6;

%%
figure('Position', [300, 300, 1000, 400]); hold on; grid on; set(gca, 'FontSize', 16);
colormap(lines(6));
scatter3(-V(:,2), V(:,3), V(:,4), 50, vY, 'Fill', 'MarkerEdgeColor', 'k');
vTicks = 5/12 + (0:5) * 5/6 ;
colorbar('Ticks', vTicks, 'TickLabels', {'Wake', 'REM', 'S1', 'S2', 'S3', 'S4'});
% colorbar('Ticks', vTicks);
% axis equal;
axis tight;
xlabel('$\psi_1$', 'Interpreter', 'latex');
ylabel('$\psi_2$', 'Interpreter', 'latex');
zlabel('$\psi_3$', 'Interpreter', 'latex');

%%
mX           = V(:,1:20);
% mX           = mX + randn(size(mX)) / 1000;
SvmTemplate  = templateSVM('KernelFunction', 'gaussian', 'kernelScale', 0.13);
% SvmTemplate  = templateSVM('KernelFunction', 'polynomial', 'PolynomialOrder', 11);
SvmMdl       = fitcecoc(mX, vY, 'Learners', SvmTemplate, 'KFold', 10);
SvmMdl.kfoldLoss
1 - SvmMdl.kfoldLoss

%%
vStr = {'Wake', 'REM', 'S1', 'S2', 'S3', 'S4'};
figure('Position', [300, 300, 500, 400]);
cm = confusionchart(vStr(vY+1), vStr(SvmMdl.kfoldPredict+1), ...
                    'ColumnSummary','column-normalized', ...
                    'RowSummary','row-normalized');
cm.Title = ['Total Accuracy = ', num2str(100 * (1 - SvmMdl.kfoldLoss)), '%'];