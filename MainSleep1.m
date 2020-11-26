close all
clear

dirPath  = '..\..\cap-sleep-database-1.0.0\cap-sleep-database-1.0.0\';

load PcaCoeffs.mat

%%
subIdx    = 2;
fId       = fopen([dirPath, 'n', num2str(subIdx), '_annotations.txt']);
line      = fgetl(fId);
line      = fgetl(fId);
line      = fgetl(fId);
line      = fgetl(fId);

vDuration = [];
vLabel    = [];
sLabel    = [];
while line ~= -1
    cSplit           = strsplit(line);
    vDuration(end+1) = str2num(cSplit{3});
    sLabel{end+1}    = cSplit{10};
    vLabel(end+1)    = ToLabel(cSplit{10});
    
    line = fgetl(fId);
end

% vTimeIdx  = 500:1000;
vDuration = vDuration(vTimeIdx);
vLabel    = vLabel(vTimeIdx);

figure; plot(vDuration, vLabel);

%%
Fs = 512;
Ts = 1 / Fs;

[hdr, mData] = edfread([dirPath, 'n2.edf']);
% vSensorIdx   = [1, 2, 3, 4, 5, 8];
mData = mData(vSensorIdx,:);

%% Notch
b = [1, -1];
a = [1, -0.9];
mData = filter(b, a, mData, [], 2);

%%
pcaDim      = 10;
[Ns, Nt]    = size(mData);
Nn          = length(vLabel);
Data{Ns,Nn} = [];
Covs{Ns,Nn} = [];
vDataLabel  = [];
vCovIdx     = [];
Nwin        = 30 * Fs;
figure;
for ss = 1 : Ns
    mUs = PcaCoeffs{ss}(:,1:pcaDim);
    for ii = 1 : Nn
        vIdx     = vDuration(ii) + (1 : Nwin) - 1;
        vXi      = mData(ss,vIdx);
        mSi      = spectrogram(vXi, 5 * Fs, 4.5 * Fs, 2^11, Fs, 'Yaxis');
%         spectrogram(vXi, 5 * Fs, 4.5 * Fs, 2^11, Fs, 'Yaxis');
        mXi      = mUs' * abs(mSi);
        n        = size(mXi, 2);

        Covs{ss,ii} = inv(cov(mXi'));
        Data{ss,ii} = mXi;
        if ss == 1
            vDataLabel(end+1:end+n) = vLabel(ii);
            vCovIdx(end+1:end+n)       = ii;
        end
    end
    
    mX = cat(2, Data{ss,:});
    subplot(3,2,ss); scatter3(mX(1,:), mX(2,:), mX(3,:), 50, vDataLabel, 'Fill', 'MarkerEdgeColor', 'k');
    drawnow;
end

%%
N    = size(mX, 2);
vIdx = 2 : 8 : N;
vY   = vDataLabel(vIdx);

Kernels{Ns} = [];
for ss = 1 : Ns
    Ns - ss
    sIdx    = ss;
    mX      = cat(2, Data{sIdx,:});
    mX      = mX(:,vIdx);
    N       = size(mX, 2);
    
    mW   = zeros(N, N);
    for ii = 1 : N
        N - ii
        iIdx = vCovIdx(vIdx(ii));
        vXi  = mX(:,ii);
        mCi  = Covs{sIdx,iIdx};
        for jj = (ii + 1) : N
            jIdx = vCovIdx(vIdx(jj));
            vXj  = mX(:,jj);
            mCj  = Covs{sIdx,jIdx};
            mC   = RM(mCi, mCj);
            
            vD   = vXi - vXj;
            
            mW(ii,jj) = vD' * mC * vD;
            
        end
    end
    
    mW  = sqrt(mW + mW');
    eps = median(mW(:));
    mK  = exp(-mW.^2 / eps^2);
    mK  = (mK + mK') / 2;
    Kernels{ss} = mK;
end

%%
save(['Sub_n', num2str(subIdx)], 'Kernels', 'vY');

%%
function label = ToLabel(state)

    label = NaN;
    if strcmp(state, 'W')  == 1, label = 0; end
    if strcmp(state, 'R')  == 1, label = 1; end
    if strcmp(state, 'S1') == 1, label = 2; end
    if strcmp(state, 'S2') == 1, label = 3; end
    if strcmp(state, 'S3') == 1, label = 4; end
    if strcmp(state, 'S4') == 1, label = 5; end
    

end

%%
function M = RM(A, B)
    
    Symm = @(A) (A + A') / 2;

    A2 = sqrtm(Symm(A));
    AI = inv(A2);
    M  = Symm(A2 * sqrtm(Symm(AI * Symm(B) * AI)) * A2);
%     M  = Symm(A * sqrtm(Symm(A \ B)));

end