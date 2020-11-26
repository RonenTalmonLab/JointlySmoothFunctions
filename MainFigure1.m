close all
clear

rng(1);

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
vCosr = cos(2 * pi * vZ);
vSinr = sin(2 * pi * vZ);
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
scatter(mX(1,:), mX(2,:), 30, 'b', 'Fill', 'MarkerEdgeColor', 'k');
xlabel('$x_{1}$',      'Interpreter', 'latex');
ylabel('$x_{2}$',      'Interpreter', 'latex');
axis equal; axis tight

figure('Position', vBox); hold on; grid on; set(gca, 'FontSize', 16);
scatter3(mY(1,:), mY(2,:), mY(3,:), 30, 'b', 'Fill', 'MarkerEdgeColor', 'k');
xlabel('$y_{1}$',      'Interpreter', 'latex');
ylabel('$y_{2}$',      'Interpreter', 'latex');
zlabel('$y_{3}$',      'Interpreter', 'latex');
axis equal; axis tight
view(vi);

%%
figure; hold on; grid on; set(gca, 'FontSize', 16);
scatter(mX(1,:), mX(2,:), 30, sin(2*pi*vZ + pi/2), 'Fill', 'MarkerEdgeColor', 'none');
xlabel('$x_{1}$',      'Interpreter', 'latex');
ylabel('$x_{2}$',      'Interpreter', 'latex');
axis equal; axis tight

figure('Position', vBox); hold on; grid on; set(gca, 'FontSize', 16);
scatter3(mY(1,:), mY(2,:), mY(3,:), 30, sin(2*pi*vZ + pi/2), 'Fill', 'MarkerEdgeColor', 'none');
xlabel('$y_{1}$',      'Interpreter', 'latex');
ylabel('$y_{2}$',      'Interpreter', 'latex');
zlabel('$y_{3}$',      'Interpreter', 'latex');
axis equal; axis tight
view(vi);

%%
figure; hold on; grid on; set(gca, 'FontSize', 16);
scatter(mX(1,:), mX(2,:), 30, vEps, 'Fill', 'MarkerEdgeColor', 'none');
xlabel('$x_{1}$',      'Interpreter', 'latex');
ylabel('$x_{2}$',      'Interpreter', 'latex');
axis equal; axis tight

figure('Position', vBox); hold on; grid on; set(gca, 'FontSize', 16);
scatter3(mY(1,:), mY(2,:), mY(3,:), 30, vEps, 'Fill', 'MarkerEdgeColor', 'none');
xlabel('$y_{1}$',      'Interpreter', 'latex');
ylabel('$y_{2}$',      'Interpreter', 'latex');
zlabel('$y_{3}$',      'Interpreter', 'latex');
axis equal; axis tight
view(vi);

%%
figure; hold on; grid on; set(gca, 'FontSize', 16);
scatter(mX(1,:), mX(2,:), 30, sin(2*pi*vEta), 'Fill', 'MarkerEdgeColor', 'none');
xlabel('$x_{1}$',      'Interpreter', 'latex');
ylabel('$x_{2}$',      'Interpreter', 'latex');
axis equal; axis tight

figure('Position', vBox); hold on; grid on; set(gca, 'FontSize', 16);
scatter3(mY(1,:), mY(2,:), mY(3,:), 30, sin(2*pi*vEta), 'Fill', 'MarkerEdgeColor', 'none');
xlabel('$y_{1}$',      'Interpreter', 'latex');
ylabel('$y_{2}$',      'Interpreter', 'latex');
zlabel('$y_{3}$',      'Interpreter', 'latex');
axis equal; axis tight
view(vi);

