close all
clear

%%
M = [-2, 1
     -1, -1];

g = @(x1, x2) [x1-x2^2;...
               -x1^2+x2+2*x1*x2^2-x2^4];
           
g2 = @(x1, x2) [x1+x1^4+2*x1^2*x2+x2^2
               x1^2+x2];
   
y = g(-4,2);
g2(y(1),y(2))


if 1
    p1      = .1;
    p2      = -.2;
    p3      = .2;
else
    p1      = .2;
    p2      = .3;
    p3      = -.1;
end

g(p1 + p2, p3)

%%
g    = @(mX) [mX(1,:) - mX(2,:).^2; ...
              -mX(1,:).^2 + mX(2,:) + 2 * mX(1,:) .* mX(2,:).^2 - mX(2,:).^4];
       
gInv = @(mX) [mX(1,:) + mX(1,:).^4 + 2*mX(1,:).^2.*mX(2,:) + mX(2,:).^2
              mX(1,:).^2 + mX(2,:)];

          
          
mX = 1 * (rand(2, 1000) - 1/2);
% mX = [1;1]
mY = g(mX);
mZ = gInv(mY);

%%
Nt         = 301;
t          = linspace(0, 2, Nt)';

figure; hold on; grid on; set(gca, 'FontSize', 22);
Nv         = 7;
v          = linspace(-1, 1, Nv);
[XX1, XX2] = meshgrid(v, v);
XX         = [XX1(:), XX2(:)]';
N          = size(XX, 2);
for ii = 1 :  N
    vX0 = XX(:,ii);
    if max(abs(vX0)) < .5
        continue
    end
    [t2, mX] = ode45(@(t,x) DxDt(t, x, p1, p2 ,p3), t, vX0);

    plot(mX(:,1), mX(:,2), 'r', 'LineWidth', 1.5);
    arrowh(mX(:,1), mX(:,2), 'r', 200, [10, 35]);
end

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis equal;
axis tight;
axis([-1, 1, -1, 1]);

%=========================================================================%
%%
function Dx = DxDt(~, vX, p1, p2, p3)

    M = [-2, 1
         -1, -1];
               
    gInv = @(vX) [vX(1)+vX(1)^4+2*vX(1)^2*vX(2)+vX(2)^2
                    vX(1)^2+vX(2)];
               
    Jg = @(vX) [1,                    -2*vX(2);
                -2*vX(1) + 2*vX(2)^2, 1 + 4*vX(1)*vX(2) - 4*vX(2)^3];
            
    f  = @(vX) M * [vX(1) - (p1 + p2); ...
                    vX(2) - p3];
    x1 = vX(1);
    x2 = vX(2);
    Dx = Jg(gInv(vX)) * f(gInv(vX)); 
                     
end
      