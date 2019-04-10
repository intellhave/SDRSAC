%% demo.m
%
% Shows a couple of sample registrations using 
% ICP - Iterative Closest Point
%
% Jakob Wilm and Martin Kjer, Technical University of Denmark, 2012

addpath('../../Utils');

close all;
m = 80; % width of grid
n = m^2; % number of points

% [X,Y] = meshgrid(linspace(-2,2,m), linspace(-2,2,m));
% X = reshape(X,1,[]);
% Y = reshape(Y,1,[]);
% Z = sin(X).*cos(Y);
[X, Y, Z] = readPointCloud('../Dataset/bunny/data/bun000.ply', 0.03);
%[X1, Y1, Z1] = readPointCloud('../Dataset/bunny/data/bun090.ply', 0.1);
X1 = X; Y1 = Y; Z1 = Z;

% Create the data point-matrix
D = [X; Y; Z];
M = [X1; Y1; Z1];
n = size(D,2);

% Translation values (a.u.):
Tx = 0.03;
Ty = -0.01;
Tz = 0.01;

% Translation vector
T = [Tx; Ty; Tz];

% Rotation values (rad.):
rx = -0.01;
ry = 0.01;
rz = 0.5;

Rx = [1 0 0;
      0 cos(rx) -sin(rx);
      0 sin(rx) cos(rx)];
  
Ry = [cos(ry) 0 sin(ry);
      0 1 0;
      -sin(ry) 0 cos(ry)];
  
Rz = [cos(rz) -sin(rz) 0;
      sin(rz) cos(rz) 0;
      0 0 1];
% Rotation matrix
R = Rx*Ry*Rz;

% Transform data-matrix plus noise into model-matrix 
M = R * M + repmat(T, 1, size(M,2));
rng(2912673);
M = M + 0.001*randn(3,size(M,2));
idx = find(M(1,:)>-0.07);
M = M(:, idx);
idx = randsample(size(M,2), round(size(M,2)*0.5));
M = M(:, idx);

plot3(D(1,:), D(2,:), D(3,:), 'b.',M(1,:), M(2,:), M(3,:), 'r.');
axis off;
% Add noise to model and data
%D = D + 0.001*randn(3,size(D,2));

%% Run ICP (standard settings)
[Ricp Ticp ER t] = icp(M, D, 1000);

% Transform data-matrix using ICP result
Dicp = Ricp * D + repmat(Ticp, 1, n);

figure;
plot3(Dicp(1,:),Dicp(2,:),Dicp(3,:),'b.',M(1,:),M(2,:),M(3,:), 'r.');
axis off;
% Plot model points blue and transformed points red
% figure;
% subplot(2,2,1);
% plot3(M(1,:),M(2,:),M(3,:),'bo',D(1,:),D(2,:),D(3,:),'r.');
% axis equal;
% xlabel('x'); ylabel('y'); zlabel('z');
% title('Red: z=sin(x)*cos(y), blue: transformed point cloud');
% 
% % Plot the results
% subplot(2,2,2);
% plot3(M(1,:),M(2,:),M(3,:),'bo',Dicp(1,:),Dicp(2,:),Dicp(3,:),'r.');
% axis equal;
% xlabel('x'); ylabel('y'); zlabel('z');
% title('ICP result');
% 
% % Plot RMS curve
% % subplot(2,2,[3 4]);
% % plot(0:15,ER,'--x');
% % xlabel('iteration#');
% % ylabel('d_{RMS}');
% % legend('bruteForce matching');
% % title(['Total elapsed time: ' num2str(t(end),2) ' s']);
% 
% %% Run ICP (fast kDtree matching and extrapolation)
% [Ricp Ticp ER t] = icp(M, D, 15, 'Matching', 'kDtree', 'Extrapolation', true);
% 
% % Transform data-matrix using ICP result
% Dicp = Ricp * D + repmat(Ticp, 1, n);
% 
% figure;
% 
% 
% % Plot model points blue and transformed points red
% % figure;
% % subplot(2,2,1);
% % plot3(M(1,:),M(2,:),M(3,:),'bo',D(1,:),D(2,:),D(3,:),'r.');
% % axis equal;
% % xlabel('x'); ylabel('y'); zlabel('z');
% % title('Red: z=sin(x)*cos(y), blue: transformed point cloud');
% % 
% % % Plot the results
% % subplot(2,2,2);
% % plot3(M(1,:),M(2,:),M(3,:),'bo',Dicp(1,:),Dicp(2,:),Dicp(3,:),'r.');
% % axis equal;
% % xlabel('x'); ylabel('y'); zlabel('z');
% % title('ICP result');
% % 
% % % Plot RMS curve
% % subplot(2,2,[3 4]);
% % plot(0:15,ER,'--x');
% % xlabel('iteration#');
% % ylabel('d_{RMS}');
% % legend('kDtree matching and extrapolation');
% % title(['Total elapsed time: ' num2str(t(end),2) ' s']);
% 
% %% Run ICP (partial data)
% 
% % Partial model point cloud
% Mp = M(:,Y>=0);
% 
% % Boundary of partial model point cloud
% b = (abs(X(Y>=0)) == 2) | (Y(Y>=0) == min(Y(Y>=0))) | (Y(Y>=0) == max(Y(Y>=0)));
% bound = find(b);
% 
% % Partial data point cloud
% Dp = D(:,X>=0);
% 
% [Ricp Ticp ER t] = icp(Mp, Dp, 50, 'EdgeRejection', true, 'Boundary', bound, 'Matching', 'kDtree');
% 
% % Transform data-matrix using ICP result
% Dicp = Ricp * Dp + repmat(Ticp, 1, size(Dp,2));
% 
% % Plot model points blue and transformed points red
% figure;
% subplot(2,2,1);
% plot3(Mp(1,not(b)),Mp(2,not(b)),Mp(3,not(b)),'bo',...
%       Mp(1,b),Mp(2,b),Mp(3,b),'go',...
%       Dp(1,:),Dp(2,:),Dp(3,:),'r.')
% axis equal;
% xlabel('x'); ylabel('y'); zlabel('z');
% title('Red: z=sin(x)*cos(y), blue: transformed point cloud');
% 
% % Plot the results
% subplot(2,2,2);
% plot3(Mp(1,not(b)),Mp(2,not(b)),Mp(3,not(b)),'bo',...
%       Mp(1,b),Mp(2,b),Mp(3,b),'go',...
%       Dicp(1,:),Dicp(2,:),Dicp(3,:),'r.');
% axis equal;
% xlabel('x'); ylabel('y'); zlabel('z');
% title('ICP result');
% 
% % Plot RMS curve
% subplot(2,2,[3 4]);
% plot(0:50,ER,'--x');
% xlabel('iteration#');
% ylabel('d_{RMS}');
% legend('partial overlap');
% title(['Total elapsed time: ' num2str(t(end),2) ' s']);