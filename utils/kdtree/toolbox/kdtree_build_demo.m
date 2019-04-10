%% KDTREE_BUILD_DEMO: performances achieved by the preprocessing
% measure tree construction time
clc, clear, close all;
startValue = 1e3;
doublings = 10;

x = zeros(doublings,1);
y = zeros(doublings,1);

for i=1:doublings
   value = startValue * 2^(i-1);
   p = rand( value, 2 );
   tic;
   tree = kdtree_build( p );
   y(i) = toc;
   x(i) = value;
   fprintf('building tree of %d points took %f seconds\n',value,y(i));
   kdtree_delete( tree );
end


figure; loglog(x, y, '.b'); title('log-log plot')
[coeff, r] = polyfit( log(x), log(y),1);
fprintf('t= %0.2g * n^(%0.2f)\n', exp(coeff(2)), coeff(1))
hold on, plot( x, exp(coeff(2))*x.^(coeff(1)), 'Color', 'red');

%% Test nan values
data = rand(100,4);
data(1:10:end) = nan;
kd = KDTree(data);