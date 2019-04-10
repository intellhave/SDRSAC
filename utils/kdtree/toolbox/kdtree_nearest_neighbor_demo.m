clc, clear, close all;

%% create data and execute query
p = rand( 50, 2 ); % input data
q = rand( 10, 2 ); % query data
tree = kdtree_build( p );
idxs = kdtree_nearest_neighbor(tree,q);

% visualize
close all;
xlim( [0 1] );
ylim( [0 1] );
hold on; axis equal; axis off;
plot(p(:,1), p(:,2), '.b');
plot(q(:,1), q(:,2),'.r');
plot(p(idxs,1), p(idxs,2),'or');
legend('database', 'query', 'query results');
for i=1:size(q,1)
    lx = [q(i,1), p(idxs(i),1)];
    ly = [q(i,2), p(idxs(i),2)];
    line(lx,ly,'color','green');
end
set(gcf,'color','white');
kdtree_delete(tree);

%% Stress test correctedness (assertion on fail)
clc, clear, close all;
N = 100000; % #dataset
D = 5;     % #dims 
M = 10000;   % #test queries
rand('twister',1); %#ok<RAND>
p = rand(N,D);
kd = kdtree_build(p);
for i=1:1
    query = rand(1,D);
    idx = kdtree_nearest_neighbor(kd,query);   
    dists = sum( (p - repmat(query,[N,1])).^2,2 );
    [~,idx2] = min(dists);
    assert( idx == idx2 );
end
kdtree_delete(kd);
