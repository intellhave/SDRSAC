function plotPointClouds(M, B, markerm, markerb)
    if (nargin < 4)
        markerm = 'r.';
    end
    if (nargin < 3)
        markerb = 'b.';
    end
    %figure;
    plot3(M(1,:), M(2,:), M(3,:), markerm, B(1,:), B(2,:), B(3,:), markerb, 'MarkerSize',0.8);
    axis off;    

end