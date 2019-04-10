function h = plotPointCloud(fileName)

    h = figure();
    [X, Y, Z] = readPointCloud(fileName);
    plot3(X, Y, Z, '.');


end