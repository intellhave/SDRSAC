function genSyntheticData()
    
    addpath('../Utils/toolbox_graph/');    
    addpath('../Utils/sample_icp/');  
    addpath('../Utils/MyCrustOpen070909/');
    addpath('../Utils/');
    
   config = readConfig();
   [X, Y, Z] = readPointCloud(config.plyPath, config.syntheticN);    
       
    %------------------------    
    X1 = X; Y1 = Y; Z1 = Z;    
    D = [X; Y; Z];
    M = [X1; Y1; Z1];
    n = size(D,2);
   
    % Translation values (a.u.):
    Tx = -1.25; Ty = -1.5;  Tz = 2.5;
    
    % Translation vector
    T = [Tx; Ty; Tz];
    
    % Rotation values (rad.):
    rx = -2.5;  ry = 2.8;  rz = -1.5;
    Rx = [1 0 0;0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)];
    Ry = [cos(ry) 0 sin(ry); 0 1 0; -sin(ry) 0 cos(ry)];
    Rz = [cos(rz) -sin(rz) 0; sin(rz) cos(rz) 0; 0 0 1];    
    R = Rx*Ry*Rz;
    
    % Transform data-matrix plus noise into model-matrix 
    M = R * M + repmat(T, 1, size(M,2));    
    
   
    idx=1:size(M,2);     
    Mperm = M(:,idx);
    M = Mperm;
       
    % Add noise to model and data
    D = D + 0.01*randn(3,size(D,2));
    
    % Generate Outliers
    nOutliers = round(n*config.OutlierRate/100.0);    
    
    idxRemove = 1:nOutliers;
    %idxRemove = randsample(n, nOutliers);
    D(:, idxRemove)=[];
      
    close all; plotPointClouds(M, D, 'b.', 'r.');   
    save(config.matPath, 'X', 'Y', 'Z', 'M', 'D');
    
    pcM = pointCloud(M'); pcD = pointCloud(D');
    pcwrite(pcM, config.plyOutA);
    pcwrite(pcD, config.plyOutB);
    
    %Write data for GoICP
    fout = fopen(config.txtOutA,'w');
    fprintf(fout, "%d\n", size(M,2));
    fprintf(fout, "%f %f %f\n", M);
    fclose(fout);
    
    fout = fopen(config.txtOutB,'w');
    fprintf(fout, "%d\n", size(D,2));
    fprintf(fout, "%f %f %f\n", D);
    fclose(fout);
    
    disp('------Data generated----');
    
    % Gen 4PCS script
    fpcs_script_path = [config.PCSFolder 'run.sh'];
    fp = fopen(fpcs_script_path, 'w');
    data_pref = '/home/intellhave/Dropbox/PointCloudReg/Rigid/';
    fprintf(fp, './Super4PCS -i %s %s -o 0.6 -d 0.01 -m -n 200 rs.txt -t 5\n', ...
                [data_pref config.plyOutA(4:end)], [data_pref config.plyOutB(4:end)]);
        
        
    fclose(fp);
    
    

end