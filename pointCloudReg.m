function [out] = pointCloudReg(M, D, config, method)

           
    if strcmp(method, 'SDRSAC')
        start_time = tic;
        out = SDRSAC(M, D, config);
        out.run_time = toc(start_time);
    
    % To be added - The case with correspondences
%     elseif strcmp(method, 'CSDRSAC')
%         
%         out = CSDRReg(M, D, [], config);        
                
    else 
        disp('=======WRONG METHOD==============');
        
    end
        
    
    
    



end