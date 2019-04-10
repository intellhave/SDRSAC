function [maps, mapeq, maps_b, mapeq_b] = genSDPConstraintMap(n,k, W)

    maps =  {}; %cell(2*n,1);     
    mapeq = {};
    eqCount = 1; ieqCount = 1;
    for i=1:n
        A = zeros(n,n);
        A(i,:) = 1;
        Av = [A(:); 0];
        
        M = sparse(n*n+1, n*n+1);
        M(n*n+1, :) = Av;
        maps{ieqCount} = M; ieqCount = ieqCount + 1;
    end
    
    for i=1:n
        A = zeros(n,n);
        A(:,i) = 1;
        Av = [A(:); 0];
        
        M = sparse(n*n+1, n*n+1);
        M(n*n+1, :) = Av;
        maps{ieqCount} = M; ieqCount = ieqCount + 1;
    end
    maps_b = ones (2*n, 1);
    
    
    %----------------------------------------------------
    
    % 1^TX1 = k
    M = sparse(n*n+1, n*n+1);
    M(n*n+1,1:n*n) = 1;
    mapeq{eqCount} = M; eqCount = eqCount+1;
    mapeq_b = k;
        
    % Sum Yqrst = k
%     M = ones(n*n, n*n);
%     M = [M; zeros(1,n*n)];
%     M = [M zeros(n*n+1,1)];
%     mapeq{eqCount} = M; eqCount = eqCount + 1;
%     mapeq_b = [mapeq_b; k*k];       

    % Wherever W(i,j) < 0, force Y(i,j) to 0
%     for i=1:size(W,1)
%         for j=1:size(W,2)
%            if (W(i,j) < 0)
%                M = sparse(n*n+1, n*n+1); M(i,j) = 1;
%                mapeq{eqCount} = M; eqCount = eqCount + 1;
%                mapeq_b = [mapeq_b; 0];    
%             end
%         end
%     end
% %     
%     
    % Force Yqrst = 0;    
%      for p=1:n
%        for q=1:n
%            for s=1:n
%                for t=1:n
%                    i = p + (q-1)*n;
%                    j = s + (t-1)*n;
%                    if (p==s && q~=t)||(p~=s && q==t)                    
%                        M = sparse(n*n+1, n*n+1); M(i,j) = 1;
%                        mapeq{eqCount} = M; eqCount = eqCount + 1;
%                        mapeq_b = [mapeq_b; 0];                                                           
%                    else
% %                        M = sparse(n*n+1, n*n+1);
% %                        M(i,j) = 1; M(i,n*n+1) = -1;
% %                        maps_b = [maps_b;0];
% %                        maps{ieqCount} = M; ieqCount = ieqCount + 1;
% %                        
% %                        if (i~=j)
% %                            M = sparse(n*n+1, n*n+1);
% %                            M(i,j) = 1; M(j,n*n+1) = -1;
% %                            maps_b = [maps_b;0];
% %                            maps{ieqCount} = M; ieqCount = ieqCount + 1;
% %                        end
%                    end
%                    
%                        
%                    
%                 end
%             end
%         end
%     end
    
    
    
end