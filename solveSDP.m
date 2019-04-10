% Solve the SDP relaxation problem
% max trace(WY)
% s.t.  Y >= xx'

function [X] = solveSDP(W, k)
    
    initSDP()
    nsquare = size(W,1);
    n = sqrt(nsquare);
    W = [W  zeros(nsquare, 1)];
    W = [W; zeros(1, nsquare+1)];
    
    [maps, mapeq, maps_b, mapeq_b] = genSDPConstraintMap(n, k ,W);
    
    
    model = ccp_model('SDP');
        Y = var_sdp(nsquare+1, nsquare+1);                 
        model.add_variable(Y);
        model.maximize(inprod(W,Y));                       
        %model.minimize(inprod(W',Y));                       
        model.add_affine_constraint(Y(nsquare+1, nsquare+1)==1)
        model.add_affine_constraint(Y(nsquare+1,:) == Y(:, nsquare+1)');
        model.add_affine_constraint(maps*Y <= maps_b);                      
        model.add_affine_constraint(mapeq*Y == mapeq_b);                           
        model.add_affine_constraint(trace(Y) == k+1);    
        
        % Early stopping
        model.setparameter('tol', 0.1);
        model.setparameter('tolADM', 0.1);
        
        model.setparameter('printlevel',0);
       
    model.solve()
    
    
    Y = get_value(Y);
    t = Y(end, 1:end-1);
    X1 = reshape(t, n, n);

    %disp('----------LINER PROJECTION------------------');
    X = linearProjection(X1(:));
    
    % Extract k elements
    idx = find(X==1);
    score = X1(idx);
    [~,sidx]= sort(score);
    idx(sidx(n-k+1:end))=[];
    X(idx) = 0;
      

end