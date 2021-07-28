function xY = pca_clean(ROI,TR,cut_off)


  
    [y  X0]  = filt_tc(ROI,TR,cut_off);

    %-Confounds
    %--------------------------------------------------------------------------
    xY.X0     = X0;

    %-Compute regional response in terms of first eigenvariate
    %--------------------------------------------------------------------------
    [m,n]   = size(y);
    if m > n
        [v,s,v] = svd(y'*y);
        s       = diag(s);
        v       = v(:,1);
        u       = y*v/sqrt(s(1));
    else
        [u,s,u] = svd(y*y');
        s       = diag(s);
        u       = u(:,1);
        v       = y'*u/sqrt(s(1));
    end
    d       = sign(sum(v));
    u       = u*d;
    v       = v*d;
    Y       = u*sqrt(s(1)/n);
    
    %-Set in structure
    %--------------------------------------------------------------------------
    xY.y    = y;
    xY.u    = Y;
    xY.v    = v;
    xY.s    = s;
    xY.m    = mean(y,2);
    xY.mp   = mean(ROI,2);
    