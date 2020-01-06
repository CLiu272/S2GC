function [w, funcVal, S] = S2GC(Data, lambda, rho, tau, ViewIdx, w0)

dimension = size(Data.X, 2);
numSample = size(Data.X, 1);
funcVal = [];

if nargin <6
    w0 = zeros(dimension, 1);
end

wz= w0;
wz_old = w0;

t = 1;
t_old = 0;
iter = 0;
gamma = 1;
gamma_inc = 2;

maxIter = 100;

while iter < maxIter
    
    alpha = (t_old - 1) /t;   
    ws = (1 + alpha) * wz - alpha * wz_old;
    
    % compute function value and gradients of the search point 
    [gws, Fs] = gradVal_eval(ws,rho,ViewIdx);
    k = 10;
    [G,S] = CalGraphWeight(Data.X,ws,k);
    gws = gws + tau * (Data.X)'* G * (Data.X) * ws;
       
    innerIter = 0;
    maxInnerIter = 1000;
    while innerIter < maxInnerIter
        
        wzp = l1_projection(ws - gws/gamma, lambda / gamma);
        Fzp = funVal_eval (wzp);
        
        delta_wzp = wzp - ws;
        r_sum = norm(delta_wzp, 'fro')^2;

        Fzp_gamma = Fs + sum(sum(delta_wzp.* gws))...
            + gamma/2 * norm(delta_wzp, 'fro')^2;
        
        if (r_sum <=1e-20)
            break;
        end
        
        if (Fzp <= Fzp_gamma)
            break;
        else
            gamma = gamma * gamma_inc;
        end
        innerIter = innerIter + 1;
    end
    
    wz_old = wz;
    wz = wzp;
    
    funcVal = cat(1, funcVal, Fzp);
    % test stop condition.
    
    if iter >= 50
        if (abs( funcVal(end) - funcVal(end-1) ) <= 10e-6)
            break;
        end
    end
    
    iter = iter + 1;
    fprintf('###### Iteration ######## %d \n\n ',iter);
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
end

w = wzp;

%% subfunction
    function [z , l1_comp_val] = l1_projection (v, beta)     
        z = zeros(size(v));
        vp = v - beta/2;
        z (v> beta/2)  = vp(v> beta/2);
        vn = v + beta/2;
        z (v< -beta/2) = vn(v< -beta/2);
        l1_comp_val = sum(sum(abs(z)));
    end
    
    function [grad_w, funcVal] = gradVal_eval(w,rho,ViewIdx)
          
        grad_w = gradandorg_neglogparlike(w,Data);
        
        grad_cs = zeros(dimension,1);
        
        for i = 1:length(ViewIdx)
            for j = 1:length(ViewIdx)
                grad_cs(ViewIdx{i}) = grad_cs(ViewIdx{i}) + Data.X(:,ViewIdx{i})'* ( Data.X(:,ViewIdx{i}) * w(ViewIdx{i}) - Data.X(:,ViewIdx{j}) * w(ViewIdx{j}) ); 
            end
        end
        
        grad_w = grad_w + rho * grad_cs;
        
        funcVal = 0;
        funcVal = funcVal + neglogparlike(w,Data);
        
    end

    function [funcVal] = funVal_eval (w)
       
        funcVal = 0;
        funcVal = funcVal + neglogparlike(w,Data);

    end


end


%% individual function

function [L]=neglogparlike(w,Data)
    % Compute log likelihood L    
    X=Data.X;
    freq=Data.freq;
    cens=Data.cens;
    atrisk=Data.atrisk;
    obsfreq = freq .* ~cens;
    Xb = X*w;
    r = exp(Xb);
    risksum = flipud(cumsum(flipud(freq.*r)));
    risksum = risksum(atrisk);
    L = obsfreq'*(Xb - log(risksum));
    L = -L;
    
end

function [dl]=gradandorg_neglogparlike(w,Data)

    % Compute log likelihood L
    X=Data.X;
    freq=Data.freq;
    cens=Data.cens;
    atrisk=Data.atrisk;

    obsfreq = freq .* ~cens;
    Xb = X*w;
    r = exp(Xb);
    risksum = flipud(cumsum(flipud(freq.*r)));
    risksum = risksum(atrisk);

    [~,p] = size(X);
    Xr = X .* repmat(r.*freq,1,p);
    Xrsum = flipud(cumsum(flipud(Xr)));
    Xrsum = Xrsum(atrisk,:);
    A = Xrsum ./ repmat(risksum,1,p);
    dl = obsfreq' * (X-A);
    dl = -dl';
    
end