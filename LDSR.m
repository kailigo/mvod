function [AUC] = LDSR( X, label, outLabel, ...
                              lambda1, alpha, beta)

nv = length(X);


max_mu = 1e10;
mu = 1e-4;
rho = 1.3;


n_c = zeros(1,length(unique(label)));
for i = 1:length(unique(label))
    indi = find(label==i);
    n_c(i) = length(indi);
end

Ublk = ones(n_c(1));
for i = 2:length(n_c)
    Ublk = blkdiag(Ublk, ones(n_c(i)));
end

x_dims = zeros(nv, 1);
ns = size(X{1},2);
Q = cell(nv, 1);
XtX = cell(nv, 1);
E = cell(nv,1);
Ze = cell(nv, 1);
sumXtX = 0;
sumWtx = 0;
W = cell(nv, 1);
for n =  1:nv
  x_dims(n) = size(X{n},1); 
end
md = min(x_dims);

unq_dim = 1;
if length(unique(x_dims)) ~= 1
    unq_dim = 0;
end

for n = 1:nv
  for i = 1:size(X{n},2)
    if norm(X{n}(:,i)) >= 1e-8
        X{n}(:,i) = X{n}(:,i) ./ norm(X{n}(:,i));
    end
  end
%   W{n} = orth(X{n}');
%   XW{n} = X{n} * W{n}(:, 1:md);
  
%   Ze{n} = zeros(ns, x_dims(n));
  Ze{n} = zeros(ns, ns);
   
%   wtX{n} = XW{n} * X{n};
%   sumWtx = sumWtx + wtX{n};
 
  XtX{n} = X{n}'*X{n};
  sumXtX = sumXtX + XtX{n};
  E{n} = sparse(x_dims(n), ns);
  Q{n} = zeros(x_dims(n), ns);
end
% 

% rng(0);

% P = rand(ns, ns);
% Zc = rand(ns, ns);

P = zeros(ns, ns);
Zc = zeros(ns, ns);


epi = 1e-3;
iter = 0;
maxIter = 1e6;
tol = 1e-6;
while iter<maxIter
    iter = iter + 1;    
    %update J
    temp = Zc + P/mu;
    [U,sigma,V] = svd(temp,'econ');
    sigma = diag(sigma);
    svp = length(find(sigma>1/mu));
    if svp>=1
        sigma = sigma(1:svp)-1/mu;
    else
        svp = 1;
        sigma = 0;
    end
    J = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
     
%     temp1 = zeros(size(Ze{1}));
%     temp1 = zeros(size(wtX{1}));
%     temp2 = temp1;
%     temp3 = temp1;
%     for i = 1:nv
%        temp1 = temp1 + wtX{i} * Ze{i};
%        temp2 = temp2 + XW{i}' * E{i};
%        temp3 = temp3 + XW{i}' * Q{i};
%     end
%     temp4 = -P+mu*J+mu*sumWtx-mu*temp1-mu*temp2+temp3;
%     temp5 = mu*eye(ns)+mu*sumWtX;
%     Zc = temp4 \ temp5;
  
    temp1 = zeros(size(Ze{1}));
    temp2 = temp1;
    temp3 = temp1;
    for i = 1:nv
       temp1 = temp1 + XtX{i} * Ze{i};
       temp2 = temp2 + X{i}' * E{i};
       temp3 = temp3 + X{i}' * Q{i};
    end
    temp4 = -P+mu*J+mu*sumXtX-mu*temp1-mu*temp2+temp3;
    temp5 = mu*eye(ns)+mu*sumXtX;
    Zc = temp5 \ temp4;
    
    break_flag = 0;
    for i = 1:nv

      r = 0.5 ./ sqrt(diag(Ze{i}*Ze{i}')+epi);
%        r = 0.5 ./ sqrt(diag(Ze{i}'*Ze{i})+epi);    
      R = diag(r);
      temp1 = mu*X{i}'*(X{i}*Zc+E{i}-X{i}-Q{i}/mu);
      temp2 = 2*alpha*R+mu*X{i}'*X{i};  
%       temp2 = alpha*eye(n)+mu*X{i}'*X{i};
      Ze{i} = - temp2 \ temp1;
      
      %update Ev
      temp1 = X{i}-X{i}*Zc-X{i}*Ze{i};
      temp2 = temp1+Q{i}/mu;
      E{i} = solve_l1l2(temp2, beta/mu);
      
      % update Qv
      leq1 = temp1 - E{i};
      Q{i} = Q{i} + leq1 * mu; 
      
      if max(max(abs(leq1))) < tol 
          break_flag = 1;
          break;
      end
    end
    %update P
    leq2 = Zc - J;
    P = P + leq2 * mu;

%      fprintf("leq1: %f \n", max(max(abs(leq1))));
%     disp('leq1: ');
%     disp();
%      fprintf("leq2: %f \n", max(max(abs(leq2))));
    
%     disp('leq2: ');
%     disp(max(max(abs(leq2))));
   err1 = norm(leq1, 'fro') / norm(X{i}, 'fro');
    err2 = norm(leq2, 'fro') / norm(X{i}, 'fro');
    err(iter) = max([err1, err2]);
    
    if break_flag && max(max(abs(leq2))) < tol
%        disp("break has been executed");
       break;
    end
    mu = min(max_mu,mu*rho);
    
    for i = 1:nv
        Z{i} = Ublk .* (Zc + Ze{i});
%         Z{i} = Ublk .* (Zc + Ze{i});
%        Z{i} = Zc + Ze{i};
    end
    
    err2_flag = 0;
    if any(any(E{i}))
       err2_flag = 1;
    end
        
    for i = 1:ns
      err1 = 0;
      err2 = 0;
      for j = 1:nv
        err1 = err1 + norm(Ze{j}(:, i));
        if any(E{j}(:, i))
           err2 = err2 + norm(E{j}(:, i));
        else
           err2 = 0;  
        end
      end
      consistency(i) = -err1 - lambda1 * err2;  
    end
   
    if length(consistency) > 1
        [~, ~, ~, AUC] = perfcurve(outLabel, consistency, 0);
    end
    

    AUCS(iter) = AUC;

end
disp(AUC);

end

function [E] = solve_l1l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end
end

function [x] = solve_l2(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end
end
