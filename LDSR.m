function [AUC] = LDSR( X, outLabel, lambda1, alpha, beta)

nv = length(X);

max_mu = 1e10;
mu = 1e-4;
rho = 1.3;

x_dims = zeros(nv, 1);
ns = size(X{1},2);
Q = cell(nv, 1);
XtX = cell(nv, 1);
E = cell(nv,1);
Ze = cell(nv, 1);
sumXtX = 0;
for n =  1:nv
  x_dims(n) = size(X{n},1); 
end

for n = 1:nv
  for i = 1:size(X{n},2)
    if norm(X{n}(:,i)) >= 1e-8
        X{n}(:,i) = X{n}(:,i) ./ norm(X{n}(:,i));
    end
  end
  Ze{n} = zeros(ns, ns);
  XtX{n} = X{n}'*X{n};
  sumXtX = sumXtX + XtX{n};
  E{n} = sparse(x_dims(n), ns);
  Q{n} = zeros(x_dims(n), ns);
end

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
      R = diag(r);
      temp1 = mu*X{i}'*(X{i}*Zc+E{i}-X{i}-Q{i}/mu);
      temp2 = 2*alpha*R+mu*X{i}'*X{i};  
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
    
    if break_flag && max(max(abs(leq2))) < tol
       break;
    end
    mu = min(max_mu,mu*rho);   
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
