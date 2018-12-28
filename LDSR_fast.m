function [AUC] = LDSR_fast( X, label, outLabel, alpha, beta )

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
% Ze = cell(nv, 1);
sumXtX = 0;
% sumWtx = 0;
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
  
  W{n} = orth(X{n}');
  XW{n} = X{n} * W{n}(:, 1:md);
  Z{n} = zeros(ns, ns);
   
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

U = Zc * 0;

epi = 1e-3;
iter = 0;
maxIter = 1e6;
tol = 1e-6;

eta = 1e9;
tau = 10;
% rho0 = 15;
rho0 = 1;
rho = rho0;
gamma = 0;
lambda = 1;
    
while iter<maxIter
    iter = iter + 1;    
  
    sumZv = 0;
    for i = 1:nv
      temp1 = 2*alpha*eye(ns)+mu*XtX{i};
      temp2 = 2*alpha*Zc + mu*XtX{i} - mu*X{i}'*E{i} + mu*X{i}'*Q{i}/mu;
      Z{i} =  temp1 \ temp2;
      sumZv = sumZv + Z{i};
    end
    
%%

    [u, d, v] = lansvd(U, 1, 'L', struct('tol',1e-4)); % use U
%     [u, d, v] = svds(U, 1);
    if (d > 1)
        % This trick can handle the situation when rho is set too large.
        % The algorithm should converge when rho is larger than a unkown
        % constant. 
        derivative = (u * v');
        rho = min(rho0, sqrt(iter) / tau * d);
    else
        derivative = zeros(size(U));
    end
    % step 1: gradient descent
    A_new = Zc - eta /sqrt(iter) * (2*alpha*(nv*Zc-sumZv) + lambda * U);
    % step 2: shrinkage proximal ell_1 norm
    A_new = sign(A_new) .* max(abs(A_new) - eta /sqrt(iter) * gamma, 0);
    % step 3: update U
    U = U + tau /sqrt(iter)* (lambda * Zc - rho * derivative );
    Zc = A_new;
%%    
    
    %update J
%     temp = Zc + P/mu;
%     [U,sigma,V] = svd(temp,'econ');
%     sigma = diag(sigma);
%     svp = length(find(sigma>1/mu));
%     if svp>=1
%         sigma = sigma(1:svp)-1/mu;
%     else
%         svp = 1;
%         sigma = 0;
%     end
%     J = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
    %
    %update Zv

    %update Zc
%     Zc = 1/(mu+2*alpha*nv) * (-P+mu*J+2*alpha*sumZv);
    
    break_flag = 0;
    for i = 1:nv
      temp1 = X{i}-X{i}*Z{i};
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
%     leq2 = Zc - J;
%     P = P + leq2 * mu;

%      fprintf("leq1: %f \n", max(max(abs(leq1))));
%     disp('leq1: ');
%     disp();
%      fprintf("leq2: %f \n", max(max(abs(leq2))));
    
%     disp('leq2: ');
%     disp(max(max(abs(leq2))));
    
    if break_flag
%         && max(max(abs(leq2))) < tol
%        disp("break has been executed");
       break;
    end
    mu = min(max_mu,mu*rho);
end


% for i = 1:nv
%    Z{i} = Ublk .* Z{i};
% end

consistency = calc_consistency(Z, E, unq_dim, lambda1, lambda2, lambda3);
[~, ~, ~, AUC] = perfcurve(outLabel, consistency, 0);
disp(AUC);

end



function consistency = calc_consistency(Z, E, unq_dim, lambda1, lambda2, lambda3)

err1 = 0;
err2 = 0;
err3 = 0;
err4 = 0;

nv = length(Z);

for i = 1:nv
  for j = i+1:nv
    if i ~= j
%        s1 = sum(Ze{i}.*Ze{j}); 
       s1 = sum(Z{i}.*Z{j}); 
       err1 = err1 + abs(s1 / norm(s1));
       
       if unq_dim
          s2 = sum(E{i}.*E{j});
       else
          s2 = sum(E{i}) .* sum(E{j});
       end
       
       if norm(s2) == 0
          err2 = 0;
       else
          err2 = err2 + abs(s2 / norm(s2));
       end
      
%        err1 = 0;
%        err2 = 0;
       %s3 = sum(Ze{i}.*E{j}+Ze{j}.*E{i});
       %err3 = err3 + s3 / norm(s3);
       
       temp1 = Z{i} - Z{j};
       s3 = diag(temp1'*temp1)';
       err3 = err3 + abs(s3 / norm(s3));
%        err3 = 0;
       if unq_dim        
         temp2 = E{i} - E{j};
         s4 = diag(temp2'*temp2)';
       else
           
         temp2 = sum(E{i}) - sum(E{j});
         s4 = temp2;
%          s4 = diag(temp2'*temp2)';
       end
       
       if norm(s4) == 0
          err4 = 0;
%           err4+s3;
       else
          err4 = err4 + abs(s4 / norm(s4));
       end
       
       if norm(s4) == 0
         temp = 0; 
       else
         temp = err3 .* err4;
         temp = temp + temp / norm(temp);
       end
       
       err3 = temp;
       err4 = 0;
       
%         err3 = 0;
%        err1 = 0;
       
    end
  end 
end
%consistency = err2 - lambda1 * err1 - lambda2 * err3;

consistency = err1 - lambda1 * err2 - lambda2 * err3 - lambda3 * err4;

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
