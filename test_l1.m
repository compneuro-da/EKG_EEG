load TRUMP

y=data(:,1);

% set maximum regularization parameter lambda_max
lambda_max = l1tf_lambdamax(y);

%----------------------------------------------------------------------
% 	l1 trend filtering
%----------------------------------------------------------------------
[z1,status] = l1tf(y, 0.001*lambda_max);
[z2,status] = l1tf(y, 0.005*lambda_max);
[z3,status] = l1tf(y, 0.025*lambda_max);

% uncomment line below to solve l1 trend filtering problem using CVX
% [z1,status] = l1tf_cvx(y, 0.001*lambda_max);
% [z2,status] = l1tf_cvx(y, 0.005*lambda_max);
% [z3,status] = l1tf_cvx(y, 0.025*lambda_max);

%----------------------------------------------------------------------
% 	plot results
%----------------------------------------------------------------------
% xyzs = [x y z1 z2 z3];
% maxx = max(max(xyzs));
% minx = min(min(xyzs));

figure(1);
plot(y, 'k');
hold on; plot(z1,'r');
hold on; plot(z2,'b');
hold on; plot(z3,'g');

