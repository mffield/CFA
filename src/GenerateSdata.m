i = (1.5*pi/60:1.5*pi/30:1.5*pi)';
x = [cos(i); sin(i)];
y = [sin(i)+1; cos(i)-1]/2;
n = 20;
z = (0:1/(n-1):1)'*2-1;
xx = [kron(x,ones(n,1)) kron(y,ones(n,1)) kron(ones(length(x),1),z)]';
cc = kron(length(x):-1:1,ones(n,1))/length(x);
xx = xx([1 3 2],:);
% add some noise to the data, XX is the data we used for training
X = xx + .03*randn(3,n*length(x));
