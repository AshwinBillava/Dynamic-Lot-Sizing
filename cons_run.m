function z = cons_run(x)
% for input: 1     0     0     1     1     0     1     0     1     0
% returns:   1     0     0     2     1     0     1     0     1     0
m = size(x,1);
y = [reshape([zeros(1,m);x.'],[],1);0];
z = y;
p = find(~y);
d = 1-diff(p);
y(p) = [0;d];
y = reshape(cumsum(y(1:end-1)),[],m).';
y(:,1) = [];
z(p) = [d;0];
z = reshape(cumsum(-z(1:end-1)),[],m).';
z(:,end) = [];