function [P,Q] = prod_cum(D,alpha)
% for given demand events and alpha function will return
% cumulative production P and batch size Q
Q = [];
for j = 1:length(D)
    sigma = [];
    for k = j:length(D)
        big_pi = [];
        for l = j+1:k
            big_pi = [big_pi (1 - alpha(l))];
        end
        big_pi = prod(big_pi);
        sigma = [sigma D(k)*big_pi];
    end
    sigma = sum(sigma);
    Q = [Q alpha(j)*sigma];
end
P = cumsum(Q);


