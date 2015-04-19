function [K_cost] = NPV_K(q,rho,Q,D,t,alpha,K,beta)
% for given input: prod rate, interest, demand events, batch size, time,
% alpha, set-up cost and beta, function will return NPV equation inventory cost
K_cost = K*alpha.*exp(-rho*(t-D/q+(beta/q)*Q));
K_cost = sum(K_cost);