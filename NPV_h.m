function [h] = NPV_h(c,q,rho,Q,D,t)
% for given input: unit cost, prod rate, interest, demand events, batch size
% and time function will return
% NPV equation inventory cost
h = (c*q/rho)*(ones(1,length(Q))-exp(-(rho/q)*Q)).*exp(-rho*(t-D/q));
h = sum(h);