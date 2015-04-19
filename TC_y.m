function [y] = TC_y(D,Q,q,t)
% for given demand events, batch size, prod rate and time function will return
% TC equation time-weighted inventory, y
y = Q.*D/q - Q.*t - (Q.^2)/(2*q) + D.*t;
y = sum(y);