function [m,i] = prr(time,dem_c,prodrate)
% for given time and cum demand this function returns
% demand time when minimum start time is required
% demand time returned as i
% minimum start time returned as m
t_prime = time - dem_c/prodrate;
[m,ii] = min(t_prime);

% if two consecutive minimum then use the last minimum
[ma,mb] = ismember(t_prime,m);
ii = find(mb,1,'last');
i = time(ii);