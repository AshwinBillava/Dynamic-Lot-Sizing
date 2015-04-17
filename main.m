close all;
clear all;
clc;

% generate demand
% n = 10;
% t = 1:n;
t = [3 4 6 8 9 10 14 15 19 20];
% lambda = 2;
% demand = poissrnd(lambda,1,n);
demand = [8 6 8 4 6 7 8 5 9 7];
demand_cum = cumsum(demand);

% plot cummulative demand
stairs([0 t t(end)+1],[0 demand_cum demand_cum(end)])
hold on

% set production rate
% q = sum(demand_cum(end))/n
q = 5;

% minimum q for All-At-Once
%[q_min,q_min_ind] = nanmin((ones(1,n)*demand_cum(end)-demand_cum)./(ones(1,n)*n-t))
%plot([(n-(demand_cum(end)/q_min)),n],[0,(n-(n-(demand_cum(end)/q_min)))*q_min],'r')

% elemintate zero demand events
disp('Time & Cum Demand after zero demand events')
demand = [demand_cum(1) diff(demand_cum)];
t_ind = find(demand>0);
t = t(t_ind);
demand_cum = demand_cum(t_ind);

t_prime = t - demand_cum/q;
% find production restriction rate
ind = 0; %the time of the current dominated event
kk = 0;
dem_ind = [];
dem_time = [];
while ind < t(end)
    disp('find min of following')
    time = t(kk+1:end);
    dem_c = demand_cum(kk+1:end);
    [min,ind] = prr(time,dem_c,q); %minimum time t' and the time of corresponding event
    ind %not a valid statement
    kk = kk + find(time==ind); %what indix belongs to the time ind
    dem_ind = [dem_ind ind];
    dem_time = [dem_time min];
end



%plot production rate curves after PRR elemination
for j = 1:length(dem_time)
    %plot([dem_time(j),t(end)],[0,(t(end)-dem_time(j))*q],'g')
    plot([dem_time(j),dem_ind(j)],[0,(dem_ind(j)-dem_time(j))*q],'g')
end

%plot updated cumulative demand curves
[lia,locb] = ismember(dem_ind,t);
demand_cum_c = demand_cum(locb);

t_c = dem_ind;
stairs([0 t_c t_c(end)+1],[0 demand_cum_c demand_cum_c(end)],'k-.')