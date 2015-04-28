close all
clear all
clc
%%% generate demand
% n = 10;
% t = 1:n;
t = [3 4 6 8 9 10 14 15 19 20]; %example 2
% t = [1 3 6 10 15]; %example 1

% lambda = 2;
% demand = poissrnd(lambda,1,n);
demand = [8 6 8 4 6 7 8 5 9 7]; %example 2
% demand = [1 1 1 1 1]; %example 1
demand_cum = cumsum(demand);

%%% plot cumulative demand
stairs([0 t t(end)+1],[0 demand_cum demand_cum(end)])
xlabel('Time');
ylabel('Demand');
hold on

%%% set production rate
% q = sum(demand_cum(end))/n
q = 5; % example 2
% q = 1; % example 1

%%% minimum q for All-At-Once
%[q_min,q_min_ind] = nanmin((ones(1,n)*demand_cum(end)-demand_cum)./(ones(1,n)*n-t))
%plot([(n-(demand_cum(end)/q_min)),n],[0,(n-(n-(demand_cum(end)/q_min)))*q_min],'r')

%%% elemintate zero demand events
% disp('Time & Cum Demand after zero demand events')
demand = [demand_cum(1) diff(demand_cum)]; % re-create demand if not specified
t_ind = find(demand>0);
t = t(t_ind); % time array after zero demand event eliminated
demand = demand(t_ind); % demand array after zero demand event eliminated
demand_cum = demand_cum(t_ind); % cumulative demand after zero demand event eliminated

%%% find production restriction rate (PRR)
ind = 0; %the time of the current dominated event
kk = 0;
prod_end_t = [];
prod_start_t0 = [];
while ind < t(end)
    time = t(kk+1:end);
    dem_c = demand_cum(kk+1:end);
    %%% find minimum time for time & dem_c
    [min_start,ind] = prr(time,dem_c,q); %minimum time t' and the time of corresponding event
    kk = kk + find(time==ind); %what indix belongs to the time ind
    prod_end_t = [prod_end_t ind]; %prod end time
    prod_start_t0 = [prod_start_t0 min_start]; %prod start time at time zero
end
[lia,locb] = ismember(prod_end_t,t);
D_bar_j = demand_cum(locb); %cumulative demand after PRR
prod_end = D_bar_j;
prod_start = [0 D_bar_j(1:end-1)];
D_j = [prod_end(1) diff(prod_end)];
prod_start_t = prod_end_t - D_j./q;
t_j = prod_end_t;

%%% plot production rate curves after PRR elimination
for j = 1:length(prod_start_t)
    %plot([dem_time(j),t(end)],[0,(t(end)-dem_time(j))*q],'g')
    %plot([prod_start_t0(j),prod_end_t(j)],[0,(prod_end_t(j)-prod_start_t0(j))*q],'g')
    plot([prod_start_t(j),prod_end_t(j)],[prod_start(j),prod_end(j)],'g')
end

%%% plot updated cumulative demand curves
t_c = prod_end_t;
stairs([0 t_c t_c(end)+1],[0 D_bar_j D_bar_j(end)],'k-.')

%%% Specify Inputs
K = 36; % Setup-Cost
c = 10; % Unit Production Cost
rho = 0.1; % Interest Rate
h = c*rho; % Inventory Holding Cost
beta = 1; % begin setup cost = 0; end setup cost = 1

%%% generate alpha combinations
n = length(prod_start);
D = [0:2^n-1]';
B = rem(floor(D*pow2(-(n-2):0)),2);
alpha_perm = [ones(length(B)/2,1) B(1:length(B)/2,:)];
% alpha = alpha_perm(2,:)

%%% compute batch size Q for alpha
% [P_bar_j,Q_j] = prod_cum(D_j,alpha);

%%% compute x and y for TC equation
y = [];
x = [];
Q_perm = [];
for jj = 1:length(alpha_perm)
    alpha = alpha_perm(jj,:);
    [P_bar_j,Q_j] = prod_cum(D_j,alpha);
    Q_perm(jj,:) = Q_j;
    y(jj) = TC_y(D_j,Q_j,q,t_j);
    x(jj) = TC_x(alpha);
end
% adjust for PRR eliminated demand
adj_factor = cons_run(not(ismember(t,prod_end_t)));
adj = adj_factor.*demand.*[diff(t) 0];
adj = sum(adj);
y = y - adj;
y = y';
x = x';
% plot the scatter plot for inventory cost vs # setup
figure()
scatter(x,y)
xlabel('# of Setups');
ylabel('Time-weighted Inventory');
% compute TC cost = Kx + hy
TC = K*x + h*y;

%%% compute NPV cost
K_cost = [];
h_cost = [];
for kk = 1:length(alpha_perm)
    alpha = alpha_perm(kk,:);
    [P_bar_j,Q_j] = prod_cum(D_j,alpha);
    h_cost(kk) = NPV_h(c,q,rho,Q_j,D_j,t_j);
    K_cost(kk) = NPV_K(q,rho,Q_j,D_j,t_j,alpha,K,beta);
end
% discounted cost value of requirements adjusted for PRR eliminated demand
dcvr = c*demand.*exp(-rho*t);
dcvr = sum(dcvr);
h_cost = h_cost-dcvr;

NPV = h_cost + K_cost;

%%% find alpha minimizing NPV cost
[min_npv,min_ind] = min(NPV);
result = ['Alpha that minimizes NPV cost: ',num2str(alpha_perm(min_ind,:))];
disp(result)





