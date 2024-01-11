%
% sim_dichot_ostats.m
%
% run a large N simulation to obtain the order statistics for the decision
% times of a large group split into biased and unbiased.
%
% parameterized by threshold, level of bias, and probability of an agent
% being biased.

clear;

N = 1e4;    % number of agents
h = 1;      % threshold
b = -h/2;    % strength of bias
a = 0.5;    % probability of agents being biased

dt = 1e-4;  % timestep

% initialize belief value initial conditions
nb = sum(rand(N,1)<a);
Xb = b*ones(nb,1); X0 = zeros(N-nb,1);  % make biased and unbiased agents
db = zeros(nb,1); d0 = zeros(N-nb,1);   % make vector to store decision
Tb = zeros(nb,1); T0 = zeros(N-nb,1);   % make vector to store times

% simulate biased agents
for k=1:nb, j=1;
    while abs(Xb(k,j))<h, Xb(k,j+1)=Xb(k,j)+dt+sqrt(2*dt)*randn; j=j+1; end
    Tb(k) = j; db(k)=sign(Xb(k,j));
end

% simulate unbiased agents
for k=1:N-nb, j=1;
    while abs(X0(k,j))<h, X0(k,j+1)=X0(k,j)+dt+sqrt(2*dt)*randn; j=j+1; end
    T0(k) = j; d0(k)=sign(X0(k,j));
end

dall = [db;d0]; Tall=[Tb;T0]; [Tmin,kmin]=min(Tall); [Tmax,kmax]=max(Tall);
[Tbmax,junk]=max(Tb); [T0max,junk]=max(T0);
Xb=[Xb,zeros(nb,Tmax-Tbmax)]; X0=[X0,zeros(N-nb,Tmax-T0max)];
Xall = [Xb;X0];

% plot everyone on a linear scale
 Nk = 250;    % granularity of plotting (larger is fewer)
% figure(1), hold on
% for k=1:Nk:N, plot(dt*[0:Tall(k)-1],Xall(k,1:Tall(k)),'linewidth',0.5,'Color',ones(1,3)/2); end
% xlabel('time','fontsize',30,'interpreter','latex');
% ylabel('$X_k$','fontsize',30,'interpreter','latex');
% set(gca,'fontsize',30);
% axis([1e-4 dt*Tmax*1.1 -1.1*h 1.1*h])
% 
% plot(dt*[0:Tall(kmin)-1],Xall(kmin,1:Tall(kmin)),'r','linewidth',3);
% plot(dt*[0:Tall(kmax)-1],Xall(kmax,1:Tall(kmax)),'b','linewidth',3);
% plot(dt*[0:1.1*Tall(kmax)-1],0*[0:1.1*Tall(kmax)-1]+h,'k','linewidth',3);
% plot(dt*[0:1.1*Tall(kmax)-1],0*[0:1.1*Tall(kmax)-1]-h,'k','linewidth',3);


% plot everyone on a log scale
figure(2), hold on
for k=1:Nk:N, plot(dt*[0:Tall(k)-1],Xall(k,1:Tall(k)),'linewidth',0.5,'Color',ones(1,3)/2); end
xlabel('time','fontsize',30);
ylabel('$X_k$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);
set(gca,'xscale','log');
set(gca,'xtick',10.^[round(log10(dt)):floor(log10(Tmax))]);
axis([1e-4 dt*Tmax*1.1 -1.1*h 1.1*h])

plot(dt*[0:Tall(kmin)-1],Xall(kmin,1:Tall(kmin)),'r','linewidth',3);
plot(dt*[0:Tall(kmax)-1],Xall(kmax,1:Tall(kmax)),'b','linewidth',3);
plot(dt*[0:1.1*Tall(kmax)-1],0*[0:1.1*Tall(kmax)-1]+h,'k','linewidth',3);
plot(dt*[0:1.1*Tall(kmax)-1],0*[0:1.1*Tall(kmax)-1]-h,'k','linewidth',3);

% now to bin and plot the distribution of decision times
% tbin = linspace(0,Tmax,50);
% ptall = hist(Tall,tbin)/N/dt;
% ptb = hist(Tb,tbin)/nb/dt;
% pt0 = hist(T0,tbin)/(N-nb)/dt;
% 
% figure(3), hold on, plot(dt*tbin,ptall,'k','linewidth',3);
% figure(3), hold on, plot(dt*tbin,ptb,'m','linewidth',3);
% figure(3), hold on, plot(dt*tbin,pt0,'c','linewidth',3);
% set(gca,'fontsize',30);
% xlabel('$T$','fontsize',30,'interpreter','latex');
% ylabel('$p(T)$','fontsize',30,'interpreter','latex');