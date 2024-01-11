%
% sim_tri_dichot_ostats.m
%
% run a large N simulation to obtain the order statistics for the decision
% times of a large group split into biased and unbiased.
%
% here we consider a group making a decision from 3 alternatives. beliefs
% described by LLRs thus evolve on a planar domain with triangular boundary
%

N = 2E4;    % number of agents
h = 1/2/sqrt(3);      % threshold parameter (apothem of eq triangle)
b = 1/2;    % strength of bias (fraction of the apothem toward choice 3)
a = 0.5;    % probability of agents being biased
i = sqrt(-1);   % imaginary unit

dt = 1e-4;  % timestep

% s = 2*sqrt(3)*h; % side length of triangle
l = cos(pi/6); r=sin(pi/6); % long and short side of 30/60/90 tri

% drift will just be up (alt 1 is correct); now construct covariance matrix

% initialize belief value initial conditions (do in complex plane)
nb = sum(rand(N,1)<a);
Xb = h*b*(l-r*i)*ones(nb,1); X0 = zeros(N-nb,1);  % make biased and unbiased agents
db = zeros(nb,1); d0 = zeros(N-nb,1);   % make vector to store decision
Tb = zeros(nb,1); T0 = zeros(N-nb,1);   % make vector to store times

% simulate biased agents
for k=1:nb, j=1;
    x = real(Xb(k,j)); y = imag(Xb(k,j));
    while y<h && y>-2*h*(3*x+1) && y>2*h*(3*x-1)
        c = 0;%randn/sqrt(2);
        Xb(k,j+1)=Xb(k,j)+i*dt+sqrt(2*dt)*(randn+c+(randn+c)*i)/4;
        j=j+1;
        x = real(Xb(k,j)); y = imag(Xb(k,j));
    end
    Tb(k) = j;
    if y<2*h*(3*x-1), db(k)=3; end
    if y<-2*h*(3*x+1), db(k)=2; end
    if y>h, db(k)=1; end
end

% simulate unbiased agents
for k=1:N-nb, j=1;
    x = real(X0(k,j)); y = imag(X0(k,j));
    while y<h && y>-2*h*(3*x+1) && y>2*h*(3*x-1)
        c = 0;%randn/sqrt(2);
        X0(k,j+1)=X0(k,j)+i*dt+sqrt(2*dt)*(randn+c+(randn+c)*i)/4;
        j=j+1;
        x = real(X0(k,j)); y = imag(X0(k,j));
    end
    T0(k) = j;
    if y<2*h*(3*x-1), d0(k)=3; end
    if y<-2*h*(3*x+1), d0(k)=2; end
    if y>h, d0(k)=1; end
end

dall = [db;d0]; Tall=[Tb;T0]; [Tmin,kmin]=min(Tall); [Tmax,kmax]=max(Tall);
[Tbmax,junk]=max(Tb); [T0max,junk]=max(T0);
Xb=[Xb,zeros(nb,Tmax-Tbmax)]; X0=[X0,zeros(N-nb,Tmax-T0max)];
Xall = [Xb;X0];

% plot everyone on a linear scale
Nk = 1000;    % granularity of plotting (larger is fewer)
figure(1), hold on
for k=1:Nk:N, plot(real(Xall(k,1:Tall(k))),imag(Xall(k,1:Tall(k))),'linewidth',0.5,'Color',ones(1,3)/2); end
xlabel('$X_k^1$','fontsize',30,'interpreter','latex');
ylabel('$X_k^2$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30);
axis([-0.55 0.55 -0.65 0.35])

plot(real(Xall(kmax,1:Tall(kmax))),imag(Xall(kmax,1:Tall(kmax))),'c','linewidth',3);
plot(real(Xall(kmin,1:Tall(kmin))),imag(Xall(kmin,1:Tall(kmin))),'m','linewidth',3);

plot([-0.5, 0.5],[h, h],'g','linewidth',5);
plot([-0.5, 0],[h, -2*h],'r','linewidth',5);
plot([0, 0.5],[-2*h, h],'r','linewidth',5);
plot(real(Xall(kmin,Tall(kmin))),imag(Xall(kmin,Tall(kmin))),'m.','markersize',50);
plot(real(Xall(kmax,Tall(kmax))),imag(Xall(kmax,Tall(kmax))),'c.','markersize',50);
plot(h*b*l,-h*b*r,'ko','linewidth',3,'markersize',10);
plot(0,0,'ko','linewidth',3,'markersize',10);