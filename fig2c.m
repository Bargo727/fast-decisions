% Samantha Linn
% u1323207@gcloud.utah.edu
% 2023
% Code for Figure 2C

D = 1;                                   % Diffusivity
Mu = linspace(-8,8,500);                 % Rightward drift
n = [10^4 10^5 10^6];                    % Number and range of agents
th = 1;                                  % Interval half length
W = th/4;                                % Initial positions of agents
tmax = 1;                                % Cut off for short-time/long-time

kappa = 1000;                            % Number of terms in series approx
terms = 1000;                            % Number of terms in quadrature
tMax = 100;                              % Max time used in quadrature

colors = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.6350 0.0780 0.1840],...
    [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.3010 0.7450 0.9330],...
    [0.8500 0.3250 0.0980]};

prob = ones(length(n),length(Mu));  
bdy = zeros(length(n),length(Mu));

D = D/(2*th)^2;
mu = Mu./(2*th);
w = (W + th)./(2*th);
L = 1;

ts = logspace(-10,log10(tmax),terms)';
tl = logspace(log10(tmax),log10(tMax),terms)';
t = [ts; tl];

fs0 = zeros(length(t),1);
fs1 = zeros(length(t),1);
cdf = zeros(length(t),length(mu));
S = cdf;

s = (D.*ts)./L^2;
time = 1:length(ts);
for k = -floor((kappa-1)/2):ceil((kappa-1)/2)
    fs0(time) = fs0(time) + (w+2*k).*exp(-(w+2*k)^2./(4.*s));
    fs1(time) = fs1(time) + (1-w+2*k).*exp(-(1-w+2*k)^2./(4.*s));
end
fs0(time) = (D/(L^2))*fs0(time)./sqrt(4*pi.*s.^3);
fs1(time) = (D/(L^2))*fs1(time)./sqrt(4*pi.*s.^3);

s = (D.*tl)./L^2;
time = length(ts)+1:length(t);
for k = 1:kappa
    fs0(time) = fs0(time) + (D/(L^2))*2*pi*k*exp(-k^2*pi^2.*s).*sin(k*pi*w);
    fs1(time) = fs1(time) + (D/(L^2))*2*pi*k*exp(-k^2*pi^2.*s).*sin(k*pi*(1-w));
end

fs0 = fs0.*exp((-w/2/D)*ones(size(t))*mu - (t./4./D)*(mu.^2));
fs1 = fs1.*exp(((L-w)/2/D)*ones(size(t))*mu - (t./4./D)*(mu.^2));

f = fs0 + fs1;
for i = 1:size(f,2)
    for j = 2:length(t)
        cdf(j,i) = trapz(t(1:j),f(1:j,i));
    end
    S(:,i) = 1 - cdf(:,i);
end

for i = 1:length(n)
    N = n(i);
    for j = 1:length(mu)
        prob(i,j) = N*trapz(t,f(:,j).*(S(:,j).^(N-1)));
        bdy(i,j) = N*trapz(t,fs1(:,j).*(S(:,j).^(N-1)));
    end
end

figure;
for i = 1:length(n)
    plot(Mu,bdy(i,:)./prob(i,:),'Color',colors{i},'LineWidth',3);
    hold on
end
hold off
x = xlabel('$$\mu$$');
set(x,'Interpreter','latex')
xlim([min(Mu) max(Mu)]);
y = ylabel(strcat('Prob','$$(X_{n(1)}(T_1) = \theta > 0 \,|\, X_{n(1)}(0) > 0)$$'));
set(y,'Interpreter','latex')
leg = {'$$N=10^4$$','$$N=10^5$$','$$N=10^6$$'};
l = legend(leg,'Location','northeast');
set(l,'Interpreter','latex')
set(gca,'fontsize',18)

figure;
for i = 1:length(n)
    semilogy(Mu,1-(bdy(i,:)./prob(i,:)),'Color',colors{i},'LineWidth',3);
    hold on
end
hold off
x = xlabel('$$\mu$$');
set(x,'Interpreter','latex')
xlim([min(Mu) max(Mu)]);
y = ylabel(strcat('Prob','$$(X_{n(1)}(T_1) = -\theta < 0 \,|\, X_{n(1)}(0) > 0)$$'));
set(y,'Interpreter','latex')
set(gca,'fontsize',18)
