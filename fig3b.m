% Samantha Linn
% u1323207@gcloud.utah.edu
% 2023
% Code for Figure 3B

d = 1;                                   % Diffusivity
Mu = 1;                                  % Rightward drift
N = [5 10 100];                          % Number and range of agents
J = N;                                   % jth quickest agent
th = 1;                                  % Interval half length
WW = th.*linspace(-.99,.99,100);
q = ones(size(WW))./length(WW);          % Initial position probabilities
tmax = 1;                                % Cut off for short-time/long-time

kappa = 1000;                            % Number of terms in series approx
terms = 1000;                            % Number of terms in quadrature
tMax = 100;                              % Max time used in quadrature

colors = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.6350 0.0780 0.1840],...
    [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.3010 0.7450 0.9330],...
    [0.8500 0.3250 0.0980]};

prob = ones(length(WW),length(N)); 
bdy = zeros(length(WW),length(N));

D = d/(2*th)^2;
mu = Mu/(2*th);
W = (WW + th)./(2*th);
L = ceil(W(end)) - floor(W(1));
jj = J.*ones(size(N));

ts = logspace(-10,log10(tmax),terms)';
tl = logspace(log10(tmax),log10(tMax),terms)';
t = [ts; tl];

leg = cell(1,length(W));
fs0 = zeros(length(t),length(W));
fs1 = zeros(length(t),length(W));
cdf = zeros(length(t),length(W));

for i = 1:length(W)

    w = W(i)/L;
    s = (D.*ts)./L^2;
    time = 1:length(ts);
    for k = -floor((kappa-1)/2):ceil((kappa-1)/2)
        fs0(time,i) = fs0(time,i) + (w+2*k).*exp(-(w+2*k)^2./(4.*s));
        fs1(time,i) = fs1(time,i) + (1-w+2*k).*exp(-(1-w+2*k)^2./(4.*s));
    end
    fs0(time,i) = (D/(L^2))*fs0(time,i)./sqrt(4*pi.*s.^3);
    fs1(time,i) = (D/(L^2))*fs1(time,i)./sqrt(4*pi.*s.^3);

    s = (D.*tl)./L^2;
    time = length(ts)+1:length(t);
    for k = 1:kappa
        fs0(time,i) = fs0(time,i) + (D/(L^2))*2*pi*k*exp(-k^2*pi^2.*s).*sin(k*pi*w);
        fs1(time,i) = fs1(time,i) + (D/(L^2))*2*pi*k*exp(-k^2*pi^2.*s).*sin(k*pi*(1-w));
    end

    fs0(:,i) = fs0(:,i).*exp((-mu*W(i)/2/D)*ones(size(t)) - mu^2.*t./4./D).*q(i);
    fs1(:,i) = fs1(:,i).*exp((mu*(L-W(i))/2/D)*ones(size(t)) - mu^2.*t./4./D).*q(i);

    f = fs0 + fs1;
    for j = 2:length(t)
        cdf(j,i) = trapz(t(1:j),f(1:j,i));
    end
    leg{i} = strcat('$$x_0=\,\,$$',num2str(WW(i)/th),'$$\theta$$');

end
S = 1-sum(cdf,2);

for i = 1:length(W)
    for j = 1:length(N)
        prob(i,j) = jj(j)*nchoosek(N(j),jj(j))*...
            trapz(t,f(:,i).*(1-S).^(jj(j)-1).*(S.^(N(j)-jj(j))));
        bdy(i,j) = jj(j)*nchoosek(N(j),jj(j))*...
            trapz(t,fs1(:,i).*(1-S).^(jj(j)-1).*(S.^(N(j)-jj(j))));
    end
end

figure;
for i = 1:length(N)
    plot(WW,(bdy(:,i)./prob(:,i)),'Color',colors{i},'LineWidth',2);
    hold on
end
plot(WW,(1./(1+exp(-Mu*th/d))).*ones(size(WW)),'Color',colors{i+1},...
    'LineWidth',2,'LineStyle',':');
hold off
x = xlabel('$$y$$');
set(x,'Interpreter','latex')
xlim([-th th]);
y = ylabel('Prob(fastest agent decides w/ pos. DRIFT cond. on initial bias)');
y = ylabel('Prob($n(N)$ makes correct decision$\,| \,X_{n(N)}(0) = y$)');
set(y,'Interpreter','latex')
leg = {'$$N=5$$','$$N=10$$','$$N=100$$','$$\frac{1}{1+e^{-\mu\theta/D}}$$'};
l = legend(leg,'Location','southeast');
set(l,'Interpreter','latex')
xticks([-1 -.5 0 .5 1])
xticklabels({'-\theta','-\theta/2','0','\theta/2','\theta'});
xtickangle(0)
set(gca,'fontsize',15)
