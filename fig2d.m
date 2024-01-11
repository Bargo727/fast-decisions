% Samantha Linn
% u1323207@gcloud.utah.edu
% 2023
% Code for Figure 2D

d = [1 1/2 1/4 1/8];                     % Diffusivity (largest value first)
Mu = 1;                                  % Rightward drift
N = round(logspace(1,5,500));            % Number and range of agents
Th = 1;                                  % Interval half length
WW = Th./2;                              % Initial positions of agents
q = ones(size(d))./length(d);            % Initial position probabilities
tmax = 1;                                % Cut off for short-time/long-time

kappa = 1000;                            % Number of terms in series approx
terms = 1000;                            % Number of terms in quadrature
tMax = 100;                              % Max time used in quadrature

colors = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.6350 0.0780 0.1840],...
    [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.3010 0.7450 0.9330],...
    [0.8500 0.3250 0.0980]};

th = Th;
D = d./(2*Th)^2;
mu = Mu/(2*Th);
W = (WW + Th)./(2*Th);
L = ceil(W(end)) - floor(W(1));
w = W/L;

prob = ones(length(D),length(N));
prob1 = prob;
ts = logspace(-10,log10(tmax),terms)';
tl = logspace(log10(tmax),log10(tMax),terms)';
t = [ts; tl];

fs0 = zeros(length(t),length(D));
fs1 = zeros(length(t),length(D));
cdf = zeros(length(t),length(D));

for i = 1:length(D)

    dd = D(i);
    s = (dd.*ts)./L^2;
    time = 1:length(ts);
    for k = -floor((kappa-1)/2):ceil((kappa-1)/2)
        fs0(time,i) = fs0(time,i) + (w+2*k).*exp(-(w+2*k)^2./(4.*s));
        fs1(time,i) = fs1(time,i) + (1-w+2*k).*exp(-(1-w+2*k)^2./(4.*s));
    end
    fs0(time,i) = (dd/(L^2))*fs0(time,i)./sqrt(4*pi.*s.^3);
    fs1(time,i) = (dd/(L^2))*fs1(time,i)./sqrt(4*pi.*s.^3);

    s = (dd.*tl)./L^2;
    time = length(ts)+1:length(t);
    for k = 1:kappa
        fs0(time,i) = fs0(time,i) + (dd/(L^2))*2*pi*k*exp(-k^2*pi^2.*s).*sin(k*pi*w);
        fs1(time,i) = fs1(time,i) + (dd/(L^2))*2*pi*k*exp(-k^2*pi^2.*s).*sin(k*pi*(1-w));
    end

    fs0(:,i) = fs0(:,i).*exp((-mu*W/2/dd)*ones(size(t)) - mu^2.*t./4./dd).*q(i);
    fs1(:,i) = fs1(:,i).*exp((mu*(L-W)/2/dd)*ones(size(t)) - mu^2.*t./4./dd).*q(i);

    f = fs0 + fs1;
    for j = 2:length(t)
        cdf(j,i) = trapz(t(1:j),f(1:j,i));
    end

end
S = 1-sum(cdf,2);

jj = 3.*ones(size(N));
for i = 1:length(D)
    for j = 1:length(N)
        prob(i,j) = N(j)*trapz(t,f(:,i).*(S.^(N(j)-1)));
        prob1(i,j) = jj(j)*nchoosek(N(j),jj(j))*...
            trapz(t,f(:,i).*(1-S).^(jj(j)-1).*(S.^(N(j)-jj(j))));
    end
end

figure;
for i = 2:length(D)
    loglog(N,prob(i,:),'Color',colors{i},'LineWidth',3);
    hold on
    loglog(N,prob1(i,:),'-.','Color',colors{i},'LineWidth',3);
end
hold off
x = xlabel('$$N$$');
set(x,'Interpreter','latex')
xlim([min(N) max(N)]);
xticks([10 10^2 10^3 10^4 10^5 10^6]);
xticklabels({'','10^2','','10^4','','10^6'});
set(gca,'XMinorTick','off');
xtickangle(0);
y = ylabel(strcat('Prob','$$(S_{n(1)} = s)$$'));
set(y,'Interpreter','latex')
set(gca,'fontsize',18)

figure;
for i = 1:length(D)
    semilogx(N,prob(i,:),'Color',colors{i},'LineWidth',3);
    hold on
    semilogx(N,prob1(i,:),'-.','Color',colors{i},'LineWidth',3);
end
hold on
hold off
x = xlabel('$$N$$');
set(x,'Interpreter','latex')
xlim([min(N) 10^(3.5)]);
xticks([10 10^2 10^3]);
xticklabels({'10^1','10^2','10^3'});
set(gca,'XMinorTick','off');
xtickangle(0);
y = ylabel(strcat('Prob','$$(S_{n(1)} = s)$$'));
ylim([0 1]);
yticks([0 .2 .4 .6 .8 1]);
yticklabels({'0','0.2','0.4','0.6','0.8','1'})
set(gca,'YMinorTick','off');
ytickangle(0);
set(y,'Interpreter','latex')
set(gca,'fontsize',18)
clear leg;
leg = {'$$s=1/16$$','','$$s=1/8$$','','$$s=1/4$$','','$$s=1/2$$',''};
l = legend(leg,'Location','northwest');
set(l,'Interpreter','latex')
