% Samantha Linn
% u1323207@gcloud.utah.edu
% 2023
% Code for Figure 1B

d = 1;                                   % Diffusivity
Mu = 1.5;                                % Rightward drift
N = round(logspace(1,5,500));            % Number and range of agents
J = 1;                                   % jth quickest agent
Th = 1;                                  % Interval half length
WW = Th.*[-1/4 -1/8 0 1/8 1/4 1/2]';     % Initial positions of agents
q = ones(size(WW))./length(WW);          % Initial position probabilities
tmax = 1;                                % Cut off for short-time/long-time

kappa = 1000;                            % Number of terms in series approx
terms = 1000;                            % Number of terms in quadrature
tMax = 100;                              % Max time used in quadrature

colors = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.6350 0.0780 0.1840],...
    [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.3010 0.7450 0.9330],...
    [0.8500 0.3250 0.0980]};

th = Th;
D = d/(2*Th)^2;
mu = Mu/(2*Th);
W = (WW + Th)./(2*Th);
L = ceil(W(end)) - floor(W(1));
jj = J.*ones(size(N));

prob = ones(length(W),length(N));
bdy = zeros(length(W),length(N));
bdy0 = zeros(length(W),length(N));
ts = logspace(-10,log10(tmax),terms)';
tl = logspace(log10(tmax),log10(tMax),terms)';
t = [ts; tl];

leg = cell(1,length(W));
fs0 = zeros(length(t),length(W));
fs1 = zeros(length(t),length(W));
cdf = zeros(length(t),length(W));

beta = zeros(size(W));
betab = zeros(size(W));
eta = zeros(size(W));
etab = zeros(size(W));
ab = zeros(length(W),length(N));

Len = zeros(size(W));
for i = 1:length(W)
    Len(i) = min(L-W(i),W(i));
end
mb = find(Len==min(Len));

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
    leg{i} = strcat('$$y=\,\,$$',num2str(WW(i)/th),'$$\theta$$');

    if ~any(mb==i) && J(1) == 1
        beta(i) = (Len(i)/Len(mb(1)))^2;
        betab(i) = ((L-w)/Len(mb(1)))^2;
        if Len(i) == L/2
            eta(i) = q(i)*(q(mb(1)))^(-beta(i))*exp(sqrt(beta(i))*(mu*Len(mb(1))-mu*Len(i))/(2*D))*...
                sqrt(beta(i)*pi^(beta(i)-1))*gamma(beta(i)) +...
                q(i)*(q(mb(1)))^(-beta(i))*exp(sqrt(beta(i))*(mu*(-Len(mb(1))-Len(i)))/(2*D))*...
                sqrt(beta(i)*pi^(beta(i)-1))*gamma(beta(i));
            etab(i) = q(i)*(q(mb(1)))^(-betab(i))*exp(sqrt(betab(i))*(mu*Len(mb(1))-mu*(L-w))/(2*D))*...
                sqrt(betab(i)*pi^(betab(i)-1))*gamma(betab(i)) +...
                q(i)*(q(mb(1)))^(-betab(i))*exp(sqrt(betab(i))*(mu*(-Len(mb(1))-(L-w)))/(2*D))*...
                sqrt(betab(i)*pi^(betab(i)-1))*gamma(betab(i));
        else
            eta(i) = q(i)*(q(mb(1)))^(-beta(i))*exp(sqrt(beta(i))*(sign(WW(i))*mu*Len(mb(1))-mu*Len(i))/(2*D))*...
                sqrt(beta(i)*pi^(beta(i)-1))*gamma(beta(i));
            etab(i) = q(i)*(q(mb(1)))^(-betab(i))*exp(sqrt(betab(i))*(sign(WW(i))*mu*Len(mb(1))-mu*(L-w))/(2*D))*...
                sqrt(betab(i)*pi^(betab(i)-1))*gamma(betab(i));
        end
    end

end
S = 1-sum(cdf,2);

for i = 1:length(W)
    for j = 1:length(N)
        prob(i,j) = jj(j)*nchoosek(N(j),jj(j))*...
            trapz(t,f(:,i).*(1-S).^(jj(j)-1).*(S.^(N(j)-jj(j))));
        bdy(i,j) = jj(j)*nchoosek(N(j),jj(j))*...
            trapz(t,fs1(:,i).*(1-S).^(jj(j)-1).*(S.^(N(j)-jj(j))));
        bdy0(i,j) = jj(j)*nchoosek(N(j),jj(j))*...
            trapz(t,fs0(:,i).*(1-S).^(jj(j)-1).*(S.^(N(j)-jj(j))));
    end
end

figure;
for i = 1:length(W)
    if ~any(mb==i)
        loglog(N,prob(i,:),'Color',colors{i},'LineWidth',3);
        hold on
        asymp = eta(i).*(log(N)).^((beta(i)-1)/2).*N.^(1-beta(i));
        loglog(N,asymp,'--','Color',colors{i},'LineWidth',2)
        ab(i,:) = abs(prob(i,:)-asymp);
    end
end
hold off
x = xlabel('$$N$$');
set(x,'Interpreter','latex')
xlim([min(N) max(N)]);
xticks([10 10^2 10^3 10^4 10^5 10^6]);
xticklabels({'','10^2','','10^4','','10^6'});
set(gca,'XMinorTick','off');
xtickangle(0);
y = ylabel(strcat('Prob','$$(X_{n(1)}(0) = y)$$'));
set(y,'Interpreter','latex')
set(gca,'fontsize',18)

figure;
for i = 1:length(W)
    if ~any(mb==i)
        semilogx(N,prob(i,:),'Color',colors{i},'LineWidth',2);
        hold on
    end
end
semilogx(N,prob(mb,:),'Color','k','LineWidth',2);
hold on
hold off
x = xlabel('$$N$$');
set(x,'Interpreter','latex')
xlim([min(N) 10^(3.5)]);
xticks([10 10^2 10^3]);
xticklabels({'10^1','10^2','10^3'});
set(gca,'XMinorTick','off');
xtickangle(0);
y = ylabel(strcat('Prob','$$(X_{n(1)}(0) = y)$$'));
ylim([0 1]);
yticks([0 .2 .4 .6 .8 1]);
yticklabels({'0','0.2','0.4','0.6','0.8','1'})
set(gca,'YMinorTick','off');
ytickangle(0);
set(y,'Interpreter','latex')
set(gca,'fontsize',18)
clear leg;
leg = {'$$y=-\theta/4$$','$$y=-\theta/8$$','$$y=0$$',...
    '$$y=\theta/8$$','$$y=\theta/4$$','$$y=\theta/2$$'};
l = legend(leg,'Location','northwest');
set(l,'Interpreter','latex')

figure;
subplot(1,2,1);
for i = 1:length(W)
    if ~any(mb==i)
        loglog(N,prob(i,:),'Color',colors{i},'LineWidth',3);
        hold on
        asymp = eta(i).*(log(N)).^((beta(i)-1)/2).*N.^(1-beta(i));
        loglog(N,asymp,'--','Color',colors{i},'LineWidth',3)
        ab(i,:) = abs(prob(i,:)-asymp);
    end
end
hold off
x = xlabel('$$N$$');
set(x,'Interpreter','latex')
xlim([min(N) max(N)]);
xticks([10 10^2 10^3 10^4 10^5 10^6]);
xticklabels({'','10^2','','10^4','','10^6'});
set(gca,'XMinorTick','off');
xtickangle(0);
y = ylabel(strcat('Prob','$$(X_{n(1)}(0) = y)$$'));
set(y,'Interpreter','latex')
set(gca,'fontsize',18)
clear leg;
leg = {'$$y=-\theta/4$$','','$$y=-\theta/8$$','','$$y=0$$','',...
    '$$y=\theta/8$$','','$$y=\theta/4$$',''};
l = legend(leg,'Location','southwest');
set(l,'Interpreter','latex')


subplot(1,2,2)
for i = 1:length(W)
    if ~any(mb==i)
        loglog(N,ab(i,:)./prob(i,:),'Color',colors{i},'LineWidth',3);
        hold on
    end
end
hold off
x = xlabel('$$N$$');
set(x,'Interpreter','latex')
xlim([min(N) max(N)]);
xticks([10 10^2 10^3 10^4 10^5 10^6]);
xticklabels({'','10^2','','10^4','','10^6'});
set(gca,'XMinorTick','off');
xtickangle(0);
y = ylabel('Relative error');
set(y,'Interpreter','latex')
set(gca,'fontsize',18)
