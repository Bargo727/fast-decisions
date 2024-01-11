% Samantha Linn
% u1323207@gcloud.utah.edu
% 2023
% Code for Figure 2A

D = 1;                                   % Diffusivity
Mu = [-1 0 1];                           % Rightward drift
N = round(logspace(1,6,100));            % Number and range of agents
J = 1;                                   % jth quickest agent
th = 1/2;                                % Interval half length
WW = th.*[-1/10 -1/20 0 1/20 1/10]';     % Initial positions of agents
q = ones(size(WW))./length(WW);          % Initial position probabilities
tmax = 1;                                % Cut off for short-time/long-time

kappa = 1000;                            % Number of terms in series approx
terms = 1000;                            % Number of terms in quadrature
tMax = 100;                              % Max time used in quadrature

colors = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.6350 0.0780 0.1840],...
    [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.3010 0.7450 0.9330],...
    [0.8500 0.3250 0.0980]};

prob = ones(length(WW),length(N),length(Mu));
bdy = zeros(length(WW),length(N),length(Mu));
for kk = 1:length(Mu)

    mu = Mu(kk);
    D = D/(2*th)^2;
    mu = mu/(2*th);
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
            prob(i,j,kk) = jj(j)*nchoosek(N(j),jj(j))*...
                trapz(t,f(:,i).*(1-S).^(jj(j)-1).*(S.^(N(j)-jj(j))));
            bdy(i,j,kk) = jj(j)*nchoosek(N(j),jj(j))*...
                trapz(t,fs1(:,i).*(1-S).^(jj(j)-1).*(S.^(N(j)-jj(j))));
        end
    end

end

figure;
subplot(1,3,1)
for i = 1:length(W)
    loglog(N,(bdy(i,:,3)./prob(i,:,3)),'Color',colors{i},'LineWidth',2);
    hold on
end
hold off
x = xlabel('$$N$$');
set(x,'Interpreter','latex')
xlim([min(N) max(N)]);
y = ylabel('Prob(fastest agent decides w/ pos. DRIFT cond. on init bias)');
y = ylabel('Prob$(sgn(X_{n(1)}(T_1)) = sgn(\mu)\,| \,X_{n(1)}(0) = y)$');
set(y,'Interpreter','latex')
leg = {'$$y=-\theta/10$$','$$y=-\theta/20$$','$$y=0$$',...
    '$$y=\theta/20$$','$$y=\theta/10$$'};
l = legend(leg,'Location','southwest');
set(l,'Interpreter','latex')
set(gca,'fontsize',15)

subplot(1,3,2)
for i = 1:length(W)
    if WW(i) > 0
        semilogx(N,(bdy(i,:,3)./prob(i,:,3)),'Color',colors{i},'LineWidth',2);
        hold on
    elseif WW(i) < 0
        semilogx(N,1-(bdy(i,:,3)./prob(i,:,3)),'Color',colors{i},'LineWidth',2);
        hold on
    end
end
hold off
x = xlabel('$$N$$');
set(x,'Interpreter','latex')
xlim([min(N) max(N)]);
y = ylabel('Prob$(sgn(X_{n(1)}(0)) = sgn(X_{n(1)}(T_1))\,| \,X_{n(1)}(0) = y)$');
set(y,'Interpreter','latex')
leg = {'$$y=-\theta/10$$','$$y=-\theta/20$$',...
    '$$y=\theta/20$$','$$y=\theta/10$$'};
l = legend(leg,'Location','southeast');
set(l,'Interpreter','latex')
set(gca,'fontsize',15)

subplot(1,3,3)
for i = 1:length(W)
    if WW(i) > 0
        loglog(N,1-(bdy(i,:,3)./prob(i,:,3)),'Color',colors{i},'LineWidth',2);
        hold on
    elseif WW(i) < 0
        loglog(N,(bdy(i,:,3)./prob(i,:,3)),'Color',colors{i},'LineWidth',2);
        hold on
    end
end
hold off
x = xlabel('$$N$$');
set(x,'Interpreter','latex')
xlim([min(N) max(N)]);
y = ylabel('Prob$(sgn(X_{n(1)}(0)) \neq sgn(X_{n(1)}(T_1))\,| \,X_{n(1)}(0) = y)$');
set(y,'Interpreter','latex')
leg = {'$$y=-\theta/10$$','$$y=-\theta/20$$',...
    '$$y=\theta/20$$','$$y=\theta/10$$'};
l = legend(leg,'Location','southwest');
set(l,'Interpreter','latex')
set(gca,'fontsize',15)
