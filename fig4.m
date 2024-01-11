%
% sim_tri_dichot_ostats.m
%
% run a large N simulation to obtain the order statistics for the decision
% times of a large group split into biased and unbiased.
%
% here we consider a group making a decision from 3 alternatives. beliefs
% described by LLRs thus evolve on a planar domain with triangular boundary
%

NN = 1000;    % max number of agents
trials = 100; %number of trials
Probb = zeros(NN,1); %storing probability biased agent decides first
Probb1 = zeros(NN,1);
Probb2 = zeros(NN,1);
Prob0 = zeros(NN,1); %storing probability unbiased agent decides first


h = 1/2/sqrt(3);      % threshold parameter (apothem of eq triangle)
b1 = 1/8;
b2 = 1/4;
b = 1/2;    % strength of bias (fraction of the apothem toward choice 3)
a = 1/4;    % probability of agents being biased
i = sqrt(-1);   % imaginary unit

dt = 1e-4;  % timestep

% s = 2*sqrt(3)*h; % side length of triangle
l = cos(pi/6); r=sin(pi/6); % long and short side of 30/60/90 tri

% drift will just be up (alt 1 is correct); now construct covariance matrix

for N = 10:NN
    N
    

countb = 0;
countb1 = 0;
countb2 = 0;
count0 = 0;

for jj = 1:trials


% initialize belief value initial conditions (do in complex plane)
R = rand(N,1);
nb1 = sum(R<a);
nb2 = sum(R>a & R<2*a);
nb  = sum(R>2*a & R < 3*a);
Xb = h*b*(l-r*i)*ones(nb,1);Xb1 = h*b1*(l-r*i)*ones(nb1,1); Xb2 = h*b2*(l-r*i)*ones(nb2,1); X0 = zeros(N-nb-nb1-nb2,1);  % make biased and unbiased agents
db = zeros(nb,1); d0 = zeros(N-nb,1);   % make vector to store decision
Tb = zeros(nb,1); Tb1 = zeros(nb1,1); Tb2 = zeros(nb2,1); T0 = zeros(N-nb-nb1-nb2,1);   % make vector to store times

% simulate biased agents
for k=1:nb, j=1;
    x = real(Xb(k,j)); y = imag(Xb(k,j));
    while y<h && y>-2*h*(3*x+1) && y>2*h*(3*x-1)
        c = randn/sqrt(2);
        Xb(k,j+1)=Xb(k,j)+i*dt+sqrt(2*dt)*(randn+c+(randn+c)*i)/4;
        j=j+1;
        x = real(Xb(k,j)); y = imag(Xb(k,j));
    end
    Tb(k) = j;
    % if y<2*h*(3*x-1), db(k)=3; end
    % if y<-2*h*(3*x+1), db(k)=2; end
    % if y>h, db(k)=1; end
end

for k=1:nb1, j=1;
    x = real(Xb1(k,j)); y = imag(Xb1(k,j));
    while y<h && y>-2*h*(3*x+1) && y>2*h*(3*x-1)
        c = randn/sqrt(2);
        Xb1(k,j+1)=Xb1(k,j)+i*dt+sqrt(2*dt)*(randn+c+(randn+c)*i)/4;
        j=j+1;
        x = real(Xb1(k,j)); y = imag(Xb1(k,j));
    end
    Tb1(k) = j;
    % if y<2*h*(3*x-1), db(k)=3; end
    % if y<-2*h*(3*x+1), db(k)=2; end
    % if y>h, db(k)=1; end
end

for k=1:nb2, j=1;
    x = real(Xb2(k,j)); y = imag(Xb2(k,j));
    while y<h && y>-2*h*(3*x+1) && y>2*h*(3*x-1)
        c = randn/sqrt(2);
        Xb2(k,j+1)=Xb2(k,j)+i*dt+sqrt(2*dt)*(randn+c+(randn+c)*i)/4;
        j=j+1;
        x = real(Xb2(k,j)); y = imag(Xb2(k,j));
    end
    Tb2(k) = j;
    % if y<2*h*(3*x-1), db(k)=3; end
    % if y<-2*h*(3*x+1), db(k)=2; end
    % if y>h, db(k)=1; end
end

% simulate unbiased agents
for k=1:N-nb-nb1-nb2, j=1;
    x = real(X0(k,j)); y = imag(X0(k,j));
    while y<h && y>-2*h*(3*x+1) && y>2*h*(3*x-1)
        c = randn/sqrt(2);
        X0(k,j+1)=X0(k,j)+i*dt+sqrt(2*dt)*(randn+c+(randn+c)*i)/4;
        j=j+1;
        x = real(X0(k,j)); y = imag(X0(k,j));
    end
    T0(k) = j;
    % if y<2*h*(3*x-1), d0(k)=3; end
    % if y<-2*h*(3*x+1), d0(k)=2; end
    % if y>h, d0(k)=1; end
end

minb = min(Tb);
minb1 = min(Tb1);
minb2 = min(Tb2);
min0 = min(T0);

mint = min([minb,minb1,minb2,min0]);

if(mint == minb)
    countb = countb + 1;
elseif mint == min0
    count0 = count0 + 1;
elseif mint == minb1
    countb1 = countb1 + 1;
elseif mint == minb2
    countb2 = countb2 + 1;
end

end

Probb(N) = countb/trials;
Probb1(N) = countb1/trials;
Probb2(N) = countb2/trials;
Prob0(N) = count0/trials;
end

figure(1)
semilogx(10:NN,Prob0(10:end),10:NN,Probb1(10:end),10:NN,Probb2(10:end),10:NN,Probb(10:end),'LineWidth',3)
h = legend('unbiased','1/8','1/4','1/2');
set(h,'box','off')
set(gca,'fontsize',20)
xlabel('N')
ylabel('biased agent decides first')

