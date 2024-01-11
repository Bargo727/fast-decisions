%this script will produce a figure showing how accuracy of decision varies
%with order of decision. We assume a prior of 1/2 for each state so that as
%N \to \infty, the accuracy of the first decider converges to 1/2.
N = 1000000;              % Number of agents
th = 1;                  % Threshold
M = 1000000;              % Number of trials
D  = 1;                  % Diffusivity
dt = 0.001;              % Time step
count = zeros(N,1);      % Stores probabilities
IC = (th/3)*ones(N,1);   % Inital condition

for j = 1:M

    y = IC;
    t = 0;
    idx = zeros(N,1);
    decision_time = zeros(N,1);
    mu = 1 - 2*round(rand);                       % Drift

    while sum(abs(y)>=th) < N
        ind = abs(y)<th;                          % Undecided agents
        y(ind) = y(ind) + dt*mu + sqrt(2*D*dt).*randn(sum(ind),1);
        t = t + dt;
        idx = (abs(y)>=th) - (decision_time~=0);  % Newly decided agents      
        decision_time = decision_time + t*idx;
    end

    B = sortrows([y,decision_time],2);            % Sort with times ascending
    count = count + double(sign(B(:,1))==sign(mu));

end

probability = count/M;

theory = 1/(1 + exp(-th))*ones(length(probability,1));



% save('~/matlab_output/probability.mat','probability')
% save('~/matlab_output/theory.mat','theory')
%
figure;
plot(1:N,probability,'-',1:N,theory,'k--','LineWidth',3)
set(gca,'fontsize',20)
xlabel('decider')
ylabel('accuracy')
