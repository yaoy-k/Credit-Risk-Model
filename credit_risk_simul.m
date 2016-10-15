clear all;
clc
format long;

Nout  = 100000; % number of out-of-sample scenarios
Nin   = 5000;   % number of in-sample scenarios
Ns    = 5;      % number of idiosyncratic scenarios for each systemic

C = 8;          % number of credit states

% Filename to save out-of-sample scenarios
filename_save_out  = 'scen_out';

% Read and parse instrument data
instr_data = dlmread('instrum_data.csv', ',');
instr_id   = instr_data(:,1);           % ID
driver     = instr_data(:,2);           % credit driver
beta       = instr_data(:,3);           % beta (sensitivity to credit driver)
recov_rate = instr_data(:,4);           % expected recovery rate
value      = instr_data(:,5);           % value
prob       = instr_data(:,6:6+C-1);     % credit-state migration probabilities (default to A)
exposure   = instr_data(:,6+C:6+2*C-1); % credit-state migration exposures (default to A)
retn       = instr_data(:,6+2*C);       % market returns

K = size(instr_data, 1); % number of  counterparties

% Read matrix of correlations for credit drivers
rho = dlmread('credit_driver_corr.csv', '\t');
sqrt_rho = (chol(rho))'; % Cholesky decomp of rho (for generating correlated Normal random numbers)

disp('======= Credit Risk Model with Credit-State Migrations =======')
disp('============== Monte Carlo Scenario Generation ===============')
disp(' ')
disp(' ')
disp([' Number of out-of-sample Monte Carlo scenarios = ' int2str(Nout)])
disp([' Number of in-sample Monte Carlo scenarios = ' int2str(Nin)])
disp([' Number of counterparties = ' int2str(K)])
disp(' ')

% Find credit-state for each counterparty
% 8 = AAA, 7 = AA, 6 = A, 5 = BBB, 4 = BB, 3 = B, 2 = CCC, 1 = default
[Ltemp, CS] = max(prob, [], 2);
clear Ltemp

% Account for default recoveries
exposure(:, 1) = (1-recov_rate) .* exposure(:, 1);

% Compute credit-state boundaries
CS_Bdry = norminv( cumsum(prob(:,1:C-1), 2) );

% -------- Insert your code here -------- %

if(~exist('scenarios_out.mat','file'))
scenarios_out = normrnd(0,1, [Nout,50])*sqrt_rho; % Generate 100000*50 of the credit drivers
Losses_out = ones(Nout,100); % Pre-generate the size of the Losses_out
ynew = ones(K,1); % Pre-generate the size of the individual systemic error

 for s = 1:Nout
     ysc = scenarios_out(s,:); % Select the credit drvier (y) corresponding to scenario
     z = normrnd(0,1,[100,1]); % Generate 100*1 idiosyncratic Factor
     for j = 1:K
         d = driver(j);
         ynew(j,1) = ysc(d); % Assign correspoding credit driver to each counter parties
     end
CS_Bdry(:,8) = 10000; % Included the eighth boundry for any couterparty that has w > CS_boundary(:,7)
    %Determine the position of w in the given boundary 
     wout = (beta.*ynew)+((sqrt(1-(beta.^2)))).*z; % Creditworthiness indexes of counterparty 100*1 matrix

       C = bsxfun(@minus,CS_Bdry,wout); %Calcuate the differences between the CS_Bdry and Creditworthiness index
       C(C<0)=9000; % Change any negative value to a really large Number so that using minimum, will give the lowest value of the positive number.
       [temp, w_cs] = min(C,[],2); % The position of the couterparty based on the new credit state

   for j = 1:K
 Losses_out(s,j) = exposure(j,w_cs(j)); % The losses of the counterparty based on the new credit state. (100000*100)
    end
 end
    
    % Calculated out-of-sample losses (100000 x 100)
    % Losses_out

    save('scenarios_out', 'Losses_out')
else
    load('scenarios_out', 'Losses_out')
end

% Normal approximation computed from out-of-sample scenarios
mu_l = mean(Losses_out)'; % Average mean losses
var_l = cov(Losses_out); % Covariance of the losses

% Compute portfolio weights
portf_v = sum(value);     % portfolio value
w0{1} = value / portf_v;  % asset weights (portfolio 1)
w0{2} = ones(K, 1) / K;   % asset weights (portfolio 2)
x0{1} = (portf_v ./ value) .* w0{1};  % asset units (portfolio 1)
x0{2} = (portf_v ./ value) .* w0{2};  % asset units (portfolio 2)

% Quantile levels (99%, 99.9%)
alphas = [0.99 0.999];
tic
% Compute VaR and CVaR (non-Normal and Normal) for 100000 scenarios
for(portN = 1:2)
    for(q=1:length(alphas))
        alf = alphas(q);

        l{portN,q} = sort(Losses_out * x0{portN}); % Sorting the portfolio losses
        VaRout(portN,q) = l{portN,q}(ceil(Nout*alf)); % Calculate the portfolio losses VaR under Non-Normal distribution
        VaRinN (portN,q) = (x0{portN}')*mu_l +norminv(alf,0,1)*(sqrt(x0{portN}'*var_l*x0{portN})); % Calculate the portfolio losses VaR under Normal distribution
        CVaRout(portN,q) = (1/(Nout*(1-alf)))*((ceil(Nout*alf)-Nout*alf) * VaRout(portN,q) + sum(l{portN,q}(ceil(Nout*alf)+1:Nout))); % Calculate the portfolio losses CVaR under Non-Normal distribution
        CVaRinN(portN,q) = (x0{portN}')*mu_l + (normpdf(norminv(alf,0,1))/(1-alf))*(sqrt(x0{portN}'*var_l*x0{portN})); % Calculate the portfolio losses CVaR under Normal distribution
    
    end
end


% Perform 100 trials
N_trials = 100;
for(tr=1:N_trials)
    
    % Monte Carlo approximation 1

 ymc1 = normrnd(0,1, [Nin/Ns,50])*sqrt_rho; % Generate 10000*50 of the credit drivers
 ymc1new = ones(K,1); % Pre-generate the size of the individual systemic error
 Losses_out2 = ones(Ns,100); % Pre-generate the size of the Losses_out

 
 for s = 1:ceil(Nin/Ns) % systemic scenarios
        % -------- Insert your code here -------- %
         ymc1l = ymc1(s,:); % Select the credit drvier (y) corresponding to scenario
     for j = 1:K
         d = driver(j);
         ymc1new(j,1) = ymc1l(d); % Assign correspoding credit driver to each counter parties
     end
             
       for si = 1:Ns % idiosyncratic scenarios for each systemic, 5 different z for each scenario
           z = normrnd(0,1,[100,1]); % Generate 100*1 idiosyncratic Factor
           CS_Bdry(:,8) = 10000;
       wmc1 = (beta.*ymc1new)+((sqrt(1-(beta.^2)))).*z; % Creditworthiness indexes of counterparty 100*1 matrix

       C = bsxfun(@minus,CS_Bdry,wmc1);  %Calcuate the differences between the CS_Bdry and Creditworthiness index
       C(C<0)=9000;  % Change any negative value to a really large Number so that using minimum, will give the lowest value of the positive number.
       [temp, wmc1] = min(C,[],2); % The position of the couterparty based on the new credit state
        for j = 1:K
        Losses_out2(si,j) = exposure(j,wmc1(j)); % The losses of the counterparty based on the new credit state. (5*100)
        end
       end

       if s == 1
       Losses_inMC1 = Losses_out2 ; % The losses of the counterparty based on the new credit state. (5*100)
       else
           Losses_inMC1 = [Losses_inMC1 ;Losses_out2]; % The losses of the counterparty based on the new credit state. (5000*100)
       end
       
 end
 
   

    
    % Calculated losses for MC1 approximation (5000 x 100)
    % Losses_inMC1
    
    % Monte Carlo approximation 2
    
    % -------- Insert your code here -------- %
    ymc2 = normrnd(0,1, [Nin,50])*sqrt_rho; % Generate 5000*50 of the credit drivers
    Losses_inMC2 = ones(Nin,100); % Pre-generate the size of the individual systemic error
    ymc2new = ones(K,1); % Pre-generate the size of the Losses_out
    
    for s = 1:Nin % systemic scenarios (1 idiosyncratic scenario for each systemic)
        % -------- Insert your code here -------- %
         ymc2sc = ymc2(s,:); % Select the credit drvier (y) corresponding to scenario
     z = normrnd(0,1,[100,1]); % Generate 100*1 idiosyncratic Factor
     for j = 1:K
         d = driver(j);
         ymc2new(j,1) = ymc2sc(d); % Assign correspoding credit driver to each counter parties
     end
        CS_Bdry(:,8) = 10000;
        wmc2 = (beta.*ymc2new)+((sqrt(1-(beta.^2)))).*z; % Creditworthiness indexes of counterparty 100*1 matrix

       C = bsxfun(@minus,CS_Bdry,wmc2); %Calcuate the differences between the CS_Bdry and Creditworthiness index
       C(C<0)=9000;
       [temp, w_mc2cs] = min(C,[],2);

    for j = 1:K
        Losses_inMC2(s,j) = exposure(j,w_mc2cs(j));
    end

    end
        
    % Calculated losses for MC2 approximation (5000 x 100)
    % Losses_inMC2

    % Compute VaR and CVaR
    for(portN = 1:2)
        for(q=1:length(alphas))
            alf = alphas(q);
            % -------- Insert your code here -------- %            
            % Compute portfolio loss 
            portf_loss_inMC1 = sort(Losses_inMC1 * x0{portN});
            portf_loss_inMC2 = sort(Losses_inMC2 * x0{portN});
             % Compute portfolio mean loss mu_p_MC1 and portfolio standard deviation of losses sigma_p_MC1
            % Compute portfolio mean loss mu_p_MC2 and portfolio standard deviation of losses sigma_p_MC2           
            mu_MCl = mean(Losses_inMC1)';
            var_MCl = cov(Losses_inMC1);
            mu_MC2 = mean(Losses_inMC2)';
            var_MC2 = cov(Losses_inMC2);

            % Compute VaR and CVaR for the current trial
            VaRinMC1{portN,q}(tr) = portf_loss_inMC1(ceil((Nin)*alf)); %VaR of Monte Carlo 1 under non- normal distribution
            VaRinMC2{portN,q}(tr) = portf_loss_inMC2(ceil((Nin)*alf)); %VaR of Monte Carlo 2 under non- normal distribution
            VaRinN1{portN,q}(tr) = (x0{portN}')*mu_MCl +norminv(alf,0,1)*(sqrt(x0{portN}'*var_MCl*x0{portN})); %VaR of Monte Carlo 1 under normal distribution
            VaRinN2{portN,q}(tr) = (x0{portN}')*mu_MC2 +norminv(alf,0,1)*(sqrt(x0{portN}'*var_MC2*x0{portN})); %VaR of Monte Carlo 2 under normal distribution
            CVaRinMC1{portN,q}(tr) = (1/(Nin*(1-alf)))*((ceil(Nin*alf)-Nin*alf) * VaRinMC1{portN,q}(tr) + sum(portf_loss_inMC1((ceil(Nin*alf)+1:Nin)))); %CVaR of Monte Carlo 1 under non- normal distribution
            CVaRinMC2{portN,q}(tr) = (1/(Nin*(1-alf)))*((ceil(Nin*alf)-Nin*alf) * VaRinMC2{portN,q}(tr) + sum(portf_loss_inMC2(ceil(Nin*alf)+1:Nin))); %CVaR of Monte Carlo 2 under non- normal distribution
            CVaRinN1{portN,q}(tr) = (x0{portN}')* mu_MCl + (normpdf(norminv(alf,0,1))/(1-alf))*(sqrt(x0{portN}'*var_MCl*x0{portN})); %CVaR of Monte Carlo 1 under normal distribution
            CVaRinN2{portN,q}(tr) = (x0{portN}')* mu_MC2 + (normpdf(norminv(alf,0,1))/(1-alf))*(sqrt(x0{portN}'*var_MC2*x0{portN})); %CVaR of Monte Carlo 2 under normal distribution
toc
        end
    end
end


% Display portfolio VaR and CVaR
for(portN = 1:2)
fprintf('\nPortfolio %d:\n\n', portN)    
 for(q=1:length(alphas))
    alf = alphas(q);
    fprintf('Out-of-sample: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, VaRout(portN,q), 100*alf, CVaRout(portN,q))
    fprintf('In-sample MC1: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinMC1{portN,q}), 100*alf, mean(CVaRinMC1{portN,q}))
    fprintf('In-sample MC2: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinMC2{portN,q}), 100*alf, mean(CVaRinMC2{portN,q}))
    fprintf(' In-sample No: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, VaRinN(portN,q), 100*alf, CVaRinN(portN,q))
    fprintf(' In-sample N1: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinN1{portN,q}), 100*alf, mean(CVaRinN1{portN,q}))
    fprintf(' In-sample N2: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n\n', 100*alf, mean(VaRinN2{portN,q}), 100*alf, mean(CVaRinN2{portN,q}))
 end
end

% Plot results
% figure(1);
figure(1)
 for(portN = 1:2)
        
    losses = sort(Losses_out * x0{portN});
    figure(portN);
    [frequencyCounts, binLocations] = hist(losses, 1000);
    bar(binLocations, frequencyCounts);
    xlabel('Frequency')
    ylabel('Losses')
    hold on;
    for(q = 1:length(alphas))
        line([mean(VaRout(portN,q)),mean(VaRout(portN,q))], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
        line([mean(VaRinN (portN,q)),mean(VaRinN (portN,q))], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')     
        hold on;
    end
    hold off;
   
    text(mean(VaRout(portN,1)), max(frequencyCounts), 'VaRout@99%'); 
    text(0.9*mean(VaRinN (portN,1)), max(frequencyCounts), 'VaRinN@99%');
    text(mean(VaRout(portN,2)), max(frequencyCounts)/2, 'VaRout@99.9%'); 
    text(0.9*mean(VaRinN (portN,2)), max(frequencyCounts)/1.5, 'VaRinN@99.9%');    
end


% figure(2);
for(portN = 1:2)
    figure(portN+1);
    [frequencyCounts, binLocations] = hist(portf_loss_inMC1, 1000);
    bar(binLocations, frequencyCounts); 
        xlabel('Frequency')
    ylabel('Losses')
    hold on;
    for(q = 1:length(alphas))
        line([mean(VaRinMC1{portN,q}),mean(VaRinMC1{portN,q})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
        line([mean(VaRinN1{portN,q}),mean(VaRinN1{portN,q})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
        hold on;
    end
    hold off;
   
    text(mean(VaRinMC1{portN,1}), max(frequencyCounts), 'VaRinMC1@99%'); 
    text(0.9*mean(VaRinN1{portN,1}), max(frequencyCounts), 'VaRinN1@99%');
    text(mean(VaRinMC1{portN,2}), max(frequencyCounts)/2, 'VaRinMC1@99.9%'); 
    text(0.9*mean(VaRinN1{portN,2}), max(frequencyCounts)/2, 'VaRinN1@99.9%');
end

for(portN = 1:2)
    figure(portN+3);
    [frequencyCounts, binLocations] = hist( portf_loss_inMC2, 1000);
    bar(binLocations, frequencyCounts); 
            xlabel('Frequency')
    ylabel('Losses')
    hold on;
    for(q = 1:length(alphas))
        line([mean(VaRinMC2{portN,q}),mean(VaRinMC2{portN,q})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
        line([mean(VaRinN2{portN,q}),mean(VaRinN2{portN,q})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
        hold on;
    end
    hold off;
   
    text(mean(VaRinMC2{portN,1}), max(frequencyCounts), 'VaRinMC2@99%'); 
    text(0.9*mean(VaRinN2{portN,1}), max(frequencyCounts), 'VaRinN2@99%');
    text(mean(VaRinMC2{portN,2}), max(frequencyCounts)/2, 'VaRinMC2@99.9%'); 
    text(0.9*mean(VaRinN2{portN,2}), max(frequencyCounts)/2, 'VaRinN2@99.9%');
end