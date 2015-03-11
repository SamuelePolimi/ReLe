%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference: Deisenroth, M. P.; Neumann, G.; Peters, J. (2013),
% A Survey on Policy Search for Robotics, Foundations and Trends in
% Robotics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change the basis functions if you wish. Also, you can use different bases
% for the upper-level policy and for the value function.

% In this implementation the contexts are the goals.

clear all
domain = 'deep';
[n_obj, pol_low, episodes, steps, gamma, avg_rew_setting, max_obj] = settings(domain);

% Parameters of the (Gaussian) distribution (upper-level policy) that generates parameters theta
n_params = size(pol_low.theta,1);

% If the low-level policy has a learnable variance, we don't want to learn 
% it and we make it deterministic
n_params = n_params - pol_low.dim_variance_params;
pol_low = pol_low.makeDeterministic;

dim_basis_policy = CREPS_basis;
dim_basis_vfun = CREPS_basis;
W0 = zeros(n_params,dim_basis_policy);
mu0 = zeros(n_params,1);
Sigma0 = 100 * eye(n_params); % change according to the domain
% pol_high = smart_gaussian_policy(@CREPS_basis, n_params, W0, Sigma0);
% pol_high = smart_diag_gaussian_policy(@CREPS_basis, n_params, W0, sqrt(diag(Sigma0)));
pol_high = full_smart_gaussian_policy(@CREPS_basis, n_params, mu0, W0, Sigma0);

N = 50;
epsilon = 0.9;

iter = 0;

%% Learning
while true

    % Cells used to save J for each context during sampling
    results = cell(n_obj,1);
    iter = iter + 1;
    
    Phi_policy = zeros(N,dim_basis_policy);
    Phi_vfun = zeros(N,dim_basis_vfun);
    Theta = zeros(N,n_params);
    J = zeros(1,N);

    % Sample theta
    for k = 1 : N % to use the parfor you need to not save 'results{.}'
        % Generate context
        robj = randi(n_obj);
        context = [robj == 1 : n_obj]';
        Phi_policy(k,:) = CREPS_basis(context);
        Phi_vfun(k,:) = CREPS_basis(context); % use different bases if you want
        
        % Extract theta from the Gaussian
        theta = pol_high.drawAction(context);
        pol_low.theta(1:n_params) = theta; % set only the mean, not the variance
        Theta(k,:) = theta;

        % Collect samples with one theta
        [ds, uJ, dJ] = collect_samples(domain,episodes,steps,pol_low,avg_rew_setting,gamma);
        if gamma == 1
            J(k) = uJ(robj) .* max_obj(robj);
        else
            J(k) = dJ(robj) .* max_obj(robj);
        end
        results{robj} = [results{robj}, J(k)];
    end

    % Print a matrix with [r_obj, n_samples per context, mean_J of context]
    print_vector = zeros(n_obj, 3);
    for i = 1 : n_obj
        print_vector(i,:) = [i, numel(results{i}), mean(results{i})];
    end
    fprintf('\n%g) Context \t N. samples \t Mean Reward \n', iter)
    fprintf('%f \t %d \t\t\t %f\n', print_vector')
    
    % Optimization problem settings
    options = optimset('GradObj', 'on', 'Display', 'off', ...
        'MaxFunEvals', 10 * 5, 'Algorithm', 'interior-point', ...
        'TolX', 10^-8, 'TolFun', 10^-12, 'MaxIter', 10);
    lowerBound_v = -ones(dim_basis_vfun, 1) * 1e8;
    upperBound_v = ones(dim_basis_vfun, 1) * 1e8;
    lowerBound_eta = 1e-8;
    upperBound_eta = 1e8;
    v = ones(dim_basis_vfun,1); % init params
    eta = 1;

    maxIter = 100;
    validKL = false;
    validSF = false;
    numStepsNoKL = 0;
    
    J_tmp = J;
    
    % Iteratively solve fmincon for eta and theta separately
    for i = 1 : maxIter
        
        J = J_tmp;
        
        if ~validKL || numStepsNoKL < 5
            eta = fmincon(@(eta)CREPS_duals.eta(eta,v,J,epsilon,Phi_vfun), ...
                eta, [], [], [], [], lowerBound_eta, upperBound_eta, [], options);

            % Numerical trick
            advantage = J_tmp - v' * Phi_vfun';
            maxAdvantage = max(advantage);
            J_tmp = J_tmp - maxAdvantage;

            % Perform weighted ML to update high-level policy parameters
            d = exp( (J - v' * Phi_vfun') / eta )';

            % Check conditions
            qWeighting = ones(N,1);
            pWeighting = d;
            pWeighting = pWeighting / sum(pWeighting);
            divKL = getKL(pWeighting, qWeighting);
            error = divKL - epsilon;
            validKL = error < 0.1 * epsilon;
            featureDiff = sum(bsxfun(@times, Phi_vfun, pWeighting)) - mean(Phi_vfun);
            validSF = max(abs(featureDiff)) < 0.1;
            numStepsNoKL = numStepsNoKL + 1;
        end

        if ~validSF
            v = fmincon(@(v)CREPS_duals.theta(v,eta,J,epsilon,Phi_vfun), ...
                v, [], [], [], [], lowerBound_v, upperBound_v, [], options);
        
            % Numerical trick
            advantage = J - v' * Phi_vfun';
            maxAdvantage = max(advantage);
            J = J - maxAdvantage;

            % Perform weighted ML to update high-level policy parameters
            d = exp( (J - v' * Phi_vfun') / eta )';

            % Check conditions
            qWeighting = ones(N,1);
            pWeighting = d;
            pWeighting = pWeighting / sum(pWeighting);
            divKL = getKL(pWeighting, qWeighting);
            error = divKL - epsilon;
            validKL = error < 0.1 * epsilon;
            featureDiff = sum(bsxfun(@times, Phi_vfun, pWeighting)) - mean(Phi_vfun);
            validSF = max(abs(featureDiff)) < 0.1;
        end
        
        if validSF && validKL
            break
        end

    end

    % Stopping condition
    if divKL < 1e-3
        return
    end
    
    pol_high = pol_high.weightedMLUpdate(d,Theta,Phi_policy);
    
end