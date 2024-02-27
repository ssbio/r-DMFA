

% r-DMFA function call with updated variable names
% output = DMFASampler(model, measuredMetabolites, concentrationData, stdMeasurement, minX, maxX, N, lambdaOpt, sdFracScale, robustFit, tolerance)
%% PLEASE NOTE:
% THIS CODE REQUIRES a Licensed BARON SOLVER!
% BARON itself uses IBM CPLEX which also requires a license. 
% While BARON Academic license can be acquired for a fee,
% CPLEX solver does provide academic free license easily obtainable.
%%
function output = DMFASampler_constr(model, measuredMetabolites, concentrationData, stdMeasurement, minX, maxX, N, lambdaOpt, sdFracScale, robustFit, tolerance, blockedRxnsIdx)
% **r-DMFASampler*
% % **Description:**
% r-DMFASampler is a MATLAB function that performs the regularized dynamic metabolic flux analysis (r-DMFA) on a given metabolic model using concentration data. 
% DMFA is a technique to estimate metabolic fluxes over time, considering the dynamic changes in metabolite concentrations. 
% The function fits a DMFA model to the measured concentration data, identifies the optimal number of DMFA time points, and simulates the time-dependent metabolic fluxes.
% 
% **Inputs:**
% - `model`: A struct representing the metabolic model. It should contain the following fields:
%   - `rawMetabName`: A cell array of metabolite names in the model.
%   - `rxnIDs`: A cell array of reaction IDs in the model.
%   - `S`: The stoichiometry matrix of the model.
% 
% - `measuredMetabolites`: A cell array containing the names of metabolites with measured concentrations.
% 
% - `concentrationData`: A matrix representing the measured metabolite concentrations over time. Rows represent metabolites, and columns represent time points.
% 
% - `stdMeasurement`: A matrix representing the standard deviation of the measured concentrations over time. It should have the same dimensions as `concentrationData`.
% 
% - `minX` and `maxX`: Scalars representing the minimum and maximum values for scaling the concentration data during the analysis.
% 
% - `N`: An integer representing the number of spline interpolation points for curve fitting.
% 
% - `lambdaOpt`: A scalar representing the regularization parameter for DMFA.
% 
% - `sdFracScale`: A scalar representing the scale factor for adding random error to the concentration data.
% 
% - `robustFit`: A boolean indicating whether to perform robust fitting (true) or not (false).
% 
% - `tol`: A scalar representing the convergence tolerance for DMFA fitting.
% 
% **Outputs:**
% The function returns an output struct containing the following fields:
% 
% - `figNo`: A figure number for plotting SSR vs. number of DMFA time points.
% 
% - `meas_est_x_vals`: A matrix representing the estimated metabolite concentrations over time after scaling and normalization.
% 
% - `meas_est_x_vals_sd`: A matrix representing the standard deviation of the estimated metabolite concentrations over time after scaling and normalization.
% 
% - `flux_data_n`: A matrix representing the normalized metabolic fluxes over time.
% 
% - `flux_data_sd_n`: A matrix representing the standard deviation of the normalized metabolic fluxes over time.
% 
% - `flux_data`: A matrix representing the metabolic fluxes over time.
% 
% - `flux_data_sd`: A matrix representing the standard deviation of the metabolic fluxes over time.
% 
% - `EndRunTime`: The total runtime of the function in minutes.
% 
% - `model`: The input metabolic model struct.
% 
% - `simdata`: A struct containing the simulated concentration data.
% 
% - `nDMFA`: The optimal number of DMFA time points.
% 
% - `K`, `c0`, `U`, `tn`, and `tm`: Parameters of the fitted DMFA model.
% 
% - `data`: The input data struct used in the DMFA fitting.
% 
% - `fit`: The DMFA fitting results.
% 
% - `fit_ssr`: The sum of squares of the residuals (SSR) for the fitted DMFA model.
% 
% **Usage Example:**
% 
% ```matlab
% % Load the metabolic model and concentration data
% model = load('metabolic_model.mat');
% measuredMetabolites = {'Met1', 'Met2', 'Met3'};
% concentrationData = load('concentration_data.mat');
% stdMeasurement = load('std_measurement.mat');
% 
% % Set the input parameters
% minX = 0;
% maxX = 1;
% N = 100;
% lambdaOpt = 0.01;
% sdFracScale = 0.1;
% robustFit = true;
% tol = 1e-5;
% 
% % Run r-DMFASampler
% output = DMFASampler(model, measuredMetabolites, concentrationData, stdMeasurement, minX, maxX, N, lambdaOpt, sdFracScale, robustFit, tol);
% 
% % Access the output variables
% disp(output.meas_est_x_vals);
% disp(output.flux_data);
% disp(output.EndRunTime);
% ```
% 
% **Note:**
% Please ensure that you have the required metabolic model, concentration data, 
% and standard deviation data files available and correctly formatted before using the DMFASampler function. 
% Also, the function requires the DMFA fitting and simulation functions (`dmfa` and `dmfa_sim_3`) to be available in the MATLAB environment.

% Adapted and modified from (Leighty & Antoniewicz, 2011) DMFA matlab base code
    % Start timing the execution
    startTime = tic;
    % Reset the lastwarn message and id
    lastwarn('', '');
    % Get the list of metabolites and reaction IDs from the model
    metaboliteIDs = model.rawMetabName;
    reactionIDs = model.rxnIDs;
    % Get the indices of measured metabolites in the model
    measuredMetaboliteIndices = zeros(numel(measuredMetabolites), 1);
    for k = 1:numel(measuredMetabolites)
        measuredMetaboliteIndices(k) = find(cellfun(@(x)isequal(x, measuredMetabolites{k}), metaboliteIDs));
    end
    % Get the number of time points in the concentration data
    timePointsMeasured = size(concentrationData, 2);
    % Spline curve fitting
    time = 0:1:timePointsMeasured-1;
    NPoints = (time(end)*N) + 1;
    timeInterpolated = linspace(0, time(end), NPoints);
    splineData = spline(time, concentrationData, timeInterpolated);
    splineData(splineData < 0) = 0;
    splineData(:, 1) = concentrationData(:, 1);
    splineStd = spline(time, stdMeasurement, timeInterpolated);
    splineStd(splineStd <= tolerance) = tolerance;
    nonMeasuredIndices = setdiff(1:length(metaboliteIDs), measuredMetaboliteIndices);
    modelMetabolites = model.rawMetabName;
    modelTypes = repmat({'nbal'}, size(modelMetabolites));
    modelTypes(nonMeasuredIndices(1:end-1)) = {'bal'};
    timeInterpolatedSelection = false(1, length(timeInterpolated));
    [~, indicesDays] = ismember(time, timeInterpolated);
    timeInterpolatedSelection(:, indicesDays) = true(1, numel(indicesDays));
    selection = timeInterpolatedSelection;
    type = repmat({'metabolite'}, size(measuredMetabolites));
    standardDeviation = num2cell(splineStd, 2);
    noisyData = splineData + sdFracScale.*splineStd.*(-1 + (1+1)*randn(size(splineData)));
    noisyData(noisyData < 0) = 0;
    values = num2cell(noisyData, 2);
    data = struct('type', type, 'met', measuredMetabolites, 'sel', selection, 'val', values, 'sd', standardDeviation);
    % Metabolic model
    nonBalancedIndices = strcmp(modelTypes, 'nbal');
    balancedIndices = strcmp(modelTypes, 'bal');
    SnonBalanced = model.S(nonBalancedIndices, :);
    Sbalanced = model.S(balancedIndices, :);
    K = null(Sbalanced, 'r');
    model.type = modelTypes ;
    model.met = modelMetabolites ;     % list of metabolites in the model  
    % DMFA fit
    lambda = lambdaOpt;
    fit = dmfa(model, data(1:end), timeInterpolated, lambda, blockedRxnsIdx);
    % Find the best fit with the best DMFA time points
    [nDMFA, figNo] = findBestFit(fit, tolerance, robustFit);
    output.figNo = figNo;
%     fprintf('nDMFA: %d\n', nDMFA);
    % Simulate to get concentrations of DMFA analysis
    K = fit(nDMFA).K;
    c0 = fit(nDMFA).c0;
    U = fit(nDMFA).U;
    tn = fit(nDMFA).tn;
    simData = dmfa_sim_3(model, data(1:end), timeInterpolated, K, c0, U, tn);
    % Concentration profiles
    % Extracting data
    measuredMetabolites = {data(1:end).met}';
    estimatedConcentrations = cell2mat({simData(1:end).sim}');
    estimatedConcentrations(estimatedConcentrations < 0) = 0;
    estimatedConcentrationsStd = cell2mat({simData(1:end).sd}');
    % Interpolation of values
    estimatedConcentrations = spline(timeInterpolated(selection), estimatedConcentrations, timeInterpolated);
    estimatedConcentrations(estimatedConcentrations < 0) = 0;
    estimatedConcentrationsStd = spline(timeInterpolated, estimatedConcentrationsStd, timeInterpolated);
    estimatedConcentrationsStd(estimatedConcentrationsStd < 0) = 0;
    % Scaling and normalization of data
    concDataScaled = (maxX - minX) .* concentrationData + minX;
    estimatedConcentrationsScaled = (maxX - minX) .* estimatedConcentrations + minX;
    estimatedConcentrationsStdScaled = (maxX - minX) .* estimatedConcentrationsStd + minX;
    output.meas_raw_x_vals = concDataScaled ;
    output.meas_est_x_vals = estimatedConcentrationsScaled;
    output.meas_est_x_vals_sd = estimatedConcentrationsStdScaled;
    output.indicesDays = indicesDays ;
    output.measuredMetabolites = measuredMetabolites;
    output.time = time; 
    % Rxn_fluxes
    % Compute time-dependent parameters gamma and kappa
    if length(tn) > 1
        [gamma_t, kappa_t] = dynpars(timeInterpolated, tn);
    else
        gamma_t = timeInterpolated - timeInterpolated(1);
        kappa_t = ones(size(timeInterpolated));
    end
    Nr = size(reactionIDs, 1);
    % Normalize rates to the first reaction in the network
    tmpRxnFluxes = fit(nDMFA).V * kappa_t;
    tmpRxnFluxesStd = fit(nDMFA).V_sd * kappa_t;
    tmpRxnFluxesNorm = tmpRxnFluxes;
    tmpRxnFluxesStdNorm = tmpRxnFluxesStd;
    VPred = tmpRxnFluxesNorm(1:Nr, :);
    VPredNorm = VPred ./ VPred(1, :);
    VPredStd = tmpRxnFluxesStdNorm(1:Nr, :);
    VPredStd(abs(VPredStd) >= abs(VPred)) = 0.25 * abs(VPred(abs(VPredStd) >= abs(VPred)));
    VPredNormStd = VPredStd;
    VPredFitted = VPredNorm;
    VPredFittedStd = VPredNormStd;
    % Getting output variables
    output.flux_data_n = VPredNorm(:, indicesDays);
    output.flux_data_sd_n = VPredNormStd(:, indicesDays);
    output.flux_data = VPred(:, indicesDays);
    output.flux_data_sd = VPredStd(:, indicesDays);
    % End the timing and calculate the runtime
    output.EndRunTime = toc(startTime) / 60;
    % Store relevant data in the output structure
    output.model = model;
    output.simdata = simData;
    output.nDMFA = nDMFA;
    output.K = K;
    output.c0 = c0;
    output.U = U;
    output.tn = tn;
    output.tm = timeInterpolated;
    output.data = data;
    output.fit = fit;
    output.fit_ssr = fit(nDMFA).ssr;
    output.selection = selection ;
end

function [nDMFA, figNo] = findBestFit(fit, tol, robustFit)
    % Extract the necessary data from the fit structure
    SSRVec = [fit.ssr]';
    if robustFit
        SSRVec(SSRVec < tol) = NaN;
    end
    dmfaNT = [fit.nt]';
    SSRVecLB = [fit.ssr_lb95]';
    SSRVecUB = [fit.ssr_ub95]';
    % Get the SSR values within the given bounds
    SSRVecRange = SSRVec(SSRVec <= SSRVecUB & SSRVec >= SSRVecLB);
    if robustFit
        if ~isempty(SSRVecRange)
            % Sort the SSR values in descending order and get their indices
            [~, idx] = sort(SSRVecRange, 'descend');
            % Get the index corresponding to the lower quartile
            idx = floor(length(idx));
            idx = max(1, idx);  % Ensure the index is at least 1
            % Get the SSR value corresponding to this index
            SSRVal = SSRVecRange(idx);
            % Get the first index of this SSR value in the original SSR vector
            nDMFA = find(SSRVec == SSRVal, 1);
        else
            fprintf('\n Selecting the nearest to the SSR_ub...\n')
            % Get the index of the SSR value closest to the upper bound
            [~, idx] = min(abs(SSRVec - SSRVecUB));
            nDMFA = idx;
        end
    else
        % Find the index of the first converging value in SSRVec
        nDMFA = convergenceCheck(SSRVec, tol);
    end
    % Plot SSR vs. number of DMFA time points
    figNo = 0;
%     figure(figNo + 1);
%     plot(dmfaNT, SSRVec, '.-b');
%     hold on;
%     patch([dmfaNT' fliplr(dmfaNT')], [SSRVecUB' fliplr(SSRVecLB')], 'c', 'FaceAlpha', 0.05);
%     legend({'SSR', 'SSR-CI_{CI = (2.5 - 97.5)% }'});
%     legend('Box', 'off');
%     grid on;
%     xlabel('Maximum no. DMFA time points, nt');
%     ylabel('SSR');
%     hold on;
%     scatter(nDMFA, SSRVec(nDMFA), 'r*');
%     hold off;
%     figNo = figNo + 1;
end

function idx = convergenceCheck(SSRVec, tol)
    % Check for convergence in SSRVec
    n = length(SSRVec);
    converged = false(1, n);
    for i = 2:n
        if abs(SSRVec(i) - SSRVec(i-1)) <= tol
            converged(i) = true;
        end
    end
    if ~any(converged)
        idx = n;
    else
        idx = find(diff(abs(diff(SSRVec))<=tol)==1, 1, 'last')+1;
    end
end
% DMFA: Dynamic Metabolic Flux Analysis
% Performs dynamic metabolic flux analysis (DMFA) with given time points
% Dynamic metabolic flux analysis (DMFA) function
function fit = dmfa(model, data, tm, lambda, blockedRxnsIdx)
    % Dynamic metabolic flux analysis (DMFA).
    % Assign initial DMFA time points
    nm = length(tm);        % number of measurement time points
    tn(1) = tm(1);          % first DMFA time point (fixed)
    tn(nm-1) = tm(end);     % last DMFA time point (fixed)
    % Inflection time points (default)
    if nm > 3
        tn(2:nm-2) = (tm(2:end-2) + tm(3:end-1)) / 2;
    end
    % DMFA analysis with initial DMFA time points
    fit0 = dmfa_lslin(model, data, tm, tn, lambda, blockedRxnsIdx);   % DMFA analysis
    nt = length(tn);    % number of DMFA time points
    fit(nt) = fit0;     % store DMFA results
    % Iterative DMFA, determine optimal number of DMFA time points
    while nt > 1  % continue to reduce DMFA time points
        % Remove each inflection point, one-at-a-time, and perform DMFA analysis
        ssr = inf;   % initialize SSR
        for i = 2:nt-1 % loop over inflection time points
            tn = fit(nt).tn;      % restore DMFA time points
            tn(i) = [];           % remove inflection time point
            fit1 = dmfa_lslin(model, data, tm, tn, lambda, blockedRxnsIdx) ;  % DMFA analysis
            if fit1.ssr < ssr    % compare SSR to previous best fit
                ssr = fit1.ssr;    % record improved SSR value
                fit(nt-1) = fit1;   % store DMFA results
            end
        end
        % Update number of DMFA time points
        nt = nt - 1;
    end
    % Steady-state MFA analysis
    fit(1) = dmfa_lslin(model, data, tm, tm(1), lambda, blockedRxnsIdx);  % store MFA results
end
% Dynamic metabolic flux analysis with fixed DMFA time points
function fit = dmfa_lslin(model, data, tm, tn, lambda, blockedRxnsIdx)
   % Dynamic metabolic flux analysis (DMFA), with fixed DMFA time points.
    % Metabolic model
    nbal = strcmp(model.type, 'nbal'); % indices of non-balanced metabolites
    bal = strcmp(model.type, 'bal');   % indices of balanced metabolites
    Snbal = model.S(nbal, :); % stoichiometry matrix, non-balanced metabolites
    Sbal = model.S(bal, :);    % stoichiometry matrix, balanced metabolites
    K = null(Sbal, 'r');       % Kernel, null space of S for balanced metabolites
    % Size of variables
    nt = length(tn);        % number of DMFA time points
    nc = sum(nbal);         % number of non-balanced metabolites
    [nv, nu] = size(K);      % number of fluxes and free fluxes
    np = nc + nu * nt;        % number of parameters to be estimated
    % Compute time-dependent parameters gamma and kappa
    if nt > 1  % DMFA analysis
        [g, k] = dynpars(tm, tn);
    else  % MFA analysis
        g = tm - tm(1);
        k = ones(size(tm));
    end
    % Initialize matrix H and vector J
    H = zeros(np, np);        % init matrix H (Hessian)
    J = zeros(np, 1);         % init vector J (Jacobian)
    % Update matrix H and vector J with measurements
    for i = 1:length(data) % loop over all measurements
        % Find non-balanced metabolite in the model
        j = strcmp(data(i).met, model.met(nbal));   % index of nbal-metabolite
        % Collect measurement data
        sel = data(i).sel;        % selected measurement time points (true/false)
        nsel = sum(sel);          % number of selected measurements
        val = data(i).val(sel)';  % measured values
        W = diag(data(i).sd(sel).^(-2));   % weighting matrix
        % Update H and J
        switch data(i).type
            case 'metabolite'
                D = sparse(repmat(1:nt, nu, 1), 1:nu*nt, repmat(Snbal(j, :) * K, 1, nt));
                E = sparse(1:nsel, find(j), ones(1, nsel), nsel, nc);
                dcdp = [E, g(:, sel)' * D];
                H = H + dcdp' * W * dcdp;
                J = J + dcdp' * W * val;
            case 'rate'
                D = sparse(repmat(1:nt, nu, 1), 1:nu*nt, repmat(Snbal(j, :) * K, 1, nt));
                E = sparse(nsel, nc);
                drdp = [E, k(:, sel)' * D];
                H = H + drdp' * W * drdp;
                J = J + drdp' * W * val;
        end
    end
    % Check if every non-balanced metabolite has at least one conc measurement
    % If not, then c0 is non-observable and should be fixed (at 0 by default)
    nbalmet = model.met(nbal);
    for i = 1:length(nbalmet)
        % Find measurements for this non-balanced metabolite
        z = strcmp(nbalmet{i}, {data.met});
        % Check if measurements include at least one concentration measurement
        j = strcmp('metabolite', {data(z).type});
        % Fix metabolite initial concentration
        if ~any(j)
            H(i, i) = H(i, i) + 1;
        end
    end
    % Add regularization term to the Hessian
    H = H + lambda * eye(size(H));
    % Estimate DMFA parameters: initial concentrations and free fluxes
    %jCONSTRAINED OPTIMIZATION using baron to obtain global solution
    %% Constraints
    if isempty(blockedRxnsIdx)
        constrainedRxnIdx =  [1] ;
        Rxnbounds_lb = [0] ;
        Rxnbounds_ub = [Inf] ;
    else
        constrainedRxnIdx =  [1 blockedRxnsIdx] ;
        Rxnbounds_lb = [0 zeros(1, numel(blockedRxnsIdx))] ;
        Rxnbounds_ub = [inf zeros(1, numel(blockedRxnsIdx))] ;
    end
    baron_opts = baronset('LPSol', -1, 'LPAlg', 2, 'NLPSol', -1, 'PrLevel', 0, 'WantDual', 1 , 'threads', 30,...
                            'CplexLibName', 'libcplex2010.so', 'OutGrid', 100000, ...
                            'BoxTol', 1e-12, 'AbsIntFeasTol', 1e-12, 'AbsConFeasTol', 1e-12, ...
                            'RelIntFeasTol', 1e-12, 'RelConFeasTol', 1e-12, ...
                            'EpsA', 1e-12, 'EpsR', 1e-12) ;
    objM = @(x) (x'*H*x)/2 - J'*x; % quadratic term
    if isempty(Rxnbounds_lb)
        nlcon_cl = [] ; 
        nlcon = [] ;
    else 
        nlcon_cl = nlcon_bounds(Rxnbounds_lb, nt) ;
        nlcon = @(x)  nlcon_fun(x, constrainedRxnIdx , K, nc, nu, nt) ;
    end
    if isempty(Rxnbounds_ub)
        nlcon_cu = [] ;
    else
        nlcon_cu = nlcon_bounds(Rxnbounds_ub, nt) ;
    end
    % Call BARON to solve the problem
    [p,fval,exitflag,output] = baron(objM, [], [], [], [], [], nlcon, nlcon_cl, nlcon_cu, [] , NaN(size(H, 2), 1), baron_opts) ; %,  [], [], [], [], options);
function col_vector = nlcon_fun(x, j_vec, K, nc, nu, nt)
        x_matrix = reshape(x(nc + 1:end), nu, nt); 
        for idx = 1:length(j_vec)
            if idx== 1
                col_vector  = 	reshape(K(j_vec(idx), :) * x_matrix, nt, 1) ; 
            else 
                col_vector  = 	vertcat(col_vector ,  reshape(K(j_vec(idx), :) * x_matrix, nt, 1));
            end
        end
end
    function bound_vector = nlcon_bounds(bounds, nt)
    % This function returns a column vector based on the vector of 'j' indices and input 'bounds'
        % Calculate the column vector for each 'j' index and store in the array
        for idx = 1:length(bounds)
            % Initialize an array to store the resulting column vectors for each 'j' index
            if idx== 1
                bound_vector  = 	bounds(idx).*ones(nt, 1)  ; 
            else 
                bound_vector  = 	vertcat(bound_vector ,  bounds(idx).*ones(nt, 1) );
            end
        end
   end
  %
    fprintf('\n Optimization Exitflag is %d',  exitflag)
    pcov = inv(H);
    % Extract initial concentrations and free fluxes
    c0  = p(1:nc);
    U   = reshape(p(nc+1:end), nu, nt);
    % Return estimated model parameters and standard deviations
    fit.tn  = tn;     % DMFA time points
    fit.K   = K;      % null space of stoichiometry matrix
    fit.c0  = c0;     % estimated initial concentrations
    fit.U   = U;      % estimated free fluxes
    fit.V   = K * U;  % estimated fluxes
    fit.c0_sd = sqrt(diag(pcov(1:nc, 1:nc)));    % stdev of initial conc.
    fit.U_sd  = reshape(sqrt(diag(pcov(nc+1:end, nc+1:end))), nu, nt);
    fit.V_sd  = zeros(size(fit.V)); % init, stdev of fluxes
    for i = 1:nt  % loop over DMFA time points
        j = nc + nu * (i-1) + 1:nc + nu * i;  % indices for free fluxes
        fit.V_sd(:, i) = sqrt(diag(K * pcov(j, j) * K'));   % stdev of fluxes
    end
    % Calculate variance-weighted sum of squared residuals (SSR)
    nn = 0;    % init, number of fitted measurements
    ssr = 0;   % init, SSR
    for i = 1:length(data)  % loop over all measurements
        % Find non-balanced metabolite in the model
        j = strcmp(data(i).met, model.met(nbal));   % index of nbal-metabolite
        sel = data(i).sel;        % selected measurement time points (true/false)
        % Simulate measurements
        switch data(i).type
            case 'metabolite'
                simdata = c0(j) + Snbal(j, :) * K * U * g(:, sel);
            case 'rate'
                simdata = Snbal(j, :) * K * U * k(:, sel);
        end
        % Update SSR
        res = data(i).val(sel) - simdata;   % residual
        sd  = data(i).sd(sel);              % measurement stdev
        ssr = ssr + (res./sd) * (res./sd)';   % update SSR
        % Update number of fitted measurements
        nn = nn + sum(sel);
    end
    % Return goodness-of-fit analysis
    fit.ssr = ssr;         % variance-weighted sum of squared residuals SSR
    fit.nt  = nt;          % number of DMFA time points
    fit.n = nn;            % number of fitted measurements
    fit.p = np;            % number of estimated parameters
    fit.dof = nn - np;     % degrees of freedom
    % Return expected range of SSR values at 95% confidence level
    fit.ssr_lb95 = chi2inv(0.025, nn - np);
    fit.ssr_ub95 = chi2inv(0.975, nn - np);
end
% Simulate dynamic metabolic flux analysis measurements
function data = dmfa_sim_3(model, data, tm, K, c0, U, tn)
    % Simulation of dynamic metabolic flux analysis measurements.
    % Metabolic model
    nbal = strcmp(model.type, 'nbal'); % indices of non-balanced metabolites
    Snbal = model.S(nbal, :); % stoichiometry matrix, non-balanced metabolites
    % Compute time-dependent parameters gamma and kappa
    if length(tn) > 1  % DMFA analysis
        [g, k] = dynpars(tm, tn);
    else  % MFA analysis
        g = tm - tm(1);
        k = ones(size(tm));
    end
    % Simulate measurements
    for i = 1:length(data)
        % Find non-balanced metabolite in the model
        j = strcmp(data(i).met, model.met(nbal));   % index of nbal-metabolite
        sel = data(i).sel;        % selected measurement time points (true/false)
        % Return simulated measurements
        switch data(i).type
            case 'metabolite'
                data(i).sim = c0(j) + Snbal(j, :) * K * U * g(:, sel);
            case 'rate'
                data(i).sim = Snbal(j, :) * K * U * k(:, sel);
        end
    end
end

% Calculate time-dependent parameters gamma, kappa, and lambda
function [g, k, l, dg, dk, dg2, dk2] = dynpars(t, tn)
    % DYNPARS Calculate time-dependent parameters gamma, kappa, and lambda.
    % These parameters are used for dynamic metabolic flux analysis (DMFA).
    %
    %   [G, K, L, DG, DK, DG2, DK2] = DYNPARS(T, TN) computes time-dependent
    %   parameters gamma (G), kappa (K), lambda (L), and the first and second
    %   order derivatives of G and K with respect to TN, at the requested time(s) T. 
    %   TN are the DMFA time points. T should be within the time domain
    %   TN(1)...TN(end), which includes N-2 inflection points:
    %   TN = (t_begin, t_infection_1, t_infection_2, ..., t_end). 
    % Number of requested time points (m), and number of DMFA time points (n)
    m = length(t);   % number of requested time points
    n = length(tn);  % number of DMFA time points
    % Initialize matrices for parameters
    g = zeros(n, m);
    k = zeros(n, m);
    l = zeros(n, m);
    dg = zeros(n, m, n);  % first order derivatives
    dk = zeros(n, m, n);
    dg2 = zeros(n, m, n);  % second order derivatives
    dk2 = zeros(n, m, n);
    % Compute parameter values
    for j = 1:m  % loop over requested time points
        for i = 1:n  % loop over DMFA time points
            [g(i, j), k(i, j), l(i, j), dg(i, j, :), dk(i, j, :), dg2(i, j, :), dk2(i, j, :)] = ...
                calcpars(t(j), tn, i);
        end
    end
end
% Calculate time-dependent parameters gamma, kappa, and lambda for a specific DMFA time point
function [g, k, l, dg, dk, dg2, dk2] = calcpars(t, tn, i)
    % Calculate time-dependent parameters g, k, l and also first and second
    % order derivatives of g and k with respect to DMFA time points (tn).
    % This function calculates the values corresponding to the ith DMFA time point
    % at the requested time point t.  
    % The total number of DMFA time points
    n = length(tn);
    % Initialize
    g = 0;
    k = 0;
    l = 0;
    dg = zeros(1, n);
    dk = zeros(1, n);
    dg2 = zeros(1, n);
    dk2 = zeros(1, n);
    % Calculate parameter values for time point t
    if i > 1 && t >= tn(1) && t < tn(i-1)  % 1st row of table 1
        g = 0;
        k = 0;
        l = 0;
        dg = zeros(1, n);
        dk = zeros(1, n);
        dg2 = zeros(1, n);
        dk2 = zeros(1, n);
    elseif i > 1 && t >= tn(i-1) && t < tn(i)  % 2nd row of table 1
        g = (t - tn(i-1))^2 / (tn(i) - tn(i-1)) / 2;
        k = (t - tn(i-1)) / (tn(i) - tn(i-1));
        l = 1 / (tn(i) - tn(i-1));
        dg(i-1) = (t - tn(i))^2 / (2 * (tn(i) - tn(i-1))^2) - 1/2;
        dg(i) = -(t - tn(i-1))^2 / (2 * (tn(i) - tn(i-1))^2);    
        dg2(i-1) = -(t - tn(i))^2 / (tn(i-1) - tn(i))^3;
        dg2(i) = -(t - tn(i-1))^2 / (tn(i-1) - tn(i))^3;
        dk(i-1) = (t - tn(i)) / (tn(i) - tn(i-1))^2;
        dk(i) = -(t - tn(i-1)) / (tn(i) - tn(i-1))^2;
        dk2(i-1) = 2 * (t - tn(i)) / (tn(i) - tn(i-1))^3;
        dk2(i) = 2 * (tn(i-1) - t) / (tn(i-1) - tn(i))^3;
    elseif i < n && t >= tn(i) && t < tn(i+1)  % 3rd row of table 1
        g = (t - tn(i)) - (t - tn(i))^2 / (tn(i+1) - tn(i)) / 2 + (tn(i) - tn(max(1, i-1))) / 2;
        k = 1 - (t - tn(i)) / (tn(i+1) - tn(i));
        l = -1 / (tn(i+1) - tn(i));
        dg(i) = -(t - tn(i+1))^2 / (2 * (tn(i+1) - tn(i))^2);
        dg(i+1) = (t - tn(i))^2 / (2 * (tn(i+1) - tn(i))^2);
        dg(max(1, i-1)) = dg(max(1, i-1)) - 1/2;
        dg2(i) = -(t - tn(i+1))^2 / (tn(i+1) - tn(i))^3;
        dg2(i+1) = -(t - tn(i))^2 / (tn(i+1) - tn(i))^3;
        dk(i) = -(t - tn(i+1)) / (tn(i+1) - tn(i))^2;
        dk(i+1) = (t - tn(i)) / (tn(i+1) - tn(i))^2;
        dk2(i) = 2 * (tn(i+1) - t) / (tn(i+1) - tn(i))^3;
        dk2(i+1) = -2 * (t - tn(i)) / (tn(i+1) - tn(i))^3;
    elseif i < n && t >= tn(i+1)  % 4th row of table 1
        g = (tn(i+1) - tn(max(1, i-1))) / 2;
        k = 0;
        l = 0;
        dg(i+1) = 1/2;
        dg(max(1, i-1)) = dg(max(1, i-1)) - 1/2;
        if i+1 == n
            dk(n) = 1 / (tn(n) - tn(n-1));
            dk2(n) = -1 / (tn(n) - tn(n-1))^2;
        end
    elseif t == tn(n)  % 5th row of table 1
        g = (tn(n) - tn(n-1)) / 2;
        k = 1;
        l = 1 / (tn(n) - tn(n-1));
        dg(n) = -1/2;
        dg(n-1) = -1/2;
        dk(n) = -1 / (tn(n) - tn(n-1));
        dk2(n) = 1 / (tn(n) - tn(n-1))^2;
    else  % time point outside time domain
        errmsg = 'Requested time point is outside time domain.';
        disp(sprintf('%s: %s', mfilename, errmsg)); % display error message
    end
end
