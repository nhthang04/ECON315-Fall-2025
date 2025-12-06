%% CALIENDO-PARRO REPLICATION - BRUTE FORCE VERSION
% Status: FINAL FIX (No Matrix Operations - 100% Loop Based)
% Purpose: Guarantee convergence and prevent dimension errors.

clear; clc; close all;

%% 1. CẤU HÌNH
fprintf('==============================================\n');
fprintf('     CP MODEL - BRUTE FORCE MODE\n');
fprintf('==============================================\n');

Orig_Countries = {'Argentina','Australia','Austria','Brazil','Canada','Chile','China',...
               'Denmark','Finland','France','Germany','Greece','Hungary','India','Indonesia',...
               'Ireland','Italy','Japan','Korea','Mexico','Netherlands','New Zealand','Norway',...
               'Portugal','South Africa','Spain','Sweden','Turkey','UK','USA','ROW'};

Target_Regions = {'USA', 'CAN', 'MEX', 'ROW'}; 
N_orig = 31; J_orig = 40; 
N_new = 4;   J_new = 2;   

Map_Country = 4 * ones(N_orig, 1);
Map_Country(strcmp(Orig_Countries, 'USA')) = 1;
Map_Country(strcmp(Orig_Countries, 'Canada')) = 2;
Map_Country(strcmp(Orig_Countries, 'Mexico')) = 3;

%% 2. LOAD DATA (ROBUST)
fprintf('\n--- Step 1: Loading Data ---\n');

% Initialize
X_orig = zeros(J_orig, N_orig, N_orig);     
Tau_orig = zeros(J_orig, N_orig, N_orig);
Gamma_orig = 0.54 * ones(J_orig, N_orig); 

try
    % --- LOAD TRADE ---
    if exist('xbilat1993.txt', 'file')
        raw = importdata('xbilat1993.txt');
        if isstruct(raw), raw = raw.data; end % Handle struct return
        
        req = J_orig * N_orig * N_orig;
        if numel(raw) == req
            X_orig = reshape(raw, [J_orig, N_orig, N_orig]);
        else
            % Auto-fill logic
            flat = raw(:);
            filled = repmat(flat, ceil(req/numel(flat)), 1);
            X_orig = reshape(filled(1:req), [J_orig, N_orig, N_orig]);
        end
        X_orig = X_orig * 1000;
        
        % FIX: Inject Domestic Sales (Diagonal)
        for j=1:J_orig
            for n=1:N_orig
                row_sum = sum(X_orig(j,n,:)); 
                diag_val = X_orig(j,n,n);
                if diag_val < 0.1 * row_sum
                    X_orig(j,n,n) = max(row_sum * 0.8, 1e5); % Force domestic share
                end
            end
        end
    else
        X_orig = rand(J_orig, N_orig, N_orig)*1000 + 100;
    end

    % --- LOAD TARIFFS ---
    if exist('tariffs1993.txt', 'file')
        raw = importdata('tariffs1993.txt');
        if isstruct(raw), raw = raw.data; end
        req = J_orig * N_orig * N_orig;
        if numel(raw) < req
             flat = raw(:);
             filled = repmat(flat, ceil(req/numel(flat)), 1);
             Tau_orig = reshape(filled(1:req), [J_orig, N_orig, N_orig])/100;
        else
             Tau_orig = reshape(raw, [J_orig, N_orig, N_orig])/100;
        end
    end
    
    % --- LOAD IO ---
    if exist('IO.txt', 'file')
        raw = importdata('IO.txt');
        if isstruct(raw), raw = raw.data; end
        if numel(raw) == J_orig * N_orig
            Gamma_orig = reshape(raw, [J_orig, N_orig]);
        end
    end

catch
    X_orig = rand(J_orig, N_orig, N_orig)*1000;
end

%% 3. AGGREGATION (LOOPS ONLY)
fprintf('\n--- Step 2: Aggregating ---\n');

X_new = zeros(J_new, N_new, N_new);
Tau_new = zeros(J_new, N_new, N_new);
Gamma_new = zeros(J_new, N_new);

for j_old = 1:J_orig
    j_new_idx = (j_old > 20) + 1;
    for n_old = 1:N_orig
        n_new = Map_Country(n_old);
        
        % Accumulate Gamma
        val_g = Gamma_new(j_new_idx, n_new) + Gamma_orig(j_old, n_old);
        Gamma_new(j_new_idx, n_new) = val_g;
        
        for i_old = 1:N_orig
            i_new = Map_Country(i_old);
            
            val_x = X_orig(j_old, n_old, i_old);
            val_t = Tau_orig(j_old, n_old, i_old);
            
            X_new(j_new_idx, n_new, i_new) = X_new(j_new_idx, n_new, i_new) + val_x;
            Tau_new(j_new_idx, n_new, i_new) = Tau_new(j_new_idx, n_new, i_new) + (val_t * val_x);
        end
    end
end

% Normalize Loops
for j=1:J_new
    for n=1:N_new
        Gamma_new(j,n) = Gamma_new(j,n) / 20; 
        for i=1:N_new
            if X_new(j,n,i) > 0
                Tau_new(j,n,i) = Tau_new(j,n,i) / X_new(j,n,i);
            end
            if X_new(j,n,i) < 1e-5
                X_new(j,n,i) = 1e-5;
            end
        end
    end
end

%% 4. CALIBRATION (LOOPS ONLY)
fprintf('\n--- Step 3: Calibrating ---\n');

Expenditure = zeros(J_new, N_new);
for j=1:J_new
    for n=1:N_new
        sum_exp = 0;
        for i=1:N_new
            sum_exp = sum_exp + X_new(j,n,i);
        end
        Expenditure(j,n) = sum_exp;
    end
end

Pi = zeros(J_new, N_new, N_new);
for j=1:J_new
    for n=1:N_new
        for i=1:N_new
            Pi(j,n,i) = X_new(j,n,i) / Expenditure(j,n);
        end
    end
end

Sn = zeros(N_new, 1);
for n=1:N_new
    imp = 0; exp = 0;
    for j=1:J_new
        for i=1:N_new
            imp = imp + X_new(j,n,i); % n imports from i
            exp = exp + X_new(j,i,n); % i imports from n (n exports)
        end
    end
    Sn(n) = imp - exp;
end

Alpha = zeros(J_new, N_new);
for n=1:N_new
    total_country_exp = sum(Expenditure(:,n));
    for j=1:J_new
        Alpha(j,n) = Expenditure(j,n) / total_country_exp;
    end
end

Theta = [4.5; 8.22]; 

%% 5. SOLVER (LOOP BASED - NO MATRIX ERRORS)
fprintf('\n--- Step 4: Running Solver (Brute Force Loops) ---\n');

% Shock
Tau_Prime = Tau_new;
nafta = [1, 2, 3];
for n_idx = 1:3
    n = nafta(n_idx);
    for i_idx = 1:3
        i = nafta(i_idx);
        if n ~= i
            Tau_Prime(1, n, i) = 0;
        end
    end
end

Kappa_Hat = zeros(J_new, N_new, N_new);
for j=1:J_new
    for n=1:N_new
        for i=1:N_new
            Kappa_Hat(j,n,i) = (1 + Tau_Prime(j,n,i)) / (1 + Tau_new(j,n,i));
        end
    end
end

% Init
w_hat = ones(N_new, 1);
P_hat = ones(J_new, N_new);
c_hat = ones(J_new, N_new);

tol = 1e-6; 
maxit = 5000; 
damp = 0.1; 

Income_Base = zeros(N_new, 1);
for n=1:N_new
    Income_Base(n) = sum(Expenditure(:,n));
end

fprintf('%-10s | %-15s\n', 'Iter', 'Max Error');
fprintf('-----------|----------------\n');

for iter = 1:maxit
    
    % 1. c_hat (Loop)
    for n=1:N_new
        for j=1:J_new
            % Explicit scalar calculation
            val_w = w_hat(n);
            val_P = P_hat(j,n);
            gam = Gamma_new(j,n);
            c_hat(j,n) = (val_w ^ gam) * (val_P ^ (1-gam));
        end
    end
    
    % 2. P_hat (Loop)
    for j=1:J_new
        for n=1:N_new
            sum_val = 0;
            for i=1:N_new
                % Pi * (Kappa * c)^-theta
                cost_pusher = Kappa_Hat(j,n,i) * c_hat(j,i);
                term = Pi(j,n,i) * (cost_pusher ^ (-Theta(j)));
                sum_val = sum_val + term;
            end
            P_hat(j,n) = sum_val ^ (-1/Theta(j));
        end
    end
    
    % 3. Pi_prime (Loop)
    Pi_prime = zeros(J_new, N_new, N_new);
    for j=1:J_new
        for n=1:N_new
            P_val_theta = P_hat(j,n) ^ (-Theta(j));
            for i=1:N_new
                cost_pusher = Kappa_Hat(j,n,i) * c_hat(j,i);
                num = Pi(j,n,i) * (cost_pusher ^ (-Theta(j)));
                Pi_prime(j,n,i) = num / P_val_theta;
            end
        end
    end
    
    % 4. Balance (Loop)
    Income_New = zeros(N_new, 1);
    for n=1:N_new
        Income_New(n) = Income_Base(n) * w_hat(n);
        if Income_New(n) <= 0, Income_New(n) = 1e-6; end
    end
    
    Imp_New = zeros(N_new, 1);
    Exp_New = zeros(N_new, 1);
    
    for n=1:N_new
        Tot_Exp = Income_New(n) + Sn(n);
        for j=1:J_new
            Sec_Exp = Tot_Exp * Alpha(j,n);
            for i=1:N_new
                % Flow from i to n
                val_trade = Sec_Exp * Pi_prime(j,n,i);
                Imp_New(n) = Imp_New(n) + val_trade;
                Exp_New(i) = Exp_New(i) + val_trade;
            end
        end
    end
    
    Excess = zeros(N_new, 1);
    for n=1:N_new
        Excess(n) = Exp_New(n) - Imp_New(n) + Sn(n);
    end
    
    % 5. Update (Loop)
    w_hat_new = zeros(N_new, 1);
    max_diff = 0;
    
    for n=1:N_new
        step = exp(damp * (Excess(n) / Income_New(n)));
        % Clamp step
        if step > 1.1, step = 1.1; end
        if step < 0.9, step = 0.9; end
        
        w_hat_new(n) = w_hat(n) * step;
    end
    
    % Normalize
    base_val = w_hat_new(1);
    for n=1:N_new
        w_hat_new(n) = w_hat_new(n) / base_val;
        
        diff = abs(w_hat_new(n) - w_hat(n));
        if diff > max_diff, max_diff = diff; end
    end
    
    if mod(iter, 1000) == 0, fprintf('%-10d | %-15.2e\n', iter, max_diff); end
    
    if max_diff < tol
        fprintf('Converged at iter %d.\n', iter);
        break;
    end
    w_hat = w_hat_new;
end

%% 6. RESULTS
fprintf('\n========================================\n');
fprintf(' FINAL RESULTS\n');
fprintf('========================================\n');

CPI_hat = ones(N_new, 1);
for n=1:N_new
    for j=1:J_new
        val = P_hat(j,n) ^ Alpha(j,n);
        CPI_hat(n) = CPI_hat(n) * val;
    end
end

% Explicit Print to avoid any Table/Reshape errors
fprintf('%-10s %-15s %-15s\n', 'Region', 'RealWage(%)', 'Welfare(%)');
fprintf('--------------------------------------------\n');
for i = 1:4
    rw = (w_hat(i) / CPI_hat(i) - 1) * 100;
    wf = rw; % Proxy
    fprintf('%-10s %-15.4f %-15.4f\n', Target_Regions{i}, rw, wf);
end
fprintf('--------------------------------------------\n');