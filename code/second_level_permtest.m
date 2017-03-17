function second_level_permtest(...
    MAT_files_first_level, perm_MAT_files_first_level, ...
    perm_MAT_file_second_level)

% average across runs
n_runs = length(perm_MAT_files_first_level);
for i = 1:n_runs
    X = load(perm_MAT_files_first_level{i}, ...
        'beta_contrast_permtest', 'residual_permtest');
    Y = load(MAT_files_first_level{i}, 'beta_contrast', 'residual');
    if i == 1
        beta_contrast = Y.beta_contrast;
        counts_beta_contrast = double(~isnan(beta_contrast));
        residual = Y.residual;
        counts_residual = double(~isnan(residual));
        beta_contrast_permtest = X.beta_contrast_permtest;
        counts_beta_contrast_permtest = double(~isnan(beta_contrast_permtest));
        residual_permtest = X.residual_permtest;
        counts_residual_permtest = double(~isnan(residual_permtest));
    else
        residual = ...
            nanplus( residual, Y.residual );
        xi = ~isnan(Y.residual);
        counts_residual(xi) = counts_residual(xi)+1;
        clear xi;
        
        beta_contrast = ...
            nanplus( beta_contrast, Y.beta_contrast );
        xi = ~isnan(Y.beta_contrast);
        counts_beta_contrast(xi) = counts_beta_contrast(xi)+1;
        clear xi;
        
        beta_contrast_permtest = ...
            nanplus( beta_contrast_permtest, X.beta_contrast_permtest );
        xi = ~isnan(X.beta_contrast_permtest);
        counts_beta_contrast_permtest(xi) = counts_beta_contrast_permtest(xi)+1;
        clear xi;
        
        residual_permtest = ...
            nanplus(residual_permtest, X.residual_permtest);
        xi = ~isnan(X.residual_permtest);
        counts_residual_permtest(xi) = counts_residual_permtest(xi)+1;
        clear xi;
    end
end
residual = residual ./ counts_residual;
beta_contrast = beta_contrast ./ counts_beta_contrast;
beta_contrast_permtest = beta_contrast_permtest ./ counts_beta_contrast_permtest;
residual_permtest = residual_permtest ./ counts_residual_permtest;

% estimate P value
logP_permtest = sig_via_null_gaussfit(beta_contrast, beta_contrast_permtest); %#ok<NASGU>
logP_residual_permtest = sig_via_null_gaussfit(residual, residual_permtest, ...
    'tail', 'left'); %#ok<NASGU>

% save results
save(perm_MAT_file_second_level, 'logP_permtest', 'beta_contrast_permtest', ...
    'residual_permtest', 'logP_residual_permtest');

function Z = nanplus(X,Y)

assert(all(size(X)==size(Y)))
Z = nan(size(X));

xi = ~isnan(X) & ~isnan(Y);
Z(xi) = X(xi) + Y(xi);

xi = ~isnan(X) & isnan(Y);
Z(xi) = X(xi);

xi = isnan(X) & ~isnan(Y);
Z(xi) = Y(xi);