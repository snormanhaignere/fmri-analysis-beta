function glm_across_subjects_surf_grid(exp, usubs, runtype, ...
    fwhm, analysis_name, grid_spacing_mm, grid_roi)


%%

% expand experiment and runtype to cell array with one element per subject
if ischar(exp)
    exp = repmat({exp}, 1, length(usubs));
end
if ischar(runtype)
    runtype = repmat({runtype}, 1, length(usubs));
end







