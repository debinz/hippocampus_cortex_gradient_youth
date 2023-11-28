function h = scree_plot(lambdas)
% SCREE_PLOT   Plot of the scaled eigenvalues. 
%
%   h = SCREE_PLOT(lambdas) plots the eigenvalues corresponding to a
%   gradient. Lambdas are scaled to a sum of 1. All graphics objects
%   created are stored in structure array h. 
%
%   For more information, please consult our <a
%   href="https://brainspace.readthedocs.io/en/latest/pages/matlab_doc/visualization/scree_plot.html">ReadTheDocs</a>.


h.figure = figure('Color','White', 'Position', [100 100 700 700]);
h.axes = axes(); 
h.plot = plot(lambdas ./ sum(lambdas),'o-','Color','k');
xlabel('Component Number');
ylabel('Explained variance'); %'Scaled Eigenvalues'
set(h.axes,'box','off','FontName','DroidSans','FontSize',14)

