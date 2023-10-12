function set_plot_linewidths(linewidth,fh)
% set_plot_linewidths(linewidth,fh)
%
% Set linewidths of all lines and patches to the provided linewidth
% (default is 1.5) in the current figure, or in the provided figure handle
% fh.
%
% Author: Shreyas Kousik
% Created: 29 Oct 2019
% Updated: -

    if nargin < 1
        linewidth = 1.5 ;
    end
    
    if nargin < 2
        fh = gcf ;
    end

    % get all lines
    h = findall(fh,'Type','Line') ;
    for idx = 1:length(h)
        h(idx).LineWidth = linewidth ;
    end
    
    % get all patches
    h = findall(fh,'Type','Patch') ;
    for idx = 1:length(h)
        h(idx).LineWidth = linewidth ;
    end
end