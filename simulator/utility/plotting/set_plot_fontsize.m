function set_plot_fontsize(fontsize,fh)
% set_plot_fontsize(fontsize,fh)
%
% Set all the font sizes in the current figure, or else in the provided
% figure handle. The default font size is 15.
%
% Method is taken from here: https://www.mathworks.com/matlabcentral/answers/223344-changing-font-size-in-all-the-elements-of-figures
%
% Authors: Shreyas Kousik
% Created: 23 Apr 2021
% Updated: noap

    if nargin < 1
        fontsize = 15 ;
    end
    
    if nargin < 2
        fh = gcf ;
    end
    
    % send it
    set(findall(fh,'-property','FontSize'),'FontSize',fontsize)
end