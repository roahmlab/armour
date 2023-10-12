function out = check_if_plot_is_available(simulator_object,fieldname)
    % out = check_if_plot_is_available(simulator_object,fieldname)
    %
    % Given a simulator object (agent, world, or planner), and a fieldname
    % of object.plot_data, check if that fieldname (and the corresponding
    % plot data) is available. Also check if a figure is up in the global
    % figure root. If both the figure and the object's plot data are
    % available, return true.
    %
    % If simulator_object.plot_data.(fieldname) contains a cell array of
    % n plotted object handles, this function checks each object in the
    % cell array.
    %
    % Author: Shreyas Kousik
    % Created: shrug
    % Updated: 19 Mar 2020
    
    % if no fieldname is passed in, check the first fieldname in the
    % object's plot_data property
    if nargin < 2
        F = fieldnames(simulator_object.plot_data) ;
        fieldname = F{1} ;
    end
    
    try
        % try to get the data in the given field (if the field doesn't
        % exist, this will throw an error)
        h = simulator_object.plot_data.(fieldname) ;
        
        % get the current figure (if there is no figure up, then this will
        % be an empty GraphicsPlaceholder array)
        fh = get(groot,'CurrentFigure') ;
        
        % check if plot is available:
        if iscell(h)
            % if h is a cell array of plotted objects, check each one
            out = false(1,length(h)) ;
            
            if ~isempty(fh) % if the figure is up...
                for idx = 1:length(h)
                    % iterate through the plotted objects in h, and check
                    % that each one is nonempty and a valid plot object
                    out(idx) = ~(isempty(h{idx}) || ~all(isvalid(h{idx}))) ;
                end
            end
        else
            % otherwise, just check if h is nonempty and valid (meaning
            % that the object has been plotted correctly) and if the figure
            % is up (meaning fh is nonempty)
            out = ~(isempty(h) || ~all(isvalid(h)) || isempty(fh)) ;
        end
    catch
        % the only way anything above should throw an error is if h is an
        % invalid fieldname, hence the following warning:
        warning([fieldname,' may not be not a valid plot_data field of the ',...
            simulator_object.name, ' object'])
        out = false ;
    end
end