function [varargout] = match_trajectories(T_des,varargin)
% [Z_1,Z_2,...] = match_trajectories(T_des,T_1,Z_1,T_2,Z_2,...,interp_type)
%
% Given desired sample times T_des and any number of time vectors (1-by-n_t)
% and associated trajectories (n_states-by-n_t), reinterpolate the given
% trajectories linearly at the desired times.
%
% If T_des(1) < T_i(1), then the output Z_i is pre-padded with Z_i(:,1),
% and similarly if T_des(end) > T_i(end).
%
% Authors: Shreyas Kousik
% Created: long long ago, at a desktop far, far away
% Updated: 7 Jan 2021

    varargout = cell(1,nargout) ;
    
    if ischar(varargin{end})
        interp_type = varargin{end} ;
    else
        interp_type = 'linear' ;
    end
    
    out_idx = 1 ;
    for idx = 1:2:(length(varargin)-1)
        % get the current time and trajectory
        T = varargin{idx} ;
        Z = varargin{idx+1} ;
        
        % if T_des exceeds the bounds of T, pad T and Z accordingly
        if T_des(1) < T(1)
            T = [T_des(1), T] ;
            Z = [Z(:,1), Z] ;
        end
        
        if T_des(end) > T(end)
            T = [T, T_des(end)] ;
            Z = [Z, Z(:,end)] ;
        end
        
        if length(T) == 1 && T_des == T
            % check if no interpolation is needed (this resolves a bug
            % where interp1 reinterprets the length of Z as a vector to be
            % interpolated, and throws an error because T and Z are not of
            % the same length)
            varargout{out_idx} = Z;
        else
            % call interp1 as usual
            varargout{out_idx} = interp1(T(:),Z',T_des(:),interp_type)' ;
        end
        
        out_idx = out_idx + 1 ;
    end
end