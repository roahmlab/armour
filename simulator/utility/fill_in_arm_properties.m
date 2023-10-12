function  obj = fill_in_arm_properties(obj,arm,bound_state_limits)
% obj = fill_in_arm_properties(obj,arm)
% obj = fill_in_arm_properties(obj,arm,bound_state_limits)
%
% Given the object obj, iterate through its properties and fill in any that
% start with 'arm' according to the corresponding properties of the arm.
%
% If the third (optional) input is set to true, then any of the arm's joint
% state limits that are +Inf or -Inf are set to +pi and -pi, respectively,
% in the obj.arm_joint_state_limits property.
    
    % whether or not to bound state limits
    if nargin < 3
        bound_state_limits = false ;
    end

    % get the object's fieldnames
    F = fieldnames(obj) ;
    
    % iterate through the fieldnames and fill in the arm properties
    for idx = 1:length(F)
        f_obj = F{idx} ;
        if length(f_obj) > 4 && strcmp(f_obj(1:4),'arm_')
            f_arm = f_obj(5:end) ;
            
            if strcmp(f_arm,'joint_state_limits') && bound_state_limits
                L = arm.(f_arm) ;
                inf_log = isinf(L) ;
                L(1,inf_log(1,:)) = -pi;
                L(2,inf_log(2,:)) = +pi ;
                obj.(f_obj) = L ;
            else
                obj.(f_obj) = arm.(f_arm) ;
            end
        end
    end
end