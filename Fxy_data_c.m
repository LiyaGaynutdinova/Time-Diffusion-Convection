classdef Fxy_data_c
    % jen fce x a y
    properties
        a,b % not used now
    end
    methods
        function obj = Fxy_data_c(a,b)            
            obj.a = a;
            obj.b = b;
        end
        function z = Eval(obj,x1,x2)            
            z = [1,1];            
        end
    end
end

