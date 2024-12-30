classdef Fxy_data_a
    % jen fce x a y
    properties
        a,b  % not used now
    end
    methods
        function obj = Fxy_data_a(a,b)            
            obj.a = a;
            obj.b = b;
        end
        function z = Eval(obj,x1,x2)           
            z = [1,0;0,1]*(1.1-1*sin(x2)); 
            z = [1,0;0,1]*(1.1-1*sign(sin(x2))); 
        end
        function z = Fxy_ref_a(obj) 
            z = [1,0;0,1];
        end
    end
end