classdef Fxy_RHS
    % jen fce x a y
    properties
        a,b  % not used now
    end
    methods
        function obj = Fxy_RHS(a,b)            
            obj.a = a;
            obj.b = b;
        end
        function z = Eval(obj,x1,x2)           
            z = sin(15.1*x1+19.2*x2)+cos(x1+x2)/2;% like almost random RHS
        end
    end
end