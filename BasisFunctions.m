classdef BasisFunctions
    enumeration
        Delta
    end

    % Evaulaute the inner product of the function and constant voltage 1
    methods
        function n = magnitude(obj)
            n = 0;
            switch obj
                case Delta
                    n = 1;
            end
        end
    end

    
end