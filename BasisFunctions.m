classdef BasisFunctions
    enumeration
        Delta
    end

    % Evaulaute the inner product of the function and constant voltage 1
    methods
        function n = magnitude(obj, pt, ptl, ptr)
            n = 0;
            switch obj
                case BasisFunctions.Delta
                    n = 1;
            end
        end

        % Evaluate the inner product of the operator applied to the function
        % at sigma and the test function at w
        % 
        % sigma_pt, sigma_ptl, sigma_ptr are location of the charge
        % w_pt, w_ptl, w_ptr are the location of the test function
        function v = innerProduct(obj, w_pt, w_ptl, w_ptr, sigma_pt, sigma_ptl, sigma_ptr)
            epsilon0 = 8.85418781e-12;

            v = 0;
            switch obj
                case BasisFunctions.Delta
                    % For delta functions we test at the midpoint of the
                    % charge basis functions
                    midpoint = 0.5 * (w_pt + w_ptl);
                    % Evaluate voltage from charge at pt
                    v = 1/(4*pi*epsilon0) * 1 / norm(midpoint - sigma_pt);
            end
        end

        % evaluate the voltage at v_pt created by a basis function at sigma_ptr
        function v = evaluateVoltage(obj, v_pt, sigma_pt, sigma_ptl, sigma_ptr, min_dist)
            epsilon0 = 8.85418781e-12;
            
            v = 0;
            switch obj
                case BasisFunctions.Delta
                    dist = max(norm(v_pt - sigma_pt), min_dist);

                    v = 1/(4*pi*epsilon0) * 1 / dist;
            end

        end
    end
end