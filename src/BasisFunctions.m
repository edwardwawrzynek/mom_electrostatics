classdef BasisFunctions
    enumeration
        Delta, Pulse, Triangle
    end

    % Evaulaute the inner product of the basis function and constant voltage 1
    methods
        function n = magnitude(obj, pt, ptl, ptr)
            n = 0;
            switch obj
                case BasisFunctions.Delta
                    n = 1;
                case BasisFunctions.Pulse
                    % norm is the length of the line segment
                    n = 1;
                case BasisFunctions.Triangle
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
                case BasisFunctions.Pulse
                    % Test at the midpoint of the segments
                    midpoint = 0.5 * (w_pt + w_ptl);
                    % Evaluate voltage from charge on sigma_ptl - sigma_pt
                    v = obj.evaluateVoltage(midpoint, sigma_pt, sigma_ptl, sigma_ptr, 0.5 * norm(w_pt - w_ptl));
                case BasisFunctions.Triangle
                    % Test at the midpoint of the segments
                    midpoint = 0.5 * (w_pt + w_ptl);
                    % Evaluate voltage from charge on sigma_ptl - sigma_pt
                    v = obj.evaluateVoltage(midpoint, sigma_pt, sigma_ptl, sigma_ptr, 0.5 * norm(w_pt - w_ptl));
            end
        end
        
        % evaluate voltage of a triangle (1 at ptl, 0 at ptr)
        function v = evaluateHalfTriangleVoltage(obj, v_pt, ptl, ptr, min_dist)
            epsilon0 = 8.85418781e-12;
            axial_v_weighted = @(L,r,a,b) 1/(4*pi*epsilon0) * (a*log((L+sqrt(r^2+L^2))/(-L+sqrt(r^2+L^2))) + 2*b*r - 2*b*sqrt(r^2+L^2));

            % compute vectors from point to endpoints
            r1 = ptl - v_pt;
            r2 = ptr - v_pt;
            

            % compute normnalized vector along line
            r3 = r2 - r1;
            d = norm(r3);
            u3 = r3 / d;

            L1 = min(abs(dot(r1,u3)), abs(dot(r2, u3)));
            r = max(norm(r2 - dot(r2,u3)*u3), min_dist);

            case_sgn = dot(r1,u3)*dot(r2,u3);
            if case_sgn >= 0 % case 1
                L2 = L1 + d;

                a = L2 / (L2 - L1);
                b = 1/(L2-L1);
                
                % correct for orientation of triangle
                if dot(r2, u3) < 0
                    a = -a + 1;
                    b = -b;
                end

                v = 1/2 * (axial_v_weighted(L2,r,a,b) - axial_v_weighted(L1,r,a,b));
            else % case 2
                L2 = d - L1;
                
                a = L1 / (L2 + L1);
                b = 1 / (L2 + L1);

                % correct for orientation of triangle
                if norm(r1) < norm(r2)
                    a = L2 / (L2 + L1);
                    b = -b;
                end

                v = 1/2 * (axial_v_weighted(L1,r,a,b) + axial_v_weighted(L2,r,a,-b));
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
                case BasisFunctions.Pulse
                    axial_v = @(L,r) 1/(4*pi*epsilon0) * log((L+sqrt(r^2+L^2))/(-L+sqrt(r^2+L^2)));

                    % compute vectors from point to endpoints
                    r1 = sigma_ptl - v_pt;
                    r2 = sigma_pt - v_pt;
                    

                    % compute normnalized vector along line
                    r3 = r2 - r1;
                    d = norm(r3);
                    u3 = r3 / d;

                    L1 = min(abs(dot(r1,u3)), abs(dot(r2, u3)));
                    r = max(norm(r1 - dot(r1,u3)*u3), min_dist);

                    case_sgn = dot(r1,u3)*dot(r2,u3);
                    if case_sgn >= 0
                        L2 = L1 + d;

                        v = 1/2 * (axial_v(L2,r) - axial_v(L1,r));
                    else
                        L2 = d - L1;
                        v = 1/2 * (axial_v(L2,r) + axial_v(L1,r));
                    end
                case BasisFunctions.Triangle
                    v = obj.evaluateHalfTriangleVoltage(v_pt, sigma_pt, sigma_ptl, min_dist) + obj.evaluateHalfTriangleVoltage(v_pt, sigma_pt, sigma_ptr, min_dist);
            end

        end
    end
end