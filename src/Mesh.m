classdef Mesh
    properties
        num_pts
        points  % array of ordered points in the mesh
        basis % the type of basis and testing functions to use
        weights % solved weight of charge distribution basis
    end
    methods
        function obj = Mesh(pts, basis)
            obj.points = pts;
            obj.num_pts = length(pts);
            obj.basis = basis;
            obj.weights = zeros(obj.num_pts,1);
        end
        
        % Display the mesh points and edges
        function plotMesh(obj)
            figure;
            plot([obj.points(:,1); obj.points(1,1)], [obj.points(:,2); obj.points(1,2)], '-o');
        end

        % Display the mesh colored by charge density
        function plotCharge(obj)
            % compute the absolue maximum charge value
            max_charge = max(obj.weights);
            min_charge = min(obj.weights);
            
            figure;
            title("Surface Charge Density [C]");
            cm = colormap(jet);
            colorbar;
            caxis([min_charge max_charge]);

            plot_pts = [[obj.points(:,1); obj.points(1,1)], [obj.points(:,2); obj.points(1,2)]];
            for i = 1:obj.num_pts
                color0 = cm(floor(((obj.weights(i)-min_charge)/(max_charge-min_charge)) * (size(cm, 1)-1)) + 1,:);
                
                line([plot_pts(i,1) plot_pts(i+1,1)], [plot_pts(i,2) plot_pts(i+1,2)], ...
                    'Color', color0, ...
                    'LineWidth', 2);
            end

            grid on;
            axis square;
        end

        % Compute the analytic voltage distribution for a hoop of constant
        % charge density sigma centered at the origin with radius 1
        function [v, xleft, xright, yleft, yright] = computeVoltageHoop(obj,scale_factor, sigma, n)
            epsilon0 = 8.85418781e-12;
            
            % Setup scale to mimic the numerically solved hoop
            xleft = -scale_factor;
            xright = scale_factor;
            yleft = -scale_factor;
            yright = scale_factor;
            xdist = 2;
            ydist = xdist;
            min_dist = 0.5 * min(xdist, ydist)*scale_factor/n;

            v = zeros(n);

            for i = 1:n
                for j = 1:n
                    pt = [xleft + (j-1) * scale_factor * xdist/n, yleft + (i-1) * scale_factor * ydist/n];
                    r = sqrt(pt(1)^2 + pt(2)^2);
                    % make sure we are at least min_dist away from the
                    % surface (to avoid voltage singularity on the surface)
                    if r <= 1 && r > 1-min_dist
                        r = 1-min_dist;
                    elseif r > 1 && r < 1 + min_dist
                        r = 1 + min_dist;
                    end

                    v(i,j) = sigma / (4*pi*epsilon0) * ...
                        integral(@(theta) 1./sqrt(r^2 + 1 - 2*r*cos(theta)), 0, 2*pi, "AbsTol", 1e-16);
                end
            end
        end


        % Compute voltage over space created by the charge distribution
        function [v, xleft, xright, yleft, yright] = computeVoltage(obj, scale_factor, n)
            % Get range over which to evaluate voltage
            x = obj.points(:,1);
            y = obj.points(:,2);

            xcenter = 0.5 * (max(x)+min(x));
            xdist = max(x)-min(x);
            ycenter = 0.5 * (max(y)+min(y));
            ydist = max(y)-min(y);

            xleft = xcenter - 0.5 * scale_factor * xdist;
            xright = xcenter + 0.5 * scale_factor * xdist;
            yleft = ycenter - 0.5 * scale_factor * ydist;
            yright = ycenter + 0.5 * scale_factor * ydist;

            % evaluate voltage
            pts_neighbors = [obj.points(obj.num_pts,:); obj.points; obj.points(1,:)];

            v = zeros(n);
            
            % minimum distance at which to evaluate v
            min_dist = 0.5 * min(xdist, ydist)*scale_factor/n;

            epsilon0 = 8.85418781e-12;
            
            for i = 1:n
                for j = 1:n
                    pt = [xleft + (j-1) * scale_factor * xdist/n, yleft + (i-1) * scale_factor * ydist/n];
                    % Sum contribution from each basis
                    for k = 1:obj.num_pts
                        v(i,j) = v(i,j) + obj.weights(k) * obj.basis.evaluateVoltage(pt, pts_neighbors(k+1,:), pts_neighbors(k,:), pts_neighbors(k+2,:), min_dist);
                        
                        %dist = max(norm(pt - obj.points(k,:)), min_dist);
                        %v(i,j) = v(i,j) + obj.weights(k) / (4*pi*epsilon0) * 1 / dist;
                    end
                end
            end
        end
        
        function plotMap(obj, map, xleft, xright, yleft, yright, n)
            figure;
            colormap(jet);
            imagesc(linspace(xleft, xright, n),linspace(yleft,yright,n),map);
            set(gca, 'YDir', 'normal');
            axis square;
            colorbar;
        end

        % Display voltage and electric field
        function plotVoltage(obj, scale_factor, n)
            [v, xleft, xright, yleft, yright] = obj.computeVoltage(scale_factor, n);
            
            obj.plotMap(v, xleft, xright, yleft, yright, n);
            title("Voltage [V]");
            xlabel("x [m]");
            ylabel("y [m]");

            % Compute electric field strength
            v_upper = v(1:end-1,1:end-1);
            E_x = (v(2:end, 1:end-1) - v_upper) / ((xright-xleft)/n);
            E_y = (v(1:end-1, 2:end) - v_upper) / ((yright-yleft)/n);

            E = sqrt(E_x .^ 2 + E_y .^ 2);
            
            obj.plotMap(E, xleft, xright, yleft, yright, n);
            title("Electric Field Intensity [V/m]");
            xlabel("x [m]");
            ylabel("y [m]");
        end

        % Solve for the charge distribution created by charging the mesh to
        % V0.
        function obj = solve(obj, V0)
            % Coefficient matrix
            A = zeros(obj.num_pts);
            % weighted voltage vector
            v = zeros(obj.num_pts,1);

            % compute points with neighbors included at beginning and end
            pts_neighbors = [obj.points(obj.num_pts,:); obj.points; obj.points(1,:)];
            
            % compute voltage vector
            for i = 1:1:obj.num_pts
                v(i) = V0 * obj.basis.magnitude(pts_neighbors(i+1,:), pts_neighbors(i, :), pts_neighbors(i+2, :));
            end

            % compute coefficient matrix
            for i = 1:1:obj.num_pts
                for j = 1:1:obj.num_pts
                    A(i,j) = obj.basis.innerProduct(pts_neighbors(i+1, :), pts_neighbors(i, :), pts_neighbors(i+2, :), pts_neighbors(j+1, :), pts_neighbors(j, :), pts_neighbors(j+2, :));
                end
            end
            
            % TODO: investigate why A is often singular or very ill
            % conditioned?
            obj.weights = lsqr(A,v,1e-6,128);
        end
        
        % Generate a set of points that are at the midpoint of the mesh
        % points
        function midpoints = midpointSet(obj)
            midpoints = zeros(obj.num_pts, 2);
            for i = 1:1:obj.num_pts
                midpoints(i,:) = 0.5 * (obj.points(i,:) + obj.points(mod(i, obj.num_pts)+1,:));
            end
        end

    end
end