% Methods of Moments for Electrostatics
% APPM3310 Final Project
% Edward Wawrzynek, Max Eaton, Andrew Zirger

function mom
    set(0,'defaultTextInterpreter','latex');

    % Example: solve the CU Logo 
    mesh1 = endpointsMesh(table2array(readtable("culogo.csv")), 300, BasisFunctions.Delta);
    mesh1.points = mesh1.points / 200.0;

    mesh2 = endpointsMesh(table2array(readtable("culogo.csv")), 300, BasisFunctions.Pulse);
    mesh2.points = mesh2.points / 200.0;

    mesh3 = endpointsMesh(table2array(readtable("culogo.csv")), 300, BasisFunctions.Triangle);
    mesh3.points = mesh3.points / 200.0;
    
    mesh1 = mesh1.solve(1.0);
    mesh2 = mesh2.solve(1.0);
    mesh3 = mesh3.solve(1.0);

    f = figure(1);

    subplot(1,3,1);
    mesh1.plotCharge();
    title("Delta");

    subplot(1,3,2);
    mesh2.plotCharge();
    title("Pulse");

    subplot(1,3,3);
    mesh3.plotCharge();
    title("Triange");

    h = axes(f, "visible", "off");
    h.Title.Visible = "on";
    h.XLabel.Visible = "on";
    h.YLabel.Visible = "on";
    xlabel(h, "x [m]");
    ylabel(h, "y [m]");
    title("Surface Charge Density [C]");

    c = colorbar(h, 'Position', [0.93 0.168 0.022 0.7]);
    colormap(c,'jet');
    caxis(h, [min(mesh2.weights) max(mesh2.weights)]);
end

% construct a mesh for a circle of radius 1 from the specified number of
% points
function mesh = circleMesh(num_pts, basis)
    pts = zeros(num_pts, 2);
    for i = 1:1:num_pts
        pts(i,:) = [cos(i/num_pts * 2*pi) sin(i/num_pts * 2*pi)];
    end

    mesh = Mesh(pts, basis);
end

function mesh = endpointsMesh(endpoints, numPts, basis) %creates a mesh of equally spaced points around a path defined by endpoints of line segments
    segArr = [endpoints(1:length(endpoints),:),[endpoints(2:length(endpoints),:); endpoints(1,:)]];
    segLen = 0;
    slopeArr = zeros(length(segArr),1);
    for i=1:1:length(segArr)
        segLen = segLen + ((segArr(i,1)-segArr(i,3))^2+(segArr(i,2)-segArr(i,4))^2)^0.5; %Adding the length of each segment to the total edge length
        slopeArr(i,1) = (segArr(i,4)-segArr(i,2))/(segArr(i,3)-segArr(i,1));
        if((segArr(i,4) == segArr(i,2)) && (segArr(i,3)-segArr(i,1) < 0))
            slopeArr(i,1) = -15; %Don't think about why this is
        end
    end
    meshDist = segLen/numPts; %Distance between points on the mesh
    mesh = zeros(numPts,2);
    mesh(1,:) = endpoints(1,:);
    curSeg = 1;
    for i=2:1:numPts
        newPt = mesh(i-1,:);
        distLeft = meshDist;
        while (distLeft > 0 && curSeg < 35)
            switch slopeArr(curSeg)
                case 0
                    if(segArr(curSeg,3)-newPt(1) > distLeft) 
                        newPt(1) = newPt(1) + distLeft;
                        distLeft = 0;
                    else
                        distLeft = distLeft - (segArr(curSeg,3)-newPt(1));
                        newPt(1) = newPt(1) + (segArr(curSeg,3)-newPt(1));
                        curSeg = curSeg+1;
                    end
                case -15 % This is -0 don't worry about it
                    if(newPt(1)-segArr(curSeg,3) > distLeft) 
                        newPt(1) = newPt(1) - distLeft;
                        distLeft = 0;
                    else
                        distLeft = distLeft - (newPt(1)-segArr(curSeg,3));
                        newPt(1) = newPt(1) - (newPt(1)-segArr(curSeg,3));
                        curSeg = curSeg+1;
                    end
                case 1
                    if(segArr(curSeg,3)-newPt(1) > 0) %Moving to the right and up
                        if((2^0.5)*(segArr(curSeg,3)-newPt(1)) > distLeft) 
                            newPt(1) = newPt(1) + (1/2^0.5)*distLeft;
                            newPt(2) = newPt(2) + (1/2^0.5)*distLeft;
                            distLeft = 0;
                        else
                            distLeft = distLeft - (2^0.5)*(segArr(curSeg,3)-newPt(1));
                            newPt(1) = newPt(1) + (segArr(curSeg,3)-newPt(1));
                            newPt(2) = newPt(2) + (segArr(curSeg,4)-newPt(2));
                            curSeg = curSeg+1;
                        end
                    else %Moving to the left and down
                        if((2^0.5)*(newPt(1)-segArr(curSeg,3)) > distLeft) 
                            newPt(1) = newPt(1) - (1/2^0.5)*distLeft;
                            newPt(2) = newPt(2) - (1/2^0.5)*distLeft;
                            distLeft = 0;
                        else
                            distLeft = distLeft - (2^0.5)*(newPt(1)-segArr(curSeg,3));
                            newPt(1) = newPt(1) - (newPt(1)-segArr(curSeg,3));
                            newPt(2) = newPt(2) - (newPt(2)-segArr(curSeg,4));
                            curSeg = curSeg+1;    
                        end
                        
                    end
                case -1
                    if(segArr(curSeg,3)-newPt(1) > 0) %Moving to the right and down
                        if((2^0.5)*(segArr(curSeg,3)-newPt(1)) > distLeft) 
                            newPt(1) = newPt(1) + (1/2^0.5)*distLeft;
                            newPt(2) = newPt(2) - (1/2^0.5)*distLeft;
                            distLeft = 0;
                        else
                            distLeft = distLeft - (2^0.5)*(segArr(curSeg,3)-newPt(1));
                            newPt(1) = newPt(1) + (segArr(curSeg,3)-newPt(1));
                            newPt(2) = newPt(2) - (newPt(2)-segArr(curSeg,4));
                            curSeg = curSeg+1;
                        end
                    else %Moving to the left and up
                        if((2^0.5)*(newPt(1)-segArr(curSeg,3)) > distLeft) 
                            newPt(1) = newPt(1) - (1/2^0.5)*distLeft;
                            newPt(2) = newPt(2) + (1/2^0.5)*distLeft;
                            distLeft = 0;
                        else
                            distLeft = distLeft - (2^0.5)*(newPt(1)-segArr(curSeg,3));
                            newPt(1) = newPt(1) - (newPt(1)-segArr(curSeg,3));
                            newPt(2) = newPt(2) - (newPt(2)-segArr(curSeg,4));
                            curSeg = curSeg+1;
                        end
                        
                    end
                case Inf
                    if(segArr(curSeg,4)-newPt(2) > distLeft) 
                        newPt(2) = newPt(2) + distLeft;
                        distLeft = 0;
                    else
                        distLeft = distLeft - (segArr(curSeg,4)-newPt(2));
                        newPt(2) = newPt(2) + (segArr(curSeg,4)-newPt(2));
                        curSeg = curSeg+1;
                    end
                case -Inf
                    if(newPt(2)-segArr(curSeg,4) > distLeft) 
                        newPt(2) = newPt(2) - distLeft;
                        distLeft = 0;
                    else
                        distLeft = distLeft - (newPt(2)-segArr(curSeg,4));
                        newPt(2) = newPt(2) - (newPt(2)-segArr(curSeg,4));
                        curSeg = curSeg+1;
                    end
                otherwise
                    disp("Something went wrong.")
            end
        end
        mesh(i,:) = newPt;
    end
    mesh = Mesh(mesh, basis);
end

% construct a mesh for a square of side length 1
function mesh = squareMesh(n, basis)
    pts = [ [linspace(0,1-1/n,n).' zeros(n,1)]; ...
            [ones(n,1) linspace(0,1-1/n,n).']; ...
            [linspace(1,1/n,n).' ones(n,1)]; ...
            [zeros(n,1) linspace(1,1/n,n).']
        ];
    
    mesh = Mesh(pts, basis);
end

% construct a mesh with rounded corner
function mesh = roundCornerMesh(n, radius, basis)
    pts = [ [linspace(-1+1/n,-1/n, n).' ones(n,1) .* radius]; ...
            [cos(linspace(pi/2,0,floor(n*pi*radius/2))).' .* radius sin(linspace(pi/2,0,floor(n*pi*radius/2))).' .* radius]; ...
            [ones(n,1) .* radius linspace(-1/n,-1+1/n, n).']; ...
            [cos(linspace(0,-pi/2,floor(n*pi*radius/2))).' .* radius sin(linspace(0,-pi/2,floor(n*pi*radius/2))).' .* radius - ones(floor(n*pi*radius/2),1)]; ...
            [linspace(-1/n,-1+1/n, n).' ones(n,1) .* -(1+radius)]; ...
            [(cos(linspace(-pi/2,-pi,floor(n*pi*radius/2))).' .* radius -ones(floor(n*pi*radius/2),1)) sin(linspace(-pi/2,-pi,floor(n*pi*radius/2))).' .* radius - ones(floor(n*pi*radius/2),1)]; ...
            [ones(n,1) .* -(1+radius) linspace(-1+1/n,-1/n, n).']; ...
            [(cos(linspace(pi,pi/2,floor(n*pi*radius/2))).' .* radius -ones(floor(n*pi*radius/2),1)) sin(linspace(pi,pi/2,floor(n*pi*radius/2))).' .* radius]; ...
            ];
    mesh = Mesh(pts, basis);
end