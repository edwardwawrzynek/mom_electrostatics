% Methods of Moments for Electrostatics
% APPM3310 Final Project
% Edward Wawrzynek, Max Eaton, Andrew Zirger

function mom
    set(0,'defaultTextInterpreter','latex');
    mesh = endpointsMesh(table2array(readtable("culogo.csv")), 1000, BasisFunctions.Delta);
    mesh.points = mesh.points / 200.0;
    %mesh = circleMesh(1000, BasisFunctions.Delta);
    mesh.plotMesh();

    mesh = mesh.solve(1.0);
    mesh.plotCharge();
    mesh.plotVoltage(1.7, 350);
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