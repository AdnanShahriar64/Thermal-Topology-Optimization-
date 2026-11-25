function [F_Area,Info] = Ground_Radiation_Precompute( ...
    Mesh_Information,FE_Information_Surface,Exterior_Surface_Group,Up_Vector,Dome_Center)
%GROUND_RADIATION_PRECOMPUTE Ground-view-weighted exterior nodal areas [m^2].
% Uses the existing surface Gauss-point integration data without modifying
% Surface_Integration_3D_Flux_Thermal.
%
% F_Area(i) = integral_on_exterior N_i (1-n_out.dot(up))/2 dA.
% Assumes a convex dome above unobstructed, infinite, horizontal ground.
% Dome_Center is used to orient the exterior normal radially outwards.

if nargin<3 || isempty(Exterior_Surface_Group)
    Exterior_Surface_Group = 1;
end
if nargin<4 || isempty(Up_Vector)
    Up_Vector = [0;0;1];
end
if nargin<5 || isempty(Dome_Center)
    Dome_Center = [0;0;0];
end

up = Up_Vector(:);
center = Dome_Center(:);
if numel(up)~=3 || any(~isfinite(up)) || norm(up)==0 || ...
        numel(center)~=3 || any(~isfinite(center))
    error('GroundRadiation:InvalidGeometry','Up and center must be valid 3-vectors.');
end
up = up/norm(up);

F_Area = zeros(Mesh_Information.nnp_Model,1);
S_Type = Mesh_Information.S_Type;
Area_Exterior = 0;
Orientation_Flips = 0;
Gauss_Points = 0;

for sN = 1:size(S_Type,2)
    if S_Type(1,sN)~=Exterior_Surface_Group
        continue
    end
    e = S_Type(2,sN);
    nn = Mesh_Information.Type(2,e);
    node_ids = Mesh_Information.IEN(1:nn,e);
    coordinates = Mesh_Information.Nodes(node_ids,2:4);
    ngp = FE_Information_Surface.ngp_Surf(sN);

    for gp = 1:ngp
        dTe = FE_Information_Surface.dTeS_gp(:,gp,sN);
        magnitude = norm(dTe);
        if magnitude<=eps
            error('GroundRadiation:ZeroJacobian', ...
                'Zero surface area vector at face %d GP %d.',sN,gp);
        end
        normal = dTe/magnitude;

        % T_surf_gp already includes the Gauss weight and surface Jacobian.
        N_dA = FE_Information_Surface.T_surf_gp(1:nn,gp,sN);
        dA_weighted = sum(N_dA);
        if abs(dA_weighted)<=eps
            error('GroundRadiation:BadWeights','Zero quadrature weight at face %d.',sN);
        end

        % Recover physical Gauss-point coordinates from shape weights;
        % partition of unity makes the common dA factor cancel.
        x_gp = (N_dA.'*coordinates).' / dA_weighted;
        outward = x_gp-center;
        if norm(outward)>eps && dot(normal,outward)<0
            normal = -normal;
            Orientation_Flips = Orientation_Flips+1;
        end

        Fg = max(0,min(1,(1-dot(normal,up))/2));
        F_Area(node_ids) = F_Area(node_ids) + N_dA*Fg;
        Area_Exterior = Area_Exterior + dA_weighted;
        Gauss_Points = Gauss_Points+1;
    end
end

if Area_Exterior<=0
    error('GroundRadiation:NoExterior','No exterior surface faces found.');
end

Info.Exterior_Area_m2 = Area_Exterior;
Info.Ground_View_Weighted_Area_m2 = sum(F_Area);
Info.Mean_Ground_View_Factor = sum(F_Area)/Area_Exterior;
Info.Orientation_Flips = Orientation_Flips;
Info.Gauss_Points = Gauss_Points;
end
