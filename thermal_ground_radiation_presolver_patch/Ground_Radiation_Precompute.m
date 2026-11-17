function [F_Area, Info] = Ground_Radiation_Precompute( ...
    Element_Information, Mesh_Information, FE_Information_Surface, ...
    Exterior_Surface_Group, Up_Vector, Dome_Center)
%GROUND_RADIATION_PRECOMPUTE Consistent FE load from ground view factor.
%
%   [F_Area, Info] = Ground_Radiation_Precompute( ...
%       Element_Information, Mesh_Information, FE_Information_Surface, ...
%       Exterior_Surface_Group, Up_Vector, Dome_Center)
%
%   F_Area [m^2] is an nNodes-by-1 vector:
%      F_Area_i = integral_Gamma_outer N_i * F_ground(x) dGamma
%
%   For a convex surface above an unobstructed infinite horizontal plane:
%      F_ground(x) = (1 - dot(n_out(x), Up_Vector))/2.
%
%   This function REUSES the existing surface integration data WITHOUT
%   modifying Surface_Integration_3D_Flux_Thermal:
%       dTeS_gp(:,gp,face)  : oriented physical surface-area vector;
%       T_surf_gp(:,gp,face): w_gp * N_i(gp) * norm(dTe_gp).
%   It does not use Ke_surf directly because the ground view factor varies
%   from quadrature point to quadrature point. Both arrays are produced by
%   the same integration that produces Ke_surf.
%
%   NORMAL ORIENTATION: For the present spherical/hemispherical dome only,
%   the outward normal is chosen so dot(n, x_gp - Dome_Center)>0.
%   This is unsuitable for a general concave enclosure; use ray-traced
%   visibility/view factors for such geometry.
%
%   Units: Mesh coordinates must be in meters; normals are dimensionless.

    if nargin < 5 || isempty(Up_Vector)
        Up_Vector = [0; 0; 1];
    end
    if nargin < 6 || isempty(Dome_Center)
        Dome_Center = [0; 0; 0];
    end

    Up_Vector = Up_Vector(:);
    Dome_Center = Dome_Center(:);
    if numel(Up_Vector) ~= 3 || numel(Dome_Center) ~= 3 || ...
            any(~isfinite(Up_Vector)) || any(~isfinite(Dome_Center)) || ...
            norm(Up_Vector) == 0
        error('GroundRadiation:InvalidGeometry', ...
            'Up_Vector and Dome_Center must be finite 3-vectors; Up_Vector cannot be zero.');
    end
    Up_Vector = Up_Vector / norm(Up_Vector);

    if ~isscalar(Exterior_Surface_Group) || ...
            Exterior_Surface_Group < 1 || ...
            Exterior_Surface_Group ~= floor(Exterior_Surface_Group)
        error('GroundRadiation:InvalidSurface', ...
            'Exterior_Surface_Group must be a positive integer.');
    end

    required = {'dTeS_gp', 'T_surf_gp'};
    for k = 1:numel(required)
        if ~isfield(FE_Information_Surface, required{k})
            error('GroundRadiation:MissingField', ...
                'Surface integration result is missing %s.', required{k});
        end
    end

    S_Type = Mesh_Information.S_Type;     % row 1: group; 2: element; 3: local face
    surface_indices = find(S_Type(1,:) == Exterior_Surface_Group);
    if isempty(surface_indices)
        error('GroundRadiation:NoExteriorSurface', ...
            'No surfaces were found for group %d.', Exterior_Surface_Group);
    end

    nnp = Mesh_Information.nnp_Model;
    F_Area = zeros(nnp, 1);
    total_area = 0;
    total_view_area = 0;
    orientation_flips = 0;
    num_gp = 0;

    for i = surface_indices
        e = S_Type(2,i);
        face = S_Type(3,i);
        elem_type = Mesh_Information.Type(1,e);
        nn = Element_Information.nnp_Elemental(elem_type);
        node_ids = Mesh_Information.IEN(1:nn,e);
        xyz = Mesh_Information.Nodes(node_ids,2:4);

        nqp = Element_Information.Surf_ngp_Elemental(face,elem_type);
        for gp = 1:nqp
            dTe = FE_Information_Surface.dTeS_gp(:,gp,i);
            area_jacobian = norm(dTe);
            if ~isfinite(area_jacobian) || area_jacobian <= 0
                error('GroundRadiation:DegenerateSurface', ...
                    'Invalid dTe at S_Type column %d, element %d, GP %d.', i,e,gp);
            end

            % Evaluate position at the SAME Gauss point used by dTeS_gp.
            shape = reshape( ...
                Element_Information.Surf_Shape_Elemental( ...
                gp,1:nn,face,elem_type), 1, nn);
            x_gp = (shape * xyz).';
            radial = x_gp - Dome_Center;
            if norm(radial) <= 1e-12 || ~all(isfinite(radial))
                error('GroundRadiation:UndefinedOutward', ...
                    'Cannot identify outward direction at element %d, GP %d.',e,gp);
            end

            % Outward radial is valid for the present convex dome.
            normal = dTe / area_jacobian;
            alignment = dot(normal, radial) / norm(radial);
            if abs(alignment) < 1e-8
                error('GroundRadiation:AmbiguousNormal', ...
                    'Surface normal is nearly tangential to dome radius at element %d, GP %d.',e,gp);
            end
            if alignment < 0
                normal = -normal;
                orientation_flips = orientation_flips + 1;
            end

            F_ground = max(0, min(1, (1 - dot(normal,Up_Vector))/2));

            % T_surf_gp already contains N_i * GP_weight * surface Jacobian.
            weighted_shape = FE_Information_Surface.T_surf_gp(1:nn,gp,i);
            if any(~isfinite(weighted_shape))
                error('GroundRadiation:InvalidQuadrature', ...
                    'Nonfinite nodal quadrature weights at element %d, GP %d.',e,gp);
            end
            F_Area(node_ids) = F_Area(node_ids) + weighted_shape * F_ground;
            dA = sum(weighted_shape);
            total_area = total_area + dA;
            total_view_area = total_view_area + dA * F_ground;
            num_gp = num_gp + 1;
        end
    end

    Info.Exterior_Surface_Group = Exterior_Surface_Group;
    Info.Exterior_Area_m2 = total_area;
    Info.Ground_View_Weighted_Area_m2 = total_view_area;
    Info.Integrated_Load_Area_m2 = sum(F_Area);
    Info.Mean_Ground_View_Factor = total_view_area / total_area;
    Info.Orientation_Flips = orientation_flips;
    Info.Gauss_Points = num_gp;
    Info.Up_Vector = Up_Vector;
    Info.Dome_Center = Dome_Center;
end
