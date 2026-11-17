function test_Ground_Radiation()
%TEST_GROUND_RADIATION Synthetic tests; no external Element_Library required.
% A single triangular face at x=2 has area=2 m^2 and an outward +X normal.
% For up=+Z its unobstructed flat-ground view factor is 0.5.
% The stored dTe is deliberately reversed to test normal reorientation.

    E.nnp_Elemental = 3;
    E.Surf_ngp_Elemental = 1;
    E.Surf_Shape_Elemental = reshape([1/3,1/3,1/3],1,3,1,1);

    Mesh.nnp_Model = 3;
    Mesh.Nodes = [1 2 0 1; 2 2 1 1; 3 2 0 2];
    Mesh.IEN = [1;2;3];
    Mesh.Type = 1;
    Mesh.S_Type = [1 2; 1 1; 1 1]; % second face is NON-EXTERIOR

    FE.dTeS_gp = zeros(3,1,2);
    FE.dTeS_gp(:,1,1) = [-1;0;0];
    FE.dTeS_gp(:,1,2) = [-1;0;0];
    FE.T_surf_gp = zeros(3,1,2);
    FE.T_surf_gp(:,1,1) = 2/3;    % each node gets 2/3 m2
    FE.T_surf_gp(:,1,2) = 100/3;  % deliberately excluded surface

    [unit,diagnostics] = Ground_Radiation_Precompute( ...
        E,Mesh,FE,1,[0;0;1],[0;0;0]);
    assert(max(abs(unit - ones(3,1)/3)) < 1e-12);
    assert(abs(diagnostics.Exterior_Area_m2 - 2) < 1e-12);
    assert(abs(diagnostics.Ground_View_Weighted_Area_m2 - 1) < 1e-12);
    assert(abs(diagnostics.Mean_Ground_View_Factor - 0.5) < 1e-12);
    assert(diagnostics.Orientation_Flips == 1);

    sigma = 5.670374419e-8;
    f_ir = Ground_Radiation_IR_Load(unit,0.85,0.96,sigma,300);
    expected_ir = 0.85*0.96*sigma*300^4*1;
    assert(abs(sum(f_ir)-expected_ir) < 1e-8);

    f_albedo_noon = Ground_Radiation_Albedo_Load( ...
        unit,0.76,0.12,1450,[0;0;-1],[0;0;1]);
    assert(abs(sum(f_albedo_noon)-0.76*0.12*1450) < 1e-9);

    f_albedo_night = Ground_Radiation_Albedo_Load( ...
        unit,0.76,0.12,1450,[0;0;1],[0;0;1]);
    assert(max(abs(f_albedo_night)) == 0);

    fprintf('PASS: orientation, surface-group exclusion, area, IR, and albedo checks.\n');
end
