function [Normal_Vector_form,TimeOut,T_Original,time_Store,T,Error] = Thermal_Model_Runner_Complete( ...
    Element_Information,Mesh_Information,Material_Information, ...
    T_One_Lunar_Day_Time,Initial_Time,Final_Time,T_Initial,Alpha, ...
    Maximum_allowed_temperature,T_ground,Check_rate,Time_to_avoid)
%THERMAL_MODEL_RUNNER_COMPLETE Transient 3D dome heat-transfer FEA.
% Five sections follow the calculation from assembly to transient response.
%
% Model: M*Tdot + K*T + K_Rad*(T.^4) = Q_direct + Q_groundIR
%                                             + Q_albedo + Q_convection.
% Prescribed base temperatures are included through K(free,:)*T.
% Temperatures are in K, heat loads in W, time in seconds.
% No ground shadow, internal radiation, or air convection is modeled here.

%% 1. PREPROCESS AND ASSEMBLE THE THERMAL MATRICES
% K: conduction (W/K); M: thermal capacity (J/K);
% K_Rad: outer-wall emission coefficient; K_Conv: inner-wall area matrix.

Dm = Material_Array_Generator_Heat_Transfer(Material_Information);

FE_Information_Region = Region_Integration_3D_Heat_Transfer( ...
    Element_Information,Mesh_Information);
FE_Information_Surface = Surface_Integration_3D_Flux_Thermal( ...
    Element_Information,Mesh_Information);
FE_Information_Surface_Solid = Surface_Integration_3D_Pressure_Type_Solid( ...
    Element_Information,Mesh_Information);

Ke = Elemental_Conductivity_Matrix_Generator( ...
    Element_Information,Mesh_Information,Material_Information, ...
    FE_Information_Region,Dm);
Me = Elemental_Capacitance_Matrix_Generator( ...
    Element_Information,Mesh_Information,Material_Information, ...
    FE_Information_Region);

epsilon = Material_Information.epsilon;
sigma = Material_Information.StefBC;
Ke_surf = FE_Information_Surface.Ke_surf;
Ke_Radiation = epsilon*sigma*Ke_surf(:,:,:,1);  % exterior S4, group 1
Ke_Convection = Ke_surf(:,:,:,2);               % interior S2, group 2

% Keep full matrices intact until partitioning into free/fixed DOFs.
K = sparse(Assembly(Ke,1,Mesh_Information,Element_Information));
M = sparse(Assembly(Me,1,Mesh_Information,Element_Information));
K_Rad = sparse(Assembly(Ke_Radiation,1,Mesh_Information,Element_Information));
K_Conv = sparse(Assembly(Ke_Convection,1,Mesh_Information,Element_Information));

nNodes = Mesh_Information.nnp_Model;

%% 2. EXTERIOR NORMALS AND INCOMING RADIATION LOAD DATA
% 2a. Nodal normals for existing direct-solar and geometry-update routines.
Nodal_Pressure = ones(nNodes,1);
fS = Pressure_to_Equivalent_Nodal_Force_Nodal_Pressure_Lag_GL_2( ...
    Element_Information,Mesh_Information, ...
    FE_Information_Surface_Solid,Nodal_Pressure);

Normal_Vector_form = zeros(size(fS));
for node = 1:nNodes
    idx = (3*node-2):(3*node);
    areaNormal = fS(idx);
    if norm(areaNormal) > 1e-10
        Normal_Vector_form(idx) = areaNormal/norm(areaNormal);
    end
end

% 2b. Equivalent nodal exterior area, only surface group 1.
% Unlike the old all-surfaces Flux_to_Equivalent... call, this does not
% inadvertently count inner/base surface area in direct solar heating.
QS_Unity = zeros(nNodes,1);
outerFaceIDs = find(Mesh_Information.S_Type(1,:) == 1);
for sf = outerFaceIDs
    elem = Mesh_Information.S_Type(2,sf);
    nn = Mesh_Information.Type(2,elem);
    nodes = Mesh_Information.IEN(1:nn,elem);
    ngp = FE_Information_Surface.ngp_Surf(sf);
    for gp = 1:ngp
        QS_Unity(nodes) = QS_Unity(nodes) + ...
            FE_Information_Surface.T_surf_gp(1:nn,gp,sf);
    end
end

% 2c. Ground geometry is fixed during a single transient simulation.
% These parameters are defaults; callers may override via Material_Information.
Ground_Material = Material_Information;
if ~isfield(Ground_Material,'epsilon_gr')
    Ground_Material.epsilon_gr = 0.96;
end
if ~isfield(Ground_Material,'ground_albedo')
    Ground_Material.ground_albedo = 0.12;
end
if isfield(Ground_Material,'T_ground_IR')
    T_ground_IR = Ground_Material.T_ground_IR;
else
    T_ground_IR = T_ground;  % Placeholder, distinct from the base BC.
end

Up_Vector = [0;0;1];       % mesh is rotated to +Z up in Optimizer.m
Dome_Center = [0;0;0];
Exterior_Surface_Group = 1;

[F_Area,Ground_View_Info] = Ground_Radiation_Precompute( ...
    Mesh_Information,FE_Information_Surface, ...
    Exterior_Surface_Group,Up_Vector,Dome_Center);

% Ground longwave infrared is constant for this first uniform-ground model.
f_ground_IR = Ground_Radiation_IR_Load(F_Area,T_ground_IR,Ground_Material);

% Ground-reflected shortwave (albedo) changes with the solar zenith angle.
% It is evaluated inside the time loop using current SF and Sol_Vec.
% Ground_View_Info provides integrated area diagnostics if needed.

%% 3. PRESCRIBED-TEMPERATURE BOUNDARY CONDITIONS
% S_Nodes(:,3) are fixed-temperature base nodes (group S3).
% Do not zero columns of K/M/K_Rad: keep all coupling to the base.
fixed = unique(Mesh_Information.S_Nodes(:,3));
fixed = fixed(fixed>0);
free = setdiff((1:nNodes)',fixed);

if isempty(free)
    error('ThermalRunner:NoFreeNodes','All thermal nodes are fixed.');
end
if any(fixed>nNodes)
    error('ThermalRunner:BadFixedNodes','Fixed node index exceeds mesh size.');
end

inner_nodes = Mesh_Information.S_Nodes(:,2);
inner_nodes = inner_nodes(inner_nodes>0);

T = T_ground*ones(nNodes,1);
if ~isempty(T_Initial) && any(T_Initial(:)~=0)
    if numel(T_Initial) ~= nNodes
        error('ThermalRunner:BadInitialT','T_Initial must have one value per node.');
    end
    T = T_Initial(:);
end
T(fixed) = T_ground;

%% 4. REDUCE TO FREE THERMAL DOFS
% This partition works even if the fixed nodes are not consecutive.
% Maintain full T for fixed-node conduction and nonlinear radiation.
Mff = M(free,free);
Kff = K(free,free); %#ok<NASGU> % retained for future solver changes
Kfree = K(free,:);
KradFree = K_Rad(free,:);
KconvFree = K_Conv(free,:); %#ok<NASGU> % convection presently disabled

% Factor once. Applying Mff\RHS preserves explicit Euler but avoids inv(M).
MassSolver = decomposition(Mff,'chol');

%% 5. INITIALIZE AND RUN THE EXPLICIT TRANSIENT SOLVER
dt = 2;   % seconds, unchanged from the original solver
Time_Range = Initial_Time:dt:Final_Time;
if isempty(Time_Range)
    error('ThermalRunner:BadTimeRange','Final_Time must exceed Initial_Time.');
end
if Check_rate<=0
    error('ThermalRunner:BadCheckRate','Check_rate must be positive.');
end

TimeOut = Time_Range(end);
Error = zeros(numel(inner_nodes),1);
Detect = false;

% Avoid allocating one column per 2-second timestep.
nSnapshots = ceil((Final_Time-Initial_Time)/Check_rate) + 3;
T_Original = zeros(nNodes,nSnapshots);
time_Store = zeros(1,nSnapshots);
saved = 0;

% Same solar trajectory and direct-normal solar flux as the old solver.
Rx = [1 0 0;0 cos(Alpha) -sin(Alpha);0 sin(Alpha) cos(Alpha)];
Max_T = Maximum_allowed_temperature;
Heat_transfer_coeff = 0; %#ok<NASGU> % convection retained as future option

for t = Time_Range
    %% 5a. Sun position and incoming heat loads
    beta = pi*t/T_One_Lunar_Day_Time;
    Ry = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
    Sol_Vec = Rx*Ry*[1;0;0]; % direction sunlight travels

    SF = 1450; % W/m^2, direct-normal illumination
    if sin(beta)<0
        SF = 0;
    end

    Flux = zeros(nNodes,1);
    for node = Mesh_Information.S_Nodes(:,1)'
        if node>0
            n = Normal_Vector_form(3*node-2:3*node);
            incidence = -dot(Sol_Vec,n);
            if incidence>0
                Flux(node) = SF*incidence;
            end
        end
    end
    f_solar_direct = Material_Information.absorbability* ...
        (QS_Unity.*Flux);

    f_albedo_now = Ground_Radiation_Albedo_Load( ...
        F_Area,SF,Sol_Vec,Up_Vector,Ground_Material);

    %% 5b. Transient heat balance on free nodes
    % Fixed base temperatures are enforced through the full temperature
    % vector, so there is no separate fgr term (avoids double counting).
    f_radiation_out = -KradFree*(T.^4);
    f_conduction = -Kfree*T;

    % No wall-to-air convection in the original model (h = 0).
    RHS = f_solar_direct(free) + f_ground_IR(free) + ...
          f_albedo_now(free) + f_radiation_out + f_conduction;

    T(free) = T(free) + dt*(MassSolver\RHS);
    T(fixed) = T_ground;

    if any(~isfinite(T(free)))
        error('ThermalRunner:Unstable','Nonfinite temperatures; check dt/stability.');
    end

    %% 5c. Save results and check the inner-surface temperature threshold
    if mod(t,Check_rate)==0 && t>-dt
        saved = saved+1;
        if saved>size(T_Original,2)
            T_Original(:,saved) = 0;
            time_Store(saved) = 0;
        end
        T_Original(:,saved) = T;
        time_Store(saved) = t;

        % Re-evaluate rather than retaining old exceedances.
        Error(:) = 0;
        if ~ismember(t,Time_to_avoid)
            innerT = T(inner_nodes);
            above = innerT>Max_T;
            Error(above) = innerT(above)-Max_T;
        end
        Detect = any(Error>0);

        % Preserve the original optimizer's partial-window behavior.
        if Final_Time<T_One_Lunar_Day_Time*2 && t<Final_Time
            Detect = false;
        end
        if mod(saved,10)==0
            fprintf(' %3d hr, ',round(t/3600));
        end
    end

    if Detect
        TimeOut = t;
        break
    end
end

T_Original = T_Original(:,1:saved);
time_Store = time_Store(1:saved);
end
