function[Normal_Vector_form,TimeOut,T_Original,time_Store,T,Error]=Thermal_Model_Runner_Complete(Element_Information,Mesh_Information,Material_Information,T_One_Lunar_Day_Time,Initial_Time,Final_Time,T_Initial,Alpha,Maximum_allowed_temperature,T_ground,Check_rate,Time_to_avoid)

%% 1. Set up K, M, K_Rad, K_Conv

Dm=Material_Array_Generator_Heat_Transfer(Material_Information);

FE_Information_Region=Region_Integration_3D_Heat_Transfer(Element_Information,Mesh_Information);
FE_Information_Surface=Surface_Integration_3D_Flux_Thermal(Element_Information,Mesh_Information);
FE_Information_Surface_Solid=Surface_Integration_3D_Pressure_Type_Solid(Element_Information,Mesh_Information);

Ke=Elemental_Conductivity_Matrix_Generator(Element_Information,Mesh_Information,Material_Information,FE_Information_Region,Dm);
Me=Elemental_Capacitance_Matrix_Generator(Element_Information,Mesh_Information,Material_Information,FE_Information_Region);

Ke_surf=FE_Information_Surface.Ke_surf;
Ke_Radiation=Material_Information.epsilon*Material_Information.StefBC*Ke_surf(:,:,:,1);
Ke_Convection=Ke_surf(:,:,:,2);

K=Assembly(Ke,1,Mesh_Information,Element_Information);
M=Assembly(Me,1,Mesh_Information,Element_Information);
K_Rad=Assembly(Ke_Radiation,1,Mesh_Information,Element_Information);
K_Conv=Assembly(Ke_Convection,1,Mesh_Information,Element_Information);


%% 2. Normal vectors, solar area, ground IR and albedo

Nodal_Pressure=ones(Mesh_Information.nnp_Model,1);
fS=Pressure_to_Equivalent_Nodal_Force_Nodal_Pressure_Lag_GL_2(Element_Information,Mesh_Information,FE_Information_Surface_Solid,Nodal_Pressure);

Normal_Vector_form=zeros(size(fS));
for i=1:Mesh_Information.nnp_Model
    ii=3*i-2:3*i;
    if norm(fS(ii))>1e-10
        Normal_Vector_form(ii)=fS(ii)/norm(fS(ii));
    end
end

Nodal_Flux=ones(Mesh_Information.nnp_Model,1);
QS_Unity=Flux_to_Equivalent_Nodal_Flux_Lag_GL(Element_Information,Mesh_Information,FE_Information_Surface,Nodal_Flux);

% Uniform, unshadowed lunar ground. Use T_ground as ground IR temperature for now.
Material_Information.epsilon_gr=0.96;
Material_Information.ground_albedo=0.12;
Up_Vector=[0;0;1];

F_Area=Ground_Radiation_Precompute(Mesh_Information,FE_Information_Surface,Up_Vector);
f_ground_IR=Ground_Radiation_IR_Load(F_Area,T_ground,Material_Information);


%% 3. Boundary conditions

T_bc=ones(Mesh_Information.nnp_Model,1);
FTK=zeros(Mesh_Information.nnp_Model,size(Mesh_Information.S_Nodes,1));

for i=1:size(Mesh_Information.S_Nodes,1)
    Node_No=Mesh_Information.S_Nodes(i,3);
    if Node_No>0
        FTK(:,i)=-K(:,Node_No);

        K(:,Node_No)=0;
        K(Node_No,Node_No)=1;
        M(:,Node_No)=0;
        M(Node_No,Node_No)=1;
        K_Rad(:,Node_No)=0;
        K_Conv(:,Node_No)=0;

        T_bc(Node_No)=0;
    end
end

QS_Unity=QS_Unity.*T_bc;


%% 4. Truncate the fixed degrees of freedom

Ending_range=Mesh_Information.S_Nodes(1,3)-1;

In_M=inv(M(1:Ending_range,1:Ending_range));
In_MK=In_M*K(1:Ending_range,1:Ending_range);
K_Rad_Trun=K_Rad(1:Ending_range,1:Ending_range);
K_Conv_Trun=K_Conv(1:Ending_range,1:Ending_range); % For future convection

fgr=sum(FTK,2)*T_ground;          % Fixed-base conduction
fgr=fgr(1:Ending_range);
f_ground_IR=f_ground_IR(1:Ending_range);


%% 5. Initialize and run transient time loop (explicit Euler)

Inner_Nodes=Mesh_Information.S_Nodes(:,2);
Inner_Nodes=Inner_Nodes(Inner_Nodes>0);

Base_Nodes=Mesh_Information.S_Nodes(:,3);
Base_Nodes=Base_Nodes(Base_Nodes>0);

T=T_ground*ones(Mesh_Information.nnp_Model,1);
T(Base_Nodes)=0;                   % Base contribution enters through fgr
if norm(T_Initial)>0
    T=T_Initial;
    T(Base_Nodes)=0;
end

dt=2;
Time_Range=Initial_Time:dt:Final_Time;
TimeOut=Time_Range(end);

T_Original=zeros(length(T),ceil((Final_Time-Initial_Time)/Check_rate)+2);
time_Store=zeros(1,size(T_Original,2));
count=0;

Error=zeros(length(Inner_Nodes),1);
Detect=0;
Max_T=Maximum_allowed_temperature;

Rx=[1 0 0;0 cos(Alpha) -sin(Alpha);0 sin(Alpha) cos(Alpha)];

for t=Time_Range

    % 5a. Direct solar and ground albedo
    beta=pi*t/T_One_Lunar_Day_Time;
    Ry=[cos(beta) 0 sin(beta);0 1 0;-sin(beta) 0 cos(beta)];
    Sol_Vec=Rx*Ry*[1;0;0];

    SF=1450;
    if sin(beta)<0
        SF=0;
    end

    Flux=zeros(Mesh_Information.nnp_Model,1);
    for i=1:size(Mesh_Information.S_Nodes,1)
        n=Mesh_Information.S_Nodes(i,1);
        if n>0
            Normal=Normal_Vector_form(3*n-2:3*n);
            Component=-Sol_Vec'*Normal;
            if Component>0
                Flux(n)=Component*SF;
            end
        end
    end

    Flux=Flux.*T_bc;
    Q_SL=QS_Unity.*Flux*Material_Information.absorbability;

    f_albedo=Ground_Radiation_Albedo_Load(F_Area,SF,Sol_Vec,Up_Vector,Material_Information);

    % 5b. Explicit temperature update
    T_Trun=T(1:Ending_range);

    fsl=Q_SL(1:Ending_range);
    frd=-K_Rad_Trun*T_Trun.^4;
    f=fsl+frd+fgr+f_ground_IR+f_albedo(1:Ending_range);

    T_Trun=T_Trun+dt*(In_M*f-In_MK*T_Trun);
    T(1:Ending_range)=T_Trun;

    % 5c. Save temperatures and check temperature limit
    if mod(t,Check_rate)==0 && t>-dt
        count=count+1;

        T(Base_Nodes)=T_ground;
        T_Original(:,count)=T;
        time_Store(count)=t;

        T_Inner_Surface=T(Inner_Nodes).*T_bc(Inner_Nodes);

        if ~any(Time_to_avoid==t)
            for n=1:length(T_Inner_Surface)
                if T_Inner_Surface(n)>Max_T
                    Error(n)=T_Inner_Surface(n)-Max_T;
                    Detect=1;
                end
            end
        end

        if Final_Time<T_One_Lunar_Day_Time*2 && t<Final_Time
            Detect=0;
        end

        if mod(count,10)==0
            fprintf(' %3d, ',round(t/3600))
        end
    end

    if Detect==1
        TimeOut=t;
        break
    end
end

T_Original=T_Original(:,1:count);
time_Store=time_Store(1:count);

end
