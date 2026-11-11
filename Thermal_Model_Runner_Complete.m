function[Normal_Vector_form,TimeOut,T_Original,time_Store,T,Error]=Thermal_Model_Runner_Complete(Element_Information,Mesh_Information,Material_Information,T_One_Lunar_Day_Time,Initial_Time,Final_Time,T_Initial,Alpha,Maximum_allowed_temperature,T_ground,Check_rate,Time_to_avoid)


      [Dm]=Material_Array_Generator_Heat_Transfer(Material_Information);
    
      % [Dm_Solid]=Material_Array_Generator(Material_Information);
      
%% System preprocessor
    
    [FE_Information_Region]=Region_Integration_3D_Heat_Transfer(Element_Information,Mesh_Information);
    [FE_Information_Surface]=Surface_Integration_3D_Flux_Thermal(Element_Information,Mesh_Information);

    % [FE_Information_Region_Solid]=Region_Integration_3D_Solid(Element_Information,Mesh_Information);
    [FE_Information_Surface_Solid]=Surface_Integration_3D_Pressure_Type_Solid(Element_Information,Mesh_Information);
    
    % [IS_F]=Internal_Stress_to_Nodal_FV(Element_Information,Mesh_Information,FE_Information_Region_Solid,Dm_Solid);
    
    % fprintf('Initial System building completed--------------\n')    
    
    
 %%
 
 %% Mat property and Stiffness    
 
[Ke]=Elemental_Conductivity_Matrix_Generator(Element_Information,Mesh_Information,Material_Information,FE_Information_Region,Dm);
[Me]=Elemental_Capacitance_Matrix_Generator(Element_Information,Mesh_Information,Material_Information,FE_Information_Region);


Ke_surf_per_h=FE_Information_Surface.Ke_surf;
epsilon=Material_Information.epsilon;
sigma=Material_Information.StefBC;
Ke_Radiation=epsilon*sigma*(Ke_surf_per_h(:,:,:,1));%+Ke_surf_per_h(:,:,:,2)*0); % Internal surface omitted

Ke_Convection=Ke_surf_per_h(:,:,:,2);

% Assembly    
    [K] = Assembly(Ke,1,Mesh_Information,Element_Information);
    [M] = Assembly(Me,1,Mesh_Information,Element_Information);

    
    % sum(sum(sum(sum(Ke_Radiation))))
    
    [K_Rad] = Assembly(Ke_Radiation,1,Mesh_Information,Element_Information);

    [K_Conv] = Assembly(Ke_Convection,1,Mesh_Information,Element_Information);



% Run it with matrix gen

Nodal_Pressure=Mesh_Information.Nodes(:,1)*0+1;
[fS]=Pressure_to_Equivalent_Nodal_Force_Nodal_Pressure_Lag_GL_2(Element_Information,Mesh_Information,FE_Information_Surface_Solid,Nodal_Pressure);

Nodes=Mesh_Information.Nodes;
% Normal=zeros(Mesh_Information.nnp_Model,Mesh_Information.nnp_Model*3);
Normal_Vector_form=fS*0;
for i=1:length(Nodes)
    if norm(fS(3*i-2:3*i,1))>10^-10
        % Normal(i,3*i-2:3*i-0)=fS(3*i-2:3*i,1)'/norm(fS(3*i-2:3*i,1));
        Normal_Vector_form(3*i-2:3*i,1)=fS(3*i-2:3*i,1)/norm(fS(3*i-2:3*i,1));
    end
end


Nodal_Flux=Mesh_Information.Nodes(:,1)*0+1;
[QS_Unity]=Flux_to_Equivalent_Nodal_Flux_Lag_GL(Element_Information,Mesh_Information,FE_Information_Surface,Nodal_Flux);

 
    T_bc=ones(Mesh_Information.nnp_Model,1);
    FTK=zeros(Mesh_Information.nnp_Model,length(Mesh_Information.S_Nodes));
    FTM=FTK;
    Surf_No=3;
    for i=1:length(Mesh_Information.S_Nodes)
        
        Node_No=Mesh_Information.S_Nodes(i,Surf_No);
        if Node_No>0

            FTK(:,i)=-K(:,Node_No);
            FTM(:,i)=-M(:,Node_No);
            
            K(:,Node_No)=K(:,Node_No)*0;
            K(Node_No,Node_No)=1;
            
            M(:,Node_No)=M(:,Node_No)*0;
            M(Node_No,Node_No)=1;

            K_Rad(:,Node_No)=K_Rad(:,Node_No)*0;

            K_Conv(:,Node_No)=K_Conv(:,Node_No)*0;
        
            T_bc(Node_No,1)=0;
        end
        
    end

    QS_Unity=QS_Unity.*T_bc;

%% Solver


    % T_One_Lunar_Day_Time=3600*24*15;
    % Final_Time=T_One_Lunar_Day_Time*2;

    Ending_range=Mesh_Information.S_Nodes(1,3)-1;

    In_M=inv(M(1:Ending_range,1:Ending_range));    % Inverse of Capacitance matrix
    In_MK=In_M*K(1:Ending_range,1:Ending_range);
    K_Rad_Trun=K_Rad(1:Ending_range,1:Ending_range);
    K_Conv_Trun=K_Conv(1:Ending_range,1:Ending_range);

    FTK_Ones_Unt=FTK*ones(size(FTK,2),1);
    FTK_Ones=FTK_Ones_Unt(1:Ending_range,:);

    
    % T_ground=254.8;%Min_T;
    
    Max_T=Maximum_allowed_temperature;%T_ground+50;
    Min_T=T_ground-50;%273-25;

    
    T_inner=T_ground;%(Max_T+Min_T)/2; % Min_T;

    fgr=FTK_Ones*T_ground;

    T=zeros(Mesh_Information.nnp_Model,1)+T_ground;
    
    Got_One=0;
    for i=1:length(Mesh_Information.S_Nodes)
        if Mesh_Information.S_Nodes(i,2)==0
            Got_One=1;
            break
        end
    end
    if Got_One==1
        Inner_Node_count=i-1;    
    else
        Inner_Node_count=i;
    end
    
    T(Mesh_Information.S_Nodes(1:Inner_Node_count,2))=T_inner;
    
    for i=1:length(Mesh_Information.S_Nodes)
        if Mesh_Information.S_Nodes(i,3)==0
            break
        end
    end
    Base_Node_count=i-1;
    T(Mesh_Information.S_Nodes(1:Base_Node_count,3))=0;

    if norm(T_Initial)>0
        T=T_Initial;
    end


    dt=2;
    count=1;
    
    Time_Range=Initial_Time:dt:Final_Time;
    
    % time_Store=zeros(length(Time_Range),1);

    Detect=0;
    
    Error=zeros(length((Mesh_Information.S_Nodes(1:Inner_Node_count,2))),1);

    TimeOut=Time_Range(end);

    
    T_Original=zeros(length(T),round(length(Time_Range)/(Check_rate)));
    time_Store=zeros(1,round(length(Time_Range)/(Check_rate)));

    count=1;

    %0.7854; % 67.5*pi/180;%
    Rx=[1 0 0; 0 cos(Alpha) -sin(Alpha);0 sin(Alpha) cos(Alpha)];

    T_air=T_ground;
    Mass_air=4/6*pi*2.5^3;
    Heat_transfer_coeff=0;

    

    % tic
    for t=Time_Range   

        %% Input solar flux

        Flux=zeros(Mesh_Information.nnp_Model,1);
        beta = pi*t/T_One_Lunar_Day_Time;
        Ry=[cos(beta)   0          sin(beta)
            0           1          0
           -sin(beta)   0          cos(beta)];

        Sn=1;
        Sol_Vec=Rx*Ry*[1;0;0];

        SF=1450;
        if sin(beta)<0
            SF=0;
        end
        

        for i=1:size(Mesh_Information.S_Nodes,1)
            n=Mesh_Information.S_Nodes(i,Sn);
            if n>0
                Normal_Component=Normal_Vector_form(3*n-2:3*n);
        
                Component=Sol_Vec'*-Normal_Component;
        
                if Component>0 
                   Flux(n)=Component*SF;
                end
            end
        end
        
    
        Flux=Flux.*T_bc;
        Q_SL=QS_Unity.*Flux*Material_Information.absorbability;


        %% Explicit time integration

        T_Trun=T(1:Ending_range);

        fsl=Q_SL(1:Ending_range);

        frd=-K_Rad_Trun*T_Trun.^4;

%         T_air=T_air+dt/Mass_air*Heat_transfer_coeff*(T_Trun-T_air);
% 
%         if mod(t,3600)==0
%             T_air
%             pause
%         end


        % fcon=-Heat_transfer_coeff*K_Conv_Trun*(T_Trun-T_air);
        
        f=fsl+frd+fgr;% fcon;
        
        T_Trun=T_Trun+dt*(In_M*f-In_MK*T_Trun);

        T(1:Ending_range)=T_Trun;
        

        % norm(T)

        if mod(t,Check_rate)==0 && t>-dt

            if count==10
                fprintf('Time (Hr): ')
            end

            if mod(count,10)==0 %|| (t/(60*60)>358 && t/(60*60)<362) || (t/(60*60)>-2 && t/(60*60)<2)% && Final_Time>400*60*60
                % fprintf('Time=%5d hour, Flux norm=%f\n',t/(60*60),norm(Q_SL))
                fprintf(' %3d, ',round(t/(60*60)))
            end

            time_Store(count)=t;
            % T_Store(:,count)=T;    

            % SF_Store(count)=SF;

            % frd_Store(:,count)=frd;

            T(Mesh_Information.S_Nodes(1:Base_Node_count,3))=T_ground;

            T_Original(:,count)=T;
            time_Store(count)=t;
            count=count+1;            

            T_bc_added=T.*T_bc;
            T_Inner_Surface=T_bc_added(Mesh_Information.S_Nodes(1:Inner_Node_count,2));

            if length(intersect(Time_to_avoid,t))==0

            Pos=0;
            Neg=0;
            for n=1:length(T_Inner_Surface)
                if T_Inner_Surface(n)>0
                    if T_Inner_Surface(n)>Max_T
                        Error(n)=T_Inner_Surface(n)-Max_T;
                        Detect=1;
                        Pos=1;
                    % Activate for min T
%                     elseif T_Inner_Surface(n)<Min_T
%                         Error(n)=T_Inner_Surface(n)-Min_T;
%                         Detect=1;
%                         Neg=1;
                     end
                end
                
            end
            end
            
            
%             if sum(abs(Error))>1
%                 Detect=1;
%             end
            
            if Final_Time<T_One_Lunar_Day_Time*2 && t<Final_Time
                Detect=0;
            end
                
            % if length(find(abs(Error)>1))==length(Error)
            %    Detect=1;
            % end
            
            

        end

        if Detect==1
            % fprintf('Time=%5d hour\n',t/(60*60))
            TimeOut=t;
            % fprintf('Temperature>Max=%d, Temperature<Min=%d\n',Pos,Neg);
            break
        end    
        
        
        % for e=1:nel
        %    T_e=T(IEN(:,e),1); 
        %    T_Gp(:,e)=Shape(:,1:20,7)*T_e;
        % end
        
        

%         % T_GP_S(:,:,c+1)=T_Gp;
%         count=count+1;
%         % T
    end
end
