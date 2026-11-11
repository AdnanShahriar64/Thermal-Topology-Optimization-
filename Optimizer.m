clear
clc

% Run from this folder: the Abaqus input files are read from the current folder.
Main_Dir=fileparts(mfilename('fullpath'));
addpath(fullfile(Main_Dir,'Element_Library'), ...
        fullfile(Main_Dir,'System_Preprocessor_Files'), ...
        fullfile(Main_Dir,'Runtime_Calculation_Files'));
% Assign proper dimension
    

    Alpha=0;

    Alpha=0.3927;%0.3927;

    Alpha=0.7854;

    % [0.3927 0.7854 ]

    % Alpha=67.5000;

    % for Alpha=[0.3927 0.7854 1.2654]

    AlphaCount=0;

    for Alpha=[0 0.1963    0.3927    0.5890    0.7854    0.9817    1.1781    1.3744]

    % for Alpha=[0.3927    0.5890    0.7854    0.9817    1.1781    1.3744]

    % for Alpha=0.3927

    % for Alpha=[0.7854 1.1781]

        clearvars -except Alpha AlphaCount

    Inner_Radius=2.5;
    Intial_thickness=0.1;

    RHO=2000;
    SPCH=740;
    KX=0.7;
    EPSILON=0.85;
    ABSORBIBILITY=0.76;

    T_ground=254.5;
    Maximum_allowed_temperature=254.5+50;

    CHECK_RATE=60*60*10;

    DELX=0.001;

    LAMDA=0.9;
    
%% Input FE elements

    [Element_Information]=Element_Loader;

%% Material input      
      Material_Information.MoE=70*10^9;
      Material_Information.PoRat=0.3;
      Material_Information.alpe=10^-6;

      Material_Information.rho=RHO;
      Material_Information.SpcH=SPCH;
      Material_Information.Kx=KX;
      
      Material_Information.epsilon=EPSILON;
      Material_Information.absorbability=ABSORBIBILITY;


      Material_Information.StefBC=5.67*10^-8; 


%% Model Input File
% Read Abaqus Input file

      Name='SPL_Section_Linear_2.inp';
      % Name='SPL_Thermal-Copy-Copy2.inp';
      % Name='SPL_Thermal-Copy-Copy2-Copy.inp';
      % Name='SPL_Linear_3.inp';
      S_Names={'S4','S2','S3'};
      % S_Names={'S4'};
      Reg_Name={'Sx','S1'};
      BC_Name={'S3'};

      [Mesh_Information]=Abaqus_File_Processor_4(Name,S_Names,Reg_Name,BC_Name,Element_Information);



%%


    Outer_Radius=Inner_Radius+Intial_thickness;
    Middle_Radius=(Outer_Radius+Inner_Radius)/2;

    In_O_R=2.9;
    In_I_R=2.5;
    In_M_R=2.7;
    
    Nodes=Mesh_Information.Nodes;

    try

    for i=1:length(Nodes)
        
       if Nodes(i,3)>10^-2
          Radius=norm(Nodes(i,2:4)); 
           if abs(Radius-In_O_R)<10^-4
               Nodes(i,2:4)=Nodes(i,2:4)/Radius*Outer_Radius;
           elseif abs(Radius-In_I_R)<10^-4
               Nodes(i,2:4)=Nodes(i,2:4)/Radius*Inner_Radius;
           elseif abs(Radius-In_M_R)<10^-4
               Nodes(i,2:4)=Nodes(i,2:4)/Radius*Middle_Radius;          
           else
               problem
           end
           
       else
           
          Radius=norm(Nodes(i,[2 4]));
          
          Z=abs(Nodes(i,3));
          Zmin=abs(min(Nodes(:,3)));
          Mult=1+Z/Zmin*0;
          
          
           if abs(Radius-In_O_R)<10^-4
               Nodes(i,[2 4])=Nodes(i,[2 4])/Radius*Outer_Radius*Mult;
           elseif abs(Radius-In_I_R)<10^-4
               Nodes(i,[2 4])=Nodes(i,[2 4])/Radius*Inner_Radius;
           elseif abs(Radius-In_M_R)<10^-4
               Nodes(i,[2 4])=Nodes(i,[2 4])/Radius*Middle_Radius;          
           else
               problem
           end   

            Nodes(i,3)=Nodes(i,3)*60;

       end
       
       
    end
    Mesh_Information.Nodes=Nodes;

    end

    %% Rotation
     Nodes=Mesh_Information.Nodes;
     AlphaR=-90*pi/180;
     BetaR=0;
     GammaR=0;
     
     Rz=[cos(GammaR) -sin(GammaR) 0;sin(GammaR) cos(GammaR) 0; 0 0 1];
     Ry=[cos(BetaR) 0 sin(BetaR);0 1 0;-sin(BetaR) 0 cos(BetaR)];
     Rx=[1 0 0; 0 cos(AlphaR) -sin(AlphaR);0 sin(AlphaR) cos(AlphaR)];
     R1=Rx*Ry*Rz;
     
     for i=1:length(Nodes)
         Nodes(i,2:4)=inv(R1)*[Nodes(i,2);Nodes(i,3);Nodes(i,4)];
     end    
    
    
    Mesh_Information.Nodes=Nodes;
   
%% Mesh conversion

    [Mesh_Information]=Geometry_Model_Parameters_Preprocessor_Solid(Mesh_Information,Element_Information);
    [Mesh_Information]=Element_Converter_8_to_Any(Mesh_Information,Element_Information,10);
    [Mesh_Information]=Geometry_Model_Parameters_Preprocessor_Solid(Mesh_Information,Element_Information);
   
    
%%  Correction for boundary

try 

S_Nodes=Mesh_Information.S_Nodes;

count=1;
for i=1:length(Mesh_Information.Nodes)
    z=Mesh_Information.Nodes(i,4);
    if z>-10^-3
    else
        S_Nodes(count,3)=i;
        count=count+1;
    end
end

      Mesh_Information.S_Nodes=S_Nodes;

end

      %%
      Mesh_Information.Boundary_Cond=[1;1;1];
    
    
      [Mesh_Information]=Node_Rearranger(Element_Information,Mesh_Information,Mesh_Information.S_Nodes(:,1)); % Outer layer goes 1st
      [Mesh_Information]=Node_Rearranger(Element_Information,Mesh_Information,Mesh_Information.S_Nodes(:,2));
      [Mesh_Information]=Node_Rearranger(Element_Information,Mesh_Information,Mesh_Information.S_Nodes(:,3));
      

      Mesh_Information.BC_Nodes=Mesh_Information.S_Nodes(:,3);


%% Pair selection for boundary      
  Pair=[0 0];
  count=1;
  for i=1:length(Mesh_Information.S_Nodes(:,1))
      n1=Mesh_Information.S_Nodes(i,1);

      if n1>0

      Node1=Mesh_Information.Nodes(n1,2:3);
      for j=1:length(Mesh_Information.Nodes(:,1))
          n2=Mesh_Information.Nodes(j,1);
          Node2=Mesh_Information.Nodes(n2,2:3);
          if norm(Node1-Node2)<10^-2 && n1~=n2 && sum(abs(Node1-[0 0]))>10^-10
              Pair(count,1:2)=[n1 n2];
              count=count+1;
          end

      end
      end
  end

%%

Pair_for_Interpolate=zeros(1,3);
IEN=Mesh_Information.IEN;  
count=1;
for e=1:Mesh_Information.nel_Model

    nnp_el=Element_Information.nnp_Elemental(Mesh_Information.Type(1,e));
     
    nn=IEN(05,e);
    if sum(find(Pair_for_Interpolate(:,1)==nn))==0
    Pair_for_Interpolate(count,1:3)=[IEN(05,e) IEN(17,e) IEN(01,e)];
    count=count+1;
    end

    nn=IEN(13,e);
    if sum(find(Pair_for_Interpolate(:,1)==nn))==0
    Pair_for_Interpolate(count,1:3)=[IEN(13,e) IEN(26,e) IEN(09,e)];
    count=count+1;
    end

    nn=IEN(06,e);
    if sum(find(Pair_for_Interpolate(:,1)==nn))==0
    Pair_for_Interpolate(count,1:3)=[IEN(06,e) IEN(18,e) IEN(02,e)];
    count=count+1;
    end

    nn=IEN(14,e);
    if sum(find(Pair_for_Interpolate(:,1)==nn))==0
    Pair_for_Interpolate(count,1:3)=[IEN(14,e) IEN(21,e) IEN(10,e)];
    count=count+1;
    end

    nn=IEN(07,e);
    if sum(find(Pair_for_Interpolate(:,1)==nn))==0
    Pair_for_Interpolate(count,1:3)=[IEN(07,e) IEN(19,e) IEN(03,e)];
    count=count+1;
    end

    nn=IEN(15,e);
    if sum(find(Pair_for_Interpolate(:,1)==nn))==0
    Pair_for_Interpolate(count,1:3)=[IEN(15,e) IEN(22,e) IEN(11,e)];
    count=count+1;
    end

    nn=IEN(08,e);
    if sum(find(Pair_for_Interpolate(:,1)==nn))==0
    Pair_for_Interpolate(count,1:3)=[IEN(08,e) IEN(20,e) IEN(04,e)];
    count=count+1;
    end

    nn=IEN(16,e);
    if sum(find(Pair_for_Interpolate(:,1)==nn))==0
    Pair_for_Interpolate(count,1:3)=[IEN(16,e) IEN(27,e) IEN(12,e)];
    count=count+1;
    end

    nn=IEN(23,e);
    if sum(find(Pair_for_Interpolate(:,1)==nn))==0
    Pair_for_Interpolate(count,1:3)=[IEN(23,e) IEN(24,e) IEN(25,e)];
    count=count+1;
    end

end


Nodes=Mesh_Information.Nodes;
for nn=1:length(Pair_for_Interpolate)
   
    nnn=Pair_for_Interpolate(nn,1);
    if Nodes(nnn,4)>-10^-2

        Node_no_upper=Pair_for_Interpolate(nn,1);
        Rad_Org=norm(Nodes(Node_no_upper,2:4));
        Rad_req=Outer_Radius;
        Nodes(Node_no_upper,2:4)=Nodes(Node_no_upper,2:4)*Rad_req/Rad_Org;
    
    
        Node_no_lower=Pair_for_Interpolate(nn,3);
        Rad_Org=norm(Nodes(Node_no_lower,2:4));
        Rad_req=Inner_Radius;
        Nodes(Node_no_lower,2:4)=Nodes(Node_no_lower,2:4)*Rad_req/Rad_Org;

    else

        Node_no_upper=Pair_for_Interpolate(nn,1);
        Rad_Org=norm(Nodes(Node_no_upper,2:3));
        Rad_req=Outer_Radius;
        Nodes(Node_no_upper,2:3)=Nodes(Node_no_upper,2:3)*Rad_req/Rad_Org;
    
    
        Node_no_lower=Pair_for_Interpolate(nn,3);
        Rad_Org=norm(Nodes(Node_no_lower,2:3));
        Rad_req=Inner_Radius;
        Nodes(Node_no_lower,2:3)=Nodes(Node_no_lower,2:3)*Rad_req/Rad_Org;
    end
    
    
end

[Nodes]=Middle_Node_Account(Nodes,Pair_for_Interpolate);


Mesh_Information.Nodes=Nodes;

%% Processing      
Nodes=Mesh_Information.Nodes;
IEN=Mesh_Information.IEN;

try
for n=1:Mesh_Information.S_Nodes(:,1)
    if Mesh_Information.S_Nodes(n,1)==0
        S_Node_Length=n-1;
        break
    end
    S_Node_Length=n;
end
end
Node_Update_List=Mesh_Information.S_Nodes(1:S_Node_Length,1); 

%% Pair for update
count=1;
for i=1:length(Mesh_Information.S_Nodes(:,1))
   
    nn=Mesh_Information.S_Nodes(i,1);

    if nn>0
    
    x=Nodes(nn,2);
    y=Nodes(nn,3);
    z=Nodes(nn,4); 
    
    %hold on
    %scatter3(Nodes(nn,2),Nodes(nn,3),Nodes(nn,4),'filled')
    
    if y>10^-3
       for j=1:length(Mesh_Information.S_Nodes(:,1))
            nn2=Mesh_Information.S_Nodes(j,1); 

            if nn2>0
                x2=Nodes(nn2,2);
                y2=Nodes(nn2,3);
                z2=Nodes(nn2,4);
                
                if y2<0
                    
    %                 scatter3(Nodes(nn2,2),Nodes(nn2,3),Nodes(nn2,4),'filled')
    %                 pause
                    if abs(x2-x)<10^-5 && abs(z2-z)<10^-5
                       Node_Update_List2(count,1:2)=[nn nn2];
                       count=count+1;
                    end
                end
            end
       end
    elseif abs(y)<10^-5
           Node_Update_List2(count,1:2)=[nn 0];
           count=count+1;    
    end

    end
        
end
    
Node_Update_List2=sortrows(Node_Update_List2,2);

NC=Nodes(Node_Update_List2(:,1),2:4);
x=NC(:,1); y=NC(:,2); z=NC(:,3);
% scatter3(x,y,z)   



%% Runner
TS=500;

Mesh_Information_Org=Mesh_Information;


Nodes_Store=zeros(size(Mesh_Information.Nodes,1),size(Mesh_Information.Nodes,2),TS);
Temp_Store=zeros(size(Mesh_Information.Nodes,1),TS);

Time_Store=zeros(1,TS);
Error_Store=Time_Store;

NE=1000;

itertion=1;
TimeOut=0;
count2=1;

%%


Delta=1;

Check_rate=CHECK_RATE;

Lamda=LAMDA;

DelX=DELX;

Time_to_avoid=zeros(1,10)-100;
Avoid_Timeout=1;

% load Optimization_10_2_Quadratic_Tmax313_Tmin200_TInt100.mat

% pool = parpool('Threads');
           
%%



while NE>0
%----------------------

fprintf('------------------------------------\n');

% Set material


Mesh_Information_In=Mesh_Information_Org;


% if mod(itertion,10)==0
%     Nodes=Mesh_Information_Org.Nodes;
%     hold on
%     scatter3(Nodes(:,2),Nodes(:,3),Nodes(:,4))
%     daspect([1 1 1])
%     view(45,45)
%     axis tight
%     pause (0.01)
% end


T_One_Lunar_Day_Time=3600*24*15;

Final_Time=T_One_Lunar_Day_Time*2;

TimeOut_Old=TimeOut;

Initial_Time=-T_One_Lunar_Day_Time/4;
T_Initial=zeros(Mesh_Information_Org.neq_Model,1);

% if itertion>1
% 
% if max(time_Store)>-Initial_Time
%     Initial_Time=TimeOut-T_One_Lunar_Day_Time/4;
%     
%     Loc=find(time_Store==Initial_Time)
%     if length(Loc)==0
%         Loc=find(round(time_Store)==Initial_Time)
%         fprintf("Warning! Exact time point can't be detected for original simulation\n");
%     end
%     T_Initial=T_Original(:,Loc(1));
%     pause
% end
% end

if itertion==1
    Error_Previous=100;
    [Normal_Vector_form,TimeOut,T_Original,time_Store,T,Error]=Thermal_Model_Runner_Complete(Element_Information,Mesh_Information_In,Material_Information,T_One_Lunar_Day_Time,Initial_Time,Final_Time,T_Initial,Alpha,Maximum_allowed_temperature,T_ground,Check_rate,Time_to_avoid);
else
    Error_Previous=Error_New;
    [~,TimeOut,T_Original,time_Store,T,Error]=Thermal_Model_Runner_Complete(Element_Information,Mesh_Information_In,Material_Information,T_One_Lunar_Day_Time,Initial_Time,Final_Time,T_Initial,Alpha,Maximum_allowed_temperature,T_ground,Check_rate,Time_to_avoid);
end

Error_New=norm(Error);





% if TimeOut_Old>TimeOut && itertion>2
% 
%     Mesh_Information_Org.Nodes=Nodes_Store(:,:,count2-2);
%     % T=T_Store(:,count2-2);
% 
%     TimeOut=TimeOut_Old;
%     Error=Error_Store(itertion-2);
%     DelX=DelX/2;
% 
% end

% hold on
% view(45,45)
% daspect([1 1 1])
% scatter3(Nodes(:,2),Nodes(:,3),Nodes(:,4),'filled')
fprintf('\nIteration=%2d   Error=%6f   Time max=%3d hr   Max temp=%3d(304) Min temp=%3d(229) Check rate=%d Lambda=%f delX=%f\n',itertion,Error_New,TimeOut/(60*60),int16(max(T)),int16(min(T)),Check_rate/(60*60),Lamda,DelX)
                        

Nodes_Store(:,:,count2)=Mesh_Information_Org.Nodes; 
Temp_Store(:,count2)=T;
count2=count2+1;

Cost_Func=sum(abs(Error));

dCost_Func_Forward=eye(length(Node_Update_List));

dCost_Func_Backwards=eye(length(Node_Update_List));

%% ----------------------------

Initial_Time=TimeOut-40*60*60;
if Initial_Time<0
   Initial_Time=0;
end

Loc=find(time_Store==Initial_Time);
if length(Loc)==0
    Loc=find(round(time_Store)==Initial_Time);
    fprintf("Warning! Exact time point can't be detected\n");
end
T_Initial=T_Original(:,Loc(1));

Final_Time=TimeOut;

Nodes_Org=Mesh_Information_Org.Nodes;


%%
fprintf("Getting cost function derivetive using Forward difference...")
parfor (i=1:length(Node_Update_List),4)

    if abs(Error(i))>0
    
    Nodes=Nodes_Org;
    
    nn=Node_Update_List(i,1);
    
    if nn>0
    
    Vector=Normal_Vector_form(3*nn-2:3*nn);
    Nodes(nn,2:4)=Nodes(nn,2:4)+Vector'*DelX;

    % Account for base
    for pp=1:length(Pair)
       if Pair(pp,1)==nn
          nn_Pair=Pair(pp,2);
          VectorXY=Vector; VectorXY(3)=0;
          Nodes(nn_Pair,2:4)=Nodes(nn_Pair,2:4)+VectorXY'*DelX;
       end
    end  

    % Account for middle nodes
    [Nodes]=Middle_Node_Account(Nodes,Pair_for_Interpolate);
    
    Mesh_Information_Up=Mesh_Information_Org;
    Mesh_Information_Up.Nodes=Nodes;
    
    Mesh_Information_In=Mesh_Information_Up;
    
    % [~,~,~,Error_i]=Thermal_Model_Runner(Element_Information,Mesh_Information_In,Material_Information,T_One_Lunar_Day_Time,Final_Time);
    
    % [~,~,~,~,T_i,Error_i]=Thermal_Model_Runner_Complete(Element_Information,Mesh_Information_In,Material_Information,T_One_Lunar_Day_Time,Initial_Time,Final_Time,T_Initial,Alpha);

    [~,~,~,~,T_i,Error_i]=Thermal_Model_Runner_Complete(Element_Information,Mesh_Information_In,Material_Information,T_One_Lunar_Day_Time,Initial_Time,Final_Time,T_Initial,Alpha,Maximum_allowed_temperature,T_ground,Check_rate,Time_to_avoid);


    dCost_Func_Forward(:,i)=Error_i-Error;

%     Cost_Func_i=sum(abs(Error_i));
%     dCost_Func_Forward(i,1)=(Cost_Func_i-Cost_Func);

    end

    end
    
end

% fprintf("using Backward difference...")
% parfor (i=1:length(Node_Update_List),4)
%     
%     if abs(Error(i))>0
% 
%     Nodes=Nodes_Org;
%     
%     nn=Node_Update_List(i,1);
%     
%     if nn>0
%     
%     Vector=Normal_Vector_form(3*nn-2:3*nn);
%     Nodes(nn,2:4)=Nodes(nn,2:4)-Vector'*DelX;
% 
%     % Account for base
%     for pp=1:length(Pair)
%        if Pair(pp,1)==nn
%           nn_Pair=Pair(pp,2);
%           VectorXY=Vector; VectorXY(3)=0;
%           Nodes(nn_Pair,2:4)=Nodes(nn_Pair,2:4)-VectorXY'*DelX;
%        end
%     end  
% 
%     % Account for middle nodes
%     [Nodes]=Middle_Node_Account(Nodes,Pair_for_Interpolate);
%     
%     Mesh_Information_Up=Mesh_Information_Org;
%     Mesh_Information_Up.Nodes=Nodes;
%     
%     Mesh_Information_In=Mesh_Information_Up;
%     
%     % [~,~,~,Error_i]=Thermal_Model_Runner(Element_Information,Mesh_Information_In,Material_Information,T_One_Lunar_Day_Time,Final_Time);
%     
%     % [~,~,~,~,T_i,Error_i]=Thermal_Model_Runner_Complete(Element_Information,Mesh_Information_In,Material_Information,T_One_Lunar_Day_Time,Initial_Time,Final_Time,T_Initial,Alpha);
% 
%     [~,~,~,~,T_i,Error_i]=Thermal_Model_Runner_Complete(Element_Information,Mesh_Information_In,Material_Information,T_One_Lunar_Day_Time,Initial_Time,Final_Time,T_Initial,Alpha,Maximum_allowed_temperature,T_ground,Check_rate);
% 
%     dCost_Func_Backwards(:,i)=Error_i-Error;
% 
%     %Cost_Func_i=sum(abs(Error_i));
%     %dCost_Func_Backwards(i,1)=(Cost_Func_i-Cost_Func);
% 
%     end
% 
%     end
%     
% end


%% ------------------------

    % dF=0.5*(-dCost_Func_Backwards+dCost_Func_Forward)/DelX;

    dF=(dCost_Func_Forward)/DelX;


for i=1:length(Node_Update_List)
    if abs(Error(i))==0
        dF(i,i)=1;
    end
end

if rank(dF)<size(dF)
    fprintf('Linear dependence detected!')
    dF=dF.*eye(size(dF));
end


% dS=(dCost_Func_Backwards-2*Cost_Func+dCost_Func_Forward)/DelX^2;

N_Move=-inv(dF)*Error*Lamda;
                    
if max(N_Move)==0
    fprintf('No movement at all! Avoiding this time all together!')

    Time_to_avoid(Avoid_Timeout)=TimeOut;
    Avoid_Timeout=Avoid_Timeout+1;

    % DelX=DelX/2;
else
%     DelX=DelX*2;
%     if DelX>DELX
%         DelX=DELX;
%     end
end

if max(N_Move)>DelX*100

    fprintf('Too big movement, %d (m) Correcting!',max(N_Move))
    N_Move=N_Move/max(N_Move)*(DelX*100);
    

end

% N_Move=-(Cost_Func)./(dCost_Func_Forward);

% if max(N_Move)>DelX
% 
%     Rat=max(N_Move)/DelX;
%     N_Move=N_Move/Rat/2;
% 
% end

% Update the nodal coordinates------------------------------

Nodes=Mesh_Information_Org.Nodes;

for i=1:length(Node_Update_List)
    
    nn=Node_Update_List(i,1);
    
    Vector=Normal_Vector_form(3*nn-2:3*nn);
    
    if N_Move(i)>0 
    else
%         if TimeOut<TimeOut_Old
%         else
            N_Move(i)=N_Move(i)*0.0;
%        end
    end
    
    Nodes_Cq=Nodes(nn,2:4)+Vector'*N_Move(i);
    % Node_Prev_rad=norm(Mesh_Information.Nodes(nn,2:4));
    
        if norm(Nodes_Cq)>Inner_Radius
            Nodes(nn,2:4)=Nodes(nn,2:4)+Vector'*N_Move(i);
            
            % Account for base
            for pp=1:length(Pair)
               if Pair(pp,1)==nn
                  nn_Pair=Pair(pp,2);
                  VectorXY=Vector; VectorXY(3)=0;
                  Nodes(nn_Pair,2:4)=Nodes(nn_Pair,2:4)+VectorXY'*N_Move(i);
               end
            end            
             
            % Account for middle nodes
            [Nodes]=Middle_Node_Account(Nodes,Pair_for_Interpolate);

        end

end

fprintf('Updating model parameters...')

UpdateCheckRate=1;
if abs(Error_Previous-Error_New)<10^-10
    fprintf ('\nNot converging, updating damping\n')

    % Check_rate=Check_rate*2;
    % UpdateCheckRate=0;
    
%     if Check_rate>40*(60*60)
%        fprintf ('\nUnable to converge!\n')
%        TS=itertion;
%     end

    Lamda=Lamda/2;
    if Lamda<0.1
        fprintf ('\nUnable to converge!\n')
        TS=itertion;
    end

else
    Lamda=Lamda*2;
    if Lamda>LAMDA
        Lamda=LAMDA;
    end
end


try

    Mesh_Information_Test=Mesh_Information_Org;
    Mesh_Information_Test.Nodes=Nodes;

    [FE_Information_Region]=Region_Integration_3D_Heat_Transfer(Element_Information,Mesh_Information_Test);
    [FE_Information_Surface]=Surface_Integration_3D_Flux_Thermal(Element_Information,Mesh_Information_Test);
    [FE_Information_Surface_Solid]=Surface_Integration_3D_Pressure_Type_Solid(Element_Information,Mesh_Information_Test);

    if UpdateCheckRate==1
        Check_rate=Check_rate*2;
        if Check_rate>CHECK_RATE
            Check_rate=CHECK_RATE;
        end
    end

    Mesh_Information_Org.Nodes=Nodes;
    fprintf('Done \n')

    % NE_Store(itertion)=NE;
    Time_Store(itertion)=TimeOut;
    Error_Store(itertion)=norm(Error);
    
    NE=Cost_Func;
    
    
    % Nodes_Store(:,itertion)=Nodes;
    
    % fprintf('--------------------------------------------------------------------------------------------\n')
    % fprintf('--------------------------------------------------------------------------------------------\n')
    % fprintf('--------------------------------------------------------------------------------------------\n')
    itertion=itertion+1;
    % fprintf('--------------------------------------------------------------------------------------------\n')
    % fprintf('--------------------------------------------------------------------------------------------\n')
    % fprintf('--------------------------------------------------------------------------------------------\n')
   
    if mod(itertion,2)==0
        save Optimization_Safety_Even
    else
        save Optimization_Safety_Odd
    end

catch

    fprintf ('\nNew geometry is invalid, updating check rate and lambda\n')
    Check_rate=Check_rate/2;
    Lamda=Lamda/2;

end




%---------------------------------------------------

% if mod(itertion,10)==0
%     scatter3(Nodes(:,2),Nodes(:,3),Nodes(:,4),10)
%     daspect([1 1 1])
%     axis tight
%     zlim([0 4])
%     pause(0.01)
% end

    if itertion>TS || TimeOut>T_One_Lunar_Day_Time*1.5
        break
    end

end


fprintf('Optimization complete')
Nodes=Mesh_Information_Org.Nodes;
disp(Nodes)


save 11_23_2024_

% if Alpha==0.3927
%     save 11_20_24_TMax50_225
% elseif Alpha==0.7854
%     save 11_20_24_TMax50_45
% elseif Alpha==1.2654
%     save 11_20_24_TMax50_72
% end

if AlphaCount==0
    save Data_0_5
elseif AlphaCount==1
    save Data_1_5
elseif AlphaCount==2
    save Data_2_5
elseif AlphaCount==3 
    save Data_3_5
elseif AlphaCount==4 
    save Data_4_5
elseif AlphaCount==5
    save Data_5_5
elseif AlphaCount==6
    save Data_6_5
elseif AlphaCount==7
    save Data_7_5
elseif AlphaCount==8
    save Data_8_5
elseif AlphaCount==9
    save Data_9_5
end


AlphaCount=AlphaCount+1;

    end
