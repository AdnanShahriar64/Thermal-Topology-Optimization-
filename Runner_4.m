clear
clc

addpath('Element_Library')
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

function[Nodes]=Middle_Node_Account(Nodes,Pair_for_Interpolate)


for n=1:size(Pair_for_Interpolate,1)

    Upper_Node=Pair_for_Interpolate(n,1);
    Middle_Node=Pair_for_Interpolate(n,2);
    Lower_Node=Pair_for_Interpolate(n,3);

    Node_Upper_xyz=Nodes(Upper_Node,2:4);

    Node_Lower_xyz=Nodes(Lower_Node,2:4);

    Node_Middle_xyz=(Node_Upper_xyz+Node_Lower_xyz)/2;

    Nodes(Middle_Node,2:4)=Node_Middle_xyz;

end



end




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



function[FE_Information_Region]=Region_Integration_3D_Heat_Transfer(Element_Information,Mesh_Information)
%%
Line='-------------------------------------------------------------------------------------';
% fprintf('\n%s\n',Line)
% fprintf('%s\n','Calculating Region Integration parameters...');
% fprintf('%s\n',Line)
% fprintf('%15s   %25s    %25s\n','Element no','Integration point no','Jacobian determinant') 
% fprintf('%s\n',Line)

Nodes=Mesh_Information.Nodes;
IEN=Mesh_Information.IEN;
Type=Mesh_Information.Type;
DShape_Elemental=Element_Information.DShape_Elemental;
W_Elemental=Element_Information.W_Elemental;
Shape_Elemental=Element_Information.Shape_Elemental;
nel 	= Mesh_Information.nel_Model;

%% B matrix structure
%   [dN1/dx 0       0       dN2/dx 0      0      0...|
%    0      dN1/dy  0       0      dN1/dy 0      0...| Integration point 1
%    0      0       dN1/dz  0      0      dN2/dz 0...|

%    dN1/dx 0       0       dN2/dx 0      0      0...|
%    0      dN1/dy  0       0      dN1/dy 0      0...| Integration point 2
%    0      0       dN1/dz  0      0      dN2/dz 0...|

%% Memory Preallocation
    B=zeros(length(W_Elemental(:,end))*3,Mesh_Information.nnp_max_Model,nel);
    dJB=zeros(length(W_Elemental(:,end)),nel);    

%% Each element calculation

for e=1:nel

        nn=Element_Information.nnp_Elemental(Type(1,e));                                                         % No of Nodes
        xyz=[Nodes(IEN(1:nn,e),2) Nodes(IEN(1:nn,e),3) Nodes(IEN(1:nn,e),4)];   % Setting up xyz vector

        
        for gp=1:Element_Information.ngp_Elemental(Type(1,e))                            % Goes through all Gauss points
                % Determine B matrix ---------------------------------
                % Jacobian Matrix
                    J=DShape_Elemental(:,1:Type(2,e),gp,Type(1,e))*xyz;
                % Determinant of Jacobian Matrix
                    dJ=det(J);
                    if dJ<=0
                        error('Warning! Mesh is invalid, negative determinant of jacobian detected, element no=%d, integration point no=%d\n',e,gp);
                    end
                % Derivative matrix in global coordinate
                    CDsx3=inv(J)*DShape_Elemental(:,1:Type(2,e),gp,Type(1,e));
                % B Matrix
                    %    Be=zeros(3,nn);
                    %for k=1:Type(2,e)
                    %    Be(:,3*k-2:3*k)= [CDsx3(1,k)       0           0      ;
                    %                          0       CDsx3(2,k)       0      ;
                    %                          0             0       CDsx3(3,k)];
                    %end
                    
                    
                    if sum(abs(Shape_Elemental(gp,:,Type(1,e))))<10^-3 || dJ<=0 || sum(find(isnan(CDsx3)))>0
                        fprintf(2,'Warning! Weight is set to 0, element no=%d, integration point no=%d\n',e,gp);
                    else
                    
                        B(gp*3-2:gp*3,1:nn,e)=CDsx3;
                        dJB(gp,e)=dJ;
                    end
                    
                   % fprintf('%15d   %25d    %25f\n',e,gp,dJ)
                    
        end

end

FE_Information_Region.B=B;
FE_Information_Region.dJB=dJB;

% fprintf('%s\n\n',Line)

end
        

function[FE_Information_Surface]=Surface_Integration_3D_Flux_Thermal(Element_Information,Mesh_Information)
%%
Line='---------------------------------------------------------------------------------------------';
% fprintf('\n%s\n',Line)
% fprintf('%s\n','Calculating Surface Integration parameters');
% fprintf('%s\n',Line)

% fprintf('%15s   %15s    %25s   %25s\n','Element no','Local Surface no','Integration point no','Jacobian determinant') 
% fprintf('%s\n',Line)



S_Type=Mesh_Information.S_Type;
SS=size(S_Type);
Nodes=Mesh_Information.Nodes;
IEN=Mesh_Information.IEN;
Type=Mesh_Information.Type;
Surf_Dshape=Element_Information.Surf_Dshape_Elemental;
Surf_W=Element_Information.Surf_W_Elemental;
Surf_Shape=Element_Information.Surf_Shape_Elemental;
ndof=Mesh_Information.ndof_Model;

SNs=max(S_Type(1,:));
T_surf_gp=zeros(Mesh_Information.nnp_max_Model,Mesh_Information.ngp_surf_max_Model,SS(2));
dTeS_gp=zeros(3,Mesh_Information.ngp_surf_max_Model,SS(2));

T_surf=zeros(Mesh_Information.nnp_max_Model,max(Mesh_Information.nsurf_Model),SNs);
Ke_surf=zeros(Mesh_Information.nnp_max_Model,Mesh_Information.nnp_max_Model,Mesh_Information.nel_Model,SNs);
dTeS=zeros(3,max(Mesh_Information.nsurf_Model),SNs);


%%
SNo=0;
Scount=0;
for i=1:SS(2)

        SN=S_Type(1,i);     % Extract Surface Number of the Model

        e=S_Type(2,i);      % Extract Element no
        nn=Element_Information.nnp_Elemental(Type(1,e));       % No of Nodes in the elements
        xyz=[Nodes(IEN(1:nn,e),2) Nodes(IEN(1:nn,e),3) Nodes(IEN(1:nn,e),4)];   % Extract xyz
        Surf=S_Type(3,i);    % Extract Which surface of the element

        if SN>SNo
            Scount=0;
        end
            Scount=Scount+1;


        Surf_gp_for_Surf=Element_Information.Surf_ngp_Elemental(Surf,Type(1,e));
        

        Par_Coordinate=Element_Information.Par_Coord_Elemental(:,:,Surf,Type(1,e));

            
        dTea=[0 0 0]';
        for gp=1:Surf_gp_for_Surf
           
            J=Surf_Dshape(:,1:nn,gp,Surf,Type(1,e))*xyz;
            
            if det(J)<=0
                fprintf(2,'Warning! Mesh is invalid, negative determinant of jacobian on surface detected, Element no=%d, Surface no=%d, Integration point no=%d\n',e,Surf,gp);
            end
            
            Vectors=inv(J)*Par_Coordinate;

            dTe=(cross(Vectors(:,1),Vectors(:,2)));
            dTe=dTe/(norm(dTe))^2; 
            
            T_surf(1:nn,Scount,SN)=T_surf(1:nn,Scount,SN)+Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*norm(dTe);
            
            T_surf_gp(1:nn,gp,i)=Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*norm(dTe);


            Ke_surf(1:nn,1:nn,e,SN)=Ke_surf(1:nn,1:nn,e,SN)+(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))*Surf_W(gp,Surf,Type(1,e))*norm(dTe);
            
            dTea=dTea+dTe;
            
            dTeS_gp(:,gp,i)=dTe;

            % fprintf('%15d   %15d    %25d   %25f\n',e,SN,gp,det(J))

        end
            
        dTeS(:,Scount,SN)=dTea/norm(dTea);
        
            
        %% https://scicomp.stackexchange.com/questions/26886/evaluating-the-surface-integral-in-an-fem-finite-elements-method-procedure

        SNo=SN;

end

FE_Information_Surface.T_surf_gp=T_surf_gp;
FE_Information_Surface.dTeS_gp=dTeS_gp;
FE_Information_Surface.dTeS=dTeS;
FE_Information_Surface.Ke_surf=Ke_surf;
FE_Information_Surface.T_surf=T_surf;


% fprintf('%s\n\n',Line)

end

function[Dm]=Material_Array_Generator_Heat_Transfer(Material_Information)

Line='---------------------------------------------------------------------------------------------';
% fprintf('\n%s\n',Line)
% fprintf('%s','6. Calculating material matrix...');

Kx=Material_Information.Kx;

%% Material Matrix initialize
    Dm=zeros(3,3,length(Kx));
    
%% generate material matrix
for i=1:length(Kx)
    

    Dm=[Kx  0   0
        0   Kx  0
        0   0   Kx];
 
 
end

%%

% fprintf('%s\n','Done');
% fprintf('%s\n',Line)

end

function[FE_Information_Surface]=Surface_Integration_3D_Pressure_Type_Solid(Element_Information,Mesh_Information)
%%
% Line='---------------------------------------------------------------------------------------------';
% fprintf('\n%s\n',Line)
% fprintf('%s\n','Calculating Surface Integration parameters');
% fprintf('%s\n',Line)
% fprintf('%15s   %15s    %25s   %25s\n','Element no','Local Surface no','Integration point no','Jacobian determinant') 
% fprintf('%s\n',Line)



S_Type=Mesh_Information.S_Type;
SS=size(S_Type);
Nodes=Mesh_Information.Nodes;
IEN=Mesh_Information.IEN;
Type=Mesh_Information.Type;
Surf_Dshape=Element_Information.Surf_Dshape_Elemental;
Surf_W=Element_Information.Surf_W_Elemental;
Surf_Shape=Element_Information.Surf_Shape_Elemental;
ndof=Mesh_Information.ndof_Model;

SNs=max(S_Type(1,:));
% P_surf_gp=zeros(Mesh_Information.nnp_max_Model*Mesh_Information.ndof_Model,Mesh_Information.ngp_max_Model,max(Mesh_Information.nsurf_Model),SNs);
dTeS_gp=zeros(3,Mesh_Information.ngp_max_Model,max(Mesh_Information.nsurf_Model),SNs);

P_surf=zeros(Mesh_Information.nnp_max_Model*ndof,max(Mesh_Information.nsurf_Model),SNs);
dTeS=zeros(3,max(Mesh_Information.nsurf_Model),SNs);



%%
SNo=0;
Scount=0;
for i=1:SS(2)

        SN=S_Type(1,i);     % Extract Surface Number of the Model

        e=S_Type(2,i);      % Extract Element no
        nn=Element_Information.nnp_Elemental(Type(1,e));       % No of Nodes in the elements
        xyz=[Nodes(IEN(1:nn,e),2) Nodes(IEN(1:nn,e),3) Nodes(IEN(1:nn,e),4)];   % Extract xyz
        Surf=S_Type(3,i);    % Extract Which surface of the element

        if SN>SNo
            Scount=0;
        end
            Scount=Scount+1;


        Surf_gp_for_Surf=Element_Information.Surf_ngp_Elemental(Surf,Type(1,e));
        

        Par_Coordinate=Element_Information.Par_Coord_Elemental(:,:,Surf,Type(1,e));

            
        dTea=[0 0 0]';
        for gp=1:Surf_gp_for_Surf
           
            J=Surf_Dshape(:,1:nn,gp,Surf,Type(1,e))*xyz;
            
            if det(J)<=0
                fprintf(2,'Warning! Mesh is invalid, negative determinant of jacobian on surface detected, Element no=%d, Surface no=%d, Integration point no=%d\n',e,Surf,gp);
            end
            
            Vectors=(transpose(J))*Par_Coordinate;

            dTe=(cross(Vectors(:,1),Vectors(:,2)));
            % dTe=dTe/(norm(dTe))^2;
            
            
            % Valx=Valx+(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))*P;
            % hold on
            % plot(Valx)
            % pause
            % Valx=Valx+Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(1);
            % Valy=Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(2);
            % Valz=Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(3);
            % plot(Valx)
            % hold on
            
            P_surf(1:3:nn*ndof-2,Scount,SN)=P_surf(1:3:nn*ndof-2,Scount,SN)+Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(1);
            P_surf(2:3:nn*ndof-1,Scount,SN)=P_surf(2:3:nn*ndof-1,Scount,SN)+Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(2);
            P_surf(3:3:nn*ndof-0,Scount,SN)=P_surf(3:3:nn*ndof-0,Scount,SN)+Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(3);

            % P_surf_gp(1:3:nn*ndof-2,gp,Scount,SN)=Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(1);
            % P_surf_gp(2:3:nn*ndof-1,gp,Scount,SN)=Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(2);
            % P_surf_gp(3:3:nn*ndof-0,gp,Scount,SN)=Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(3);

            dTeS_gp(:,gp,Scount,SN)=dTe;

            dTea=dTea+dTe;
            
            % fprintf('%15d   %15d    %25d   %25f\n',e,SN,gp,det(J))
            

        end
            
        dTeS(:,Scount,SN)=dTea/norm(dTea);
        
            
        %% https://scicomp.stackexchange.com/questions/26886/evaluating-the-surface-integral-in-an-fem-finite-elements-method-procedure

        SNo=SN;

end

% FE_Information_Surface.P_surf_gp=P_surf_gp;

FE_Information_Surface.dTeS_gp=dTeS_gp;
FE_Information_Surface.dTeS=dTeS;
FE_Information_Surface.P_surf=P_surf;



% fprintf('%s\n\n',Line)

end


function [K] = Assembly(Ke,ndof,Mesh_Information,Element_Information)


Type=Mesh_Information.Type;
nel=Mesh_Information.nel_Model;
nnp=Mesh_Information.nnp_Model;
nnp_Elemental=Element_Information.nnp_Elemental;

if ndof==1
    LM=Mesh_Information.IEN;
elseif ndof==3
    LM=Mesh_Information.LM_Model;
else
    error('This assembly file is not built for any DOF except 1 and 3')
end


K 	= zeros(nnp*ndof);

%% Assembly of K and M
%fprintf('\nAssemblying...')

for e = 1:nel
    nn=nnp_Elemental(Type(1,e));
    for loop1 = 1:nn*ndof
        i = LM(loop1,e);
        for loop2 =  1:nn*ndof
            j = LM(loop2,e);
            K(i,j) = K(i,j) + Ke(loop1,loop2,e);
        end
    end
    %fprintf('Element=%d   Type=%d\n',e,Type(e));
end
%fprintf('Done\n')


end

function[Me,H_Source_per_H]=Elemental_Capacitance_Matrix_Generator(Element_Information,Mesh_Information,Material_Information,FE_Information_Region)

    SpcH=Material_Information.SpcH;
    rho=Material_Information.rho;
    B=FE_Information_Region.B;
    dJB=FE_Information_Region.dJB;
    W=Element_Information.W_Elemental;
    Shape=Element_Information.Shape_Elemental;
    Type=Mesh_Information.Type;
    nel=Mesh_Information.nel_Model;
    ndof=Mesh_Information.ndof_Model;
    neq=Mesh_Information.neq_Model;
    LM=Mesh_Information.LM_Model;
    
    nnp_Max=Mesh_Information.nnp_max_Model;
    ngp_Max=Mesh_Information.ngp_max_Model;
    ngp_Elemental=Element_Information.ngp_Elemental;
    nnp_Elemental=Element_Information.nnp_Elemental;


    Me=zeros(nnp_Max,nnp_Max,nel);   % Elemental Mass Matrix preallocation

    H_Source_per_H=zeros(nnp_Max,nnp_Max,nel);
    
    for e=1:nel
        
        gpLimit=ngp_Elemental(Type(1,e));
        nn=nnp_Elemental(Type(1,e));

        for gp=1:gpLimit
            
            N=Shape(gp,1:nn,Type(1,e));

            Me(1:nn,1:nn,e)=Me(1:nn,1:nn,e)+SpcH*rho*(N'*N)*dJB(gp,e)*W(gp,Type(1,e));
            
            H_Source_per_H(1:nn,gp,e)=N'*dJB(gp,e)*W(gp,Type(1,e));
            
        end
        
    end
    
end

function[Ke]=Elemental_Conductivity_Matrix_Generator(Element_Information,Mesh_Information,Material_Information,FE_Information_Region,Dm)
%%

% rho=Material_Information.rho;
B=FE_Information_Region.B;
dJB=FE_Information_Region.dJB;
W=Element_Information.W_Elemental;
% Shape=Element_Information.Shape_Elemental;
Type=Mesh_Information.Type;
nel=Mesh_Information.nel_Model;
% ndof=Mesh_Information.ndof_Model;
% neq=Mesh_Information.neq_Model;
% LM=Mesh_Information.LM_Model;

nnp_Max=Mesh_Information.nnp_max_Model;
ngp_Max=Mesh_Information.ngp_max_Model;
ngp_Elemental=Element_Information.ngp_Elemental;
nnp_Elemental=Element_Information.nnp_Elemental;

Ke=zeros(nnp_Max*1,nnp_Max*1,nel);   % Elemental Stiffness Matrix preallocation


for e=1:nel
    
    gpLimit=ngp_Elemental(Type(1,e));
    nn=nnp_Elemental(Type(1,e));

    for gp=1:gpLimit
        
        B_matrix=B(gp*3-2:gp*3,1:nn,e);
        
        Ke(1:nn,1:nn,e)=Ke(1:nn,1:nn,e)+B_matrix'*Dm*B_matrix*dJB(gp,e)*W(gp,Type(1,e));
      
    end

    
end


end


function[fS]=Pressure_to_Equivalent_Nodal_Force_Nodal_Pressure_Lag_GL_2(Element_Information,Mesh_Information,FE_Information_Surface,Nodal_Pressure)

%%
    SS=size(Mesh_Information.S_Type);

    ndof=Mesh_Information.ndof_Model;

    % SP=size(FE_Information_Surface.P_surf);

    fS=zeros(Mesh_Information.nnp_Model*ndof,1);

    count=1;

    SNO=1;
    for i=1:SS(2) % For each surfaces


            SN=Mesh_Information.S_Type(1,i);                            % Extract Surface no 
            if SN>SNO
                count=1;
            end

            Surf_Type=Mesh_Information.S_Type(3,i);
        
            e=Mesh_Information.S_Type(2,i);                         % Extract element no

            nn=Element_Information.nnp_Elemental(Mesh_Information.Type(1,e));

            ng=Element_Information.ngp_Elemental(Mesh_Information.Type(1,e));

            if nn~=ng || max(max(Element_Information.Surf_Shape_Elemental(:,:,Surf_Type,Mesh_Information.Type(1,e))))~=1
                error('Change the Pressure calculation file, this doens seem like a Lagrangian element with GL nodes')
            end

            % Pes=Ps(i);                                       % Extract Element Pressure

            Pes=Nodal_Pressure(Mesh_Information.IEN(1:nn,e));
            
            PesV=zeros(nn*3,1);
            PesV(1:3:end-2)=Pes;
            PesV(2:3:end-1)=Pes;
            PesV(3:3:end-0)=Pes;
            
            fe=FE_Information_Surface.P_surf(1:nn*ndof,count,SN).*PesV;                  % Determine Equivalent Nodal FV

            fS(Mesh_Information.LM_Model(1:nn*3,e))=fS(Mesh_Information.LM_Model(1:nn*3,e))+fe;     % Put into Total Force Vector

            SNO=SN;

            count=count+1;

    end
   

end
            

      
function[QS]=Flux_to_Equivalent_Nodal_Flux_Lag_GL(Element_Information,Mesh_Information,FE_Information_Surface,Nodal_Flux)

%%
    SS=size(Mesh_Information.S_Type);

    ndof=1;

    % SP=size(FE_Information_Surface.P_surf);

    QS=zeros(Mesh_Information.nnp_Model*ndof,1);

    count=1;

    SNO=1;
    for i=1:SS(2) % For each surfaces


            SN=Mesh_Information.S_Type(1,i);                            % Extract Surface no 
            if SN>SNO
                count=1;
            end

            Surf_Type=Mesh_Information.S_Type(3,i);
        
            e=Mesh_Information.S_Type(2,i);                         % Extract element no

            nn=Element_Information.nnp_Elemental(Mesh_Information.Type(1,e));

            ng=Element_Information.ngp_Elemental(Mesh_Information.Type(1,e));

            if nn~=ng || max(max(Element_Information.Surf_Shape_Elemental(:,:,Surf_Type,Mesh_Information.Type(1,e))))~=1
                error('Change the Pressure calculation file, this doens seem like a Lagrangian element with GL nodes')
            end

            % Pes=Ps(i);                                       % Extract Element Pressure
            
            QV=Nodal_Flux(Mesh_Information.IEN(1:nn,e));

            % QV=zeros(nn,1);
            % PesV(1:3:end-2)=Q_nn;
            
            Qe=FE_Information_Surface.T_surf(1:nn*ndof,count,SN).*QV;                  % Determine Equivalent Nodal FV

            QS(Mesh_Information.IEN(1:nn,e))=QS(Mesh_Information.IEN(1:nn,e))+Qe;     % Put into Total Force Vector

            SNO=SN;

            count=count+1;

    end
   

end
            

function[Mesh_Information]=Geometry_Model_Parameters_Preprocessor_Solid(Mesh_Information,Element_Information)

Line='----------------------------------------------------';
% fprintf('\n%s\n',Line)
% fprintf('%s\n','Calculating basic model parameters:');
% fprintf('%s\n',Line)


%% Initialize System  
    nsd     = 3;                        % number of space dimensions 
    ndof_Model 	= 3;                    % number of degrees-of-freedom per node
     
    nnp_Model 	= Mesh_Information.Nodes(end,1);     	% number of nodal points
    S=size(Mesh_Information.IEN);
    nel_Model 	  = S(2);                 % number of elements
    neq_Model     = ndof_Model*nnp_Model;             % number of equations    
    
%% generate LM array 
    LM_Model=zeros(max(Mesh_Information.Type(2,:))*ndof_Model,nel_Model);
    for e = 1:nel_Model
        
        nen=Mesh_Information.Type(2,e)+Mesh_Information.Type(5,e);

        for j = 1:nen
            for m = 1:ndof_Model
                ind = (j-1)*ndof_Model + m;
                 LM_Model(ind,e) = ndof_Model*Mesh_Information.IEN(j,e) - ndof_Model + m;
            end
        end
    end   

%% Generate Predefined Model Matrix
    % [IEN_Colm,IEN_Colm2]=Predefined_Matrix_Assembly(nel,Type,IEN);
    
    


%% Max properties
nn_e_O=0;
for e=1:size(Mesh_Information.Type,2)
    nn_e=Element_Information.nnp_Elemental(Mesh_Information.Type(1,e));
    if nn_e>nn_e_O
        nn_e_O=nn_e;
    end
end
nnp_max_Model=nn_e_O;

ngp_e_O=0;
for e=1:size(Mesh_Information.Type,2)
    ngn_e=Element_Information.ngp_Elemental(Mesh_Information.Type(1,e));
    if ngn_e>ngp_e_O
        ngp_e_O=ngn_e;
    end
end
ngp_max_Model=ngp_e_O;


nsurf_Model=ones(1,max(Mesh_Information.S_Type(1,:)));
nsurf_Model(1,1)=size(Mesh_Information.S_Type,2);
s_prev=1;
for ss=2:size(Mesh_Information.S_Type,2)
    if Mesh_Information.S_Type(1,ss)>Mesh_Information.S_Type(1,ss-1)
        nsurf_Model(1,Mesh_Information.S_Type(ss-1))=ss-s_prev;
        s_prev=ss;
    end
end    
if length(nsurf_Model)>1
    nsurf_Model(1,Mesh_Information.S_Type(ss))=ss-s_prev+1;
end

gp_Old=0;
for ss=2:size(Mesh_Information.S_Type,2)
    eT=Mesh_Information.Type(1,Mesh_Information.S_Type(2,ss));
    sT=Mesh_Information.S_Type(3,ss);
    gp=Element_Information.Surf_ngp_Elemental(sT,eT);
    if gp>gp_Old
        gp_Old=gp;
    end
end

ngp_surf_max_Model=gp_Old;



% [gp_X,gp_Y,gp_Z]=Gauss_Point_locations_3D(Mesh_Information.Nodes,Mesh_Information.IEN,Mesh_Information.Type,Element_Information.Shape_Elemental);

[gp_X,gp_Y,gp_Z]=Gauss_Point_locations_3D(Element_Information,Mesh_Information);

Mesh_Information.ndof_Model=ndof_Model;
Mesh_Information.nnp_Model=nnp_Model;
Mesh_Information.nel_Model=nel_Model;
Mesh_Information.neq_Model=neq_Model;
Mesh_Information.LM_Model=LM_Model;
Mesh_Information.nnp_max_Model=nnp_max_Model;
Mesh_Information.ngp_max_Model=ngp_max_Model;
Mesh_Information.nsurf_Model=nsurf_Model;
Mesh_Information.ngp_surf_max_Model=ngp_surf_max_Model;

Mesh_Information.gp_X=gp_X;
Mesh_Information.gp_Y=gp_Y;
Mesh_Information.gp_Z=gp_Z;

%%

% fprintf('Number of Nodes=%d\nNumer of Elements=%d\nNumer of Equations=%d\n',nnp_Model,nel_Model,neq_Model)
% fprintf('Total no of Gauss points=%d\n',sum(Type(3,:)))
% fprintf('%s\n',Line)

end


function[Mesh_Information]=Element_Converter_8_to_Any(Mesh_Information,Element_Information,Type_In)
%%

cd Element_Library
list=dir('*.m');
ElL=length(list);


for i=1:ElL
    fcnName=list(i).name;
    fcnName=erase(fcnName,'.m');
    
    ElNum=str2num(fcnName(end-1:end));
    if length(ElNum)==0
        ElNum=str2num(fcnName(end));
    end
    
    if ElNum==Type_In
        [INTP,~,W]=eval(fcnName);
        break
    end

end

cd ..
%%

Nodes=Mesh_Information.Nodes;
S_Type=Mesh_Information.S_Type;
IEN=Mesh_Information.IEN;
Type=Mesh_Information.Type;
BC_Nodes=Mesh_Information.BC_Nodes;


Check=0;

Type(1,:)=1;


[W,Shape]=Integration_Points_Data_3D_Brick_Quadratic_07_Local(INTP,W);

[Type,gp_X,gp_Y,gp_Z]=Gauss_Point_Update_locations_3D_L(Nodes,IEN,Type,Shape,W);

NN=gp_X*0;
c=1;
for e=1:size(gp_X,2)
    for gp=1:size(gp_X,1)

        if NN(gp,e)==0
            XYZ=[gp_X(gp,e) gp_Y(gp,e) gp_Z(gp,e)];
            Name=c;
            NN(gp,e)=c;
            
            Nodes_New(c,1:4)=[c XYZ];

            for e2=1:size(gp_X,2)
                for gp2=1:size(gp_X,1)
                    if NN(gp2,e2)==0
                        XYZ2=[gp_X(gp2,e2) gp_Y(gp2,e2) gp_Z(gp2,e2)];
                        if norm(XYZ-XYZ2)<10^-5
                           NN(gp2,e2)=c;
                        end
                    end
                end
            end


            c=c+1;
        end
    end
end


    Map=Nodes(:,1)*0;
    for i=1:length(Nodes)
        for j=1:length(Nodes_New)
             if norm(Nodes(i,2:4)-Nodes_New(j,2:4))<10^-10 
                Map(i)=j;
             end
        end
    end


%%
if Check==1
    Case=0;
    for i=1:length(Nodes_New)
        for j=1:length(Nodes_New)
             if norm(Nodes_New(i,2:4)-Nodes_New(j,2:4))<10^-10 && i~=j
             Case=Case+1;
             end
        end
    end
    if Case>0
        error
    end
end
IEN_New=NN;

Type_New=Type;
Type_New(1,:)=Type_In;
Type_New(2,:)=length(W);
Type_New(3,:)=length(W);

BC_Nodes_New=BC_Nodes;
for Bc=1:size(BC_Nodes,2)
    for n=1:size(BC_Nodes,1)
        if BC_Nodes(n,Bc)>0
            Previous_Node=BC_Nodes(n,Bc);
            New_Node=Map(Previous_Node);
            BC_Nodes_New(n,Bc)=New_Node;
        end
    end
end



S_Type=[S_Type [0;0;0]];

S_Nodes_In=zeros(Nodes_New(end,1),max(S_Type(1,:)));
SN=1;
s_count=1;
for ss=1:size(S_Type,2)-1
    
        e=S_Type(2,ss);
        eType=Type_New(1,e);
        snl=S_Type(3,ss);
        
        nsn=Element_Information.Surf_nnp_Elemental(snl,eType);
        
        n_set=Element_Information.Surf_Nodes_Elemental(1:nsn,snl,eType);
        
        Surf_nodes_e=IEN_New(n_set,e);

        for i=1:length(Surf_nodes_e)
            Nope=0;
            for j=1:size(S_Nodes_In,1)
                
                if S_Nodes_In(j,SN)==0
                    break
                end
                
                if Surf_nodes_e(i)==S_Nodes_In(j,SN)
                    Nope=1;
                    break
                end
            end
            if Nope==0
               S_Nodes_In(s_count,SN)=  Surf_nodes_e(i);
               s_count=s_count+1;
            end
            
        end
        
    if S_Type(1,ss+1)==S_Type(1,ss)        
        
    else
        s_count=1;
        SN=SN+1;
    end
    
end


for i=1:size(S_Nodes_In,1)
   
    if sum(S_Nodes_In(i,:))==0
        Limit=i-1;
        break
    end
end

S_Nodes_New=S_Nodes_In(1:Limit,:);


% S_Nodes_New=S_Nodes;
% 
% for Bc=1:size(S_Nodes,2)
%     for n=1:size(S_Nodes,1)
%         if BC_Nodes(n,Bc)>0
%             Previous_Node=BC_Nodes(n,Bc);
%             New_Node=Map(Previous_Node);
%             BC_Nodes_New(n,Bc)=New_Node;
%         end
%     end
% end


%% scatter3(Nodes_New(:,2),Nodes_New(:,3),Nodes_New(:,4))


    Mesh_Information.Nodes=Nodes_New;
    Mesh_Information.IEN=IEN_New;
    Mesh_Information.Type=Type_New;
    Mesh_Information.BC_Nodes=BC_Nodes_New;
    Mesh_Information.S_Nodes=S_Nodes_New;


%%
%                 for EN=1:10
%                 hold on
%                 IEN_New=NN;
%                 scatter3(Nodes_New(IEN_New(:,EN),2),Nodes_New(IEN_New(:,EN),3),Nodes_New(IEN_New(:,EN),4))
%                 daspect([1 1 1])
%                 
%                 rotate3d on
%                 pause
%                 %%
%                 for i=1:125
%                     
%                     scatter3(Nodes_New(IEN_New(i,EN),2),Nodes_New(IEN_New(i,EN),3),Nodes_New(IEN_New(i,EN),4),'filled')
%                     pause
%                 
%                 end
%                 
%                 
%                 end
end


function[W,shape]=Integration_Points_Data_3D_Brick_Quadratic_07_Local(INTP,W)

%% ----------- Region Integration parameters -------------------

        nnp=8;
        ngp=size(INTP,1);
        shape=zeros(ngp,nnp);      
           
    % Generating Shape and dshape matrix
        for i=1:ngp
             x=INTP(i,1);
             y=INTP(i,2);
             z=INTP(i,3);
             shape(i,1:nnp)=Shape_Functions(x,y,z);
        end
                                 
  
        
end

function[N]=Shape_Functions(x,y,z)
N=[
        (x*y)/8 - y/8 - z/8 - x/8 + (x*z)/8 + (y*z)/8 - (x*y*z)/8 + 1/8
        x/8 - y/8 - z/8 - (x*y)/8 - (x*z)/8 + (y*z)/8 + (x*y*z)/8 + 1/8
        x/8 + y/8 - z/8 + (x*y)/8 - (x*z)/8 - (y*z)/8 - (x*y*z)/8 + 1/8
        y/8 - x/8 - z/8 - (x*y)/8 + (x*z)/8 - (y*z)/8 + (x*y*z)/8 + 1/8
        z/8 - y/8 - x/8 + (x*y)/8 - (x*z)/8 - (y*z)/8 + (x*y*z)/8 + 1/8
        x/8 - y/8 + z/8 - (x*y)/8 + (x*z)/8 - (y*z)/8 - (x*y*z)/8 + 1/8
        x/8 + y/8 + z/8 + (x*y)/8 + (x*z)/8 + (y*z)/8 + (x*y*z)/8 + 1/8
        y/8 - x/8 + z/8 - (x*y)/8 - (x*z)/8 + (y*z)/8 - (x*y*z)/8 + 1/8
];
end

function[Type,gp_X,gp_Y,gp_Z]=Gauss_Point_Update_locations_3D_L(Nodes,IEN,Type,Shape,W)




%% Update Type by putting gauss point numbers
    SIEN=size(IEN);
    nel 	= SIEN(2);
    for e=1:nel

        eT=Type(1,e);

        gpn=length(W(:,1));
        for i=1:length(W(:,1))
            if W(i,eT)==0
                gpn=i-1;
                break
            end
        end

        Type(3,e)=gpn;

    end
  
    
%% Get gp locations    
% Memory Preallocation
    gp_X=zeros(max(Type(3,:)),nel);
    gp_Y=gp_X;
    gp_Z=gp_X;

% Each element calculation  

for e=1:nel

    nn=Type(2,e); % No of Nodes
    xyz=[Nodes(IEN(1:nn,e),2) Nodes(IEN(1:nn,e),3) Nodes(IEN(1:nn,e),4)];
    
    gp_X(1:Type(3,e),e)=Shape(1:Type(3,e),1:Type(2,e),Type(1,e))*xyz(:,1);
    gp_Y(1:Type(3,e),e)=Shape(1:Type(3,e),1:Type(2,e),Type(1,e))*xyz(:,2);
    gp_Z(1:Type(3,e),e)=Shape(1:Type(3,e),1:Type(2,e),Type(1,e))*xyz(:,3);
    
end  


end



function[Mesh_Information,N_map,V_map]=Node_Rearranger(Element_Information,Mesh_Information,R_Ns)


Nodes=Mesh_Information.Nodes;
IEN=Mesh_Information.IEN;
Type=Mesh_Information.Type;
S_Nodes=Mesh_Information.S_Nodes;
BC_Nodes=Mesh_Information.BC_Nodes;
Reg_Nodes=Mesh_Information.Reg_Nodes;

N=size(R_Ns);
S=size(Nodes);
nnp=S(1);

%% Set up Rearranged Nodes
R_Nodes=Nodes*0;
count=1;
for i=1:nnp
    
    if sum(intersect(Nodes(i,1),R_Ns))==0
        
        R_Nodes(count,1:4)=Nodes(i,1:4);
        
        count=count+1;
        
    end
    
end

Ex_Nn=count-1;  % Max Extempt Node number

for i=1:nnp
    if sum(intersect(Nodes(i,1),R_Ns))>0
        
        R_Nodes(count,1:4)=Nodes(i,1:4);
        
        count=count+1; 
    end
end

if nnp-count~=-1
    error('Something went Wrong!');
end

%% Get the Mapping Data     Nodes > R_Nodes
N_map=zeros(max(Nodes(:,1)),1);
for i=1:nnp
   
    N_map(R_Nodes(i,1),1)=i;

end

%% Set up Nodes and IEN
R_Nodes(:,1)=Nodes(:,1);

R_IEN=IEN*0;
S=size(IEN);
nel=S(2);

for e=1:nel

    for n=1:Element_Information.nnp_Elemental(Type(1,e))
    
        R_IEN(n,e)=N_map(IEN(n,e),1);
        
    end
    
end
    
ndof=3;

% Update LM too
    LM=zeros(max(Element_Information.nnp_Elemental(Type(1,:)))*ndof,nel);
    for e = 1:nel
        nen=Element_Information.nnp_Elemental(Type(1,e));
        for j = 1:nen
            for m = 1:ndof
                ind = (j-1)*ndof + m;
                 LM(ind,e) = ndof*R_IEN(j,e) - ndof + m;
            end
        end
    end  
    
% % Update IEN_Colm too    
%     IEN_Colm=zeros(sum(Type(2,:)),1);
%     count=1;
%     for e=1:nel
%         IEN_Colm(count:count+Type(2,e)-1,1)=IEN(1:Type(2,e),e);
%         count=count+Type(2,e);
%     end   


%% Set up BC and Surface nodes
R_BC_Nodes=BC_Nodes*0;
S=size(BC_Nodes);
for Bn=1:S(2)
    
    for n=1:S(1)
        
        if BC_Nodes(n,Bn)==0
            break
        end
        
        R_BC_Nodes(n,Bn)=N_map(BC_Nodes(n,Bn),1);
        
    end
    
end

R_S_Nodes=S_Nodes*0;
S=size(S_Nodes);
for Sn=1:S(2)

    for n=1:S(1)
        
        if S_Nodes(n,Sn)==0
            break
        end
        
        R_S_Nodes(n,Sn)=N_map(S_Nodes(n,Sn),1);
        
    end
    
end

%% Set Up Region Nodes
    R_Reg_Nodes=Reg_Nodes*0;
    S=size(Reg_Nodes);
    for j=1:S(2) 
        for i=1:S(1)
            if Reg_Nodes(i,j)~=0

                R_Reg_Nodes(i,j)=N_map(Reg_Nodes(i,j));

            end
        end
    end



%% Rearrange Vectors
V_map=zeros(length(N_map)*3,1);

for i=1:length(N_map)
    
   V_map(3*i-2,1)=N_map(i)*3-2;
   V_map(3*i-1,1)=N_map(i)*3-1;
   V_map(3*i,1)=N_map(i)*3;
    
end    


Mesh_Information.Nodes=R_Nodes;
Mesh_Information.IEN=R_IEN;
Mesh_Information.S_Nodes=R_S_Nodes;

Mesh_Information.BC_Nodes=R_BC_Nodes;
Mesh_Information.Reg_Nodes=R_Reg_Nodes;
Mesh_Information.LM_Model=LM;



end

            
function[Element_Information,Nodes_Elemental,nnp_Elemental,ngp_Elemental,Intp_Elemental,Shape,DShape,W,Surf_Nodes_Elemental,Surf_Intp_Elemental,Surf_ngp_Elemental,Surf_Shape,Surf_Dshape,Surf_W,Par_Coord]=Element_Loader
% 1 - Linear Tetra                  - 4 Node
% 2 - Quadratic Tetra               - 10 Node
% 3 - Linear Wedge                  - 6 Node    - 6  GP
% 4 - Quadratic Wedge               - 15 Node   - 21 GP
% 5 - Pyramid                       - 18 Node   - 21 GP
% 6 - Linear Brick                  - 8 Node    - 8  GP
% 7 - Quadratic Brick               - 20 Node   - 27 GP
% 8 - Enriched Quadratic Brick      - 27 Node   - 27 GP  
% 9 - Manual Element 9 Node Brick   - 9 Node    - 8  GP

%% Element Library

MaxGp=27;
MaxN=20;
MaxS=6;
MaxSGp=4;
Maxnnp=27;

cd Element_Library
list=dir('*.m');

ElL=length(list);

Nodes_Elemental=zeros(Maxnnp,3,ElL);
Intp_Elemental=zeros(MaxSGp,3,ElL);

W=zeros(MaxGp,ElL);
Shape=zeros(MaxGp,MaxN,ElL);
DShape=zeros(3,MaxN,ElL);
nnp_Elemental=zeros(1,ElL);
ngp_Elemental=zeros(1,ElL);

Surf_ngp_Elemental=zeros(6,ElL);
Surf_nnp_Elemental=zeros(6,ElL);
Par_Coord=zeros(3,2,8,ElL);
Surf_W=zeros(MaxSGp,ElL);
Surf_Shape=zeros(MaxSGp,MaxN,MaxS,ElL);
Surf_Dshape=zeros(3,MaxN,MaxSGp,MaxS,ElL);
% Surf_Dshape=zeros(3,1,1,1,ElL);

Surf_Nodes_Elemental=zeros(round(Maxnnp/3),MaxS,ElL);
Surf_Intp_Elemental=zeros(MaxSGp,3,6,ElL);



Line='------------------------------------------------------------------------------------------------------------------------------------------------------------';
fprintf('\n%s\n',Line)
fprintf('Loading elements:\n')
fprintf('%s\n',Line)
fprintf('%60s   %3s     %16s     %19s       %12s       %12s\n','Imported Function Name','DOF','Nodes in Element','Surfaces in Element','Integration points','Element Type') 
fprintf('%s\n',Line)

cd ..

for i=1:ElL
   
    fcnName=list(i).name;
    fcnName=erase(fcnName,'.m');
    
    %% For Standard FEM
    if fcnName(end)~='X'  
        
        [NodesO,INTPO,WO,ShapeO,DshapeO,Surf_Nodes_ElementalO,Surf_INTPO,Surf_WO,Surf_ShapeO,Surf_DshapeO,Par_CoordO]=eval(fcnName);

        ElNum=str2num(fcnName(end-1:end));
        if size(ElNum)==0
            ElNum=str2num(fcnName(end));
        end

        % Region properties
        SE=size(ShapeO);
        SD=size(DshapeO);
        ngp=SE(1);
        nnp=SE(2);
        ndof=SD(1);
        % ngpS=size(Surf_ShapeO,1);

        W(1:ngp,ElNum)=WO;
        Shape(1:ngp,1:nnp,ElNum)=ShapeO;
        DShape(1:ndof,1:nnp,1:ngp,ElNum)=DshapeO;
        nnp_Elemental(1,ElNum)=nnp;
        ngp_Elemental(1,ElNum)=ngp;
        Nodes_Elemental(1:nnp,:,ElNum)=NodesO;

        Intp_Elemental(1:ngp,:,ElNum)=INTPO;
        
        % Surface properties
%         L=size(Surf_DshapeO);
%         if length(L)<3
%             sn=1;
%             Surf_ngp_ElementalO=1;
%         elseif length(L)<4
%             sn=1;
%         else
%             sn=L(4);
%             Surf_ngp_ElementalO=L(3);
%         end

        sn=size(Surf_DshapeO,4); % Number of surfaces
        S=size(Surf_WO);
        for ss=1:S(2)
            gpsCount=0;
            for gps=1:S(1)
                if Surf_WO(gps,ss)>0
                   gpsCount=gpsCount+1;
                   Surf_W(gpsCount,ss,ElNum)=Surf_WO(gpsCount,ss);
                   Surf_Intp_Elemental(gpsCount,:,ss,ElNum)=Surf_INTPO(gpsCount,:,ss);
                end
                Surf_ngp_Elemental(ss,ElNum)=gpsCount;
            end
        end
        maxSgp=max(Surf_ngp_Elemental(:,ElNum));

        %Surf_Nodes_ElementalO,Surf_INTPO

        nnpS=size(Surf_Nodes_ElementalO,1);

        for ss=1:S(2)
            count=0;
            for nn=1:nnpS
                if Surf_Nodes_ElementalO(nn,ss)>0
                    Surf_nnp_Elemental(ss,ElNum)=count+1;
                    count=count+1;
                end
                
            end
        end
        Surf_Nodes_Elemental(1:nnpS,1:sn,ElNum)=Surf_Nodes_ElementalO;




        Surf_Dshape(1:ndof,1:nnp,1:maxSgp,1:sn,ElNum)=Surf_DshapeO;
        Surf_Shape(1:maxSgp,1:nnp,1:sn,ElNum)=Surf_ShapeO;
        
        Par_Coord(:,:,1:sn,ElNum)=Par_CoordO;

        fprintf('%60s   %3d     %11d     %17d       %15d         %15d\n',fcnName,ndof,nnp,sn,ngp,ElNum)

        % Clearing var
        clear WO shapeO DshapeO Surf_WO Surf_ngp_ElementalO Surf_ShapeO Par_CoordO Surf_DshapeO
    
    end
    
    
end


% Element_Information=struct('Nodes_Elemental','nnp_Elemental','ngp_Elemental','Intp_Elemental','Shape,DShape','W,Surf_Nodes_Elemental','Surf_Intp_Elemental','Surf_ngp_Elemental','Surf_Shape','Surf_Dshape','Surf_W','Par_Coord',[]);

Element_Information.Nodes_Elemental=Nodes_Elemental;
Element_Information.nnp_Elemental=nnp_Elemental;
Element_Information.ngp_Elemental=ngp_Elemental;
Element_Information.Intp_Elemental=Intp_Elemental;
Element_Information.Shape_Elemental=Shape;
Element_Information.DShape_Elemental=DShape;
Element_Information.W_Elemental=W;
Element_Information.Surf_Nodes_Elemental=Surf_Nodes_Elemental;
Element_Information.Surf_Intp_Elemental=Surf_Intp_Elemental;
Element_Information.Surf_ngp_Elemental=Surf_ngp_Elemental;
Element_Information.Surf_nnp_Elemental=Surf_nnp_Elemental;
Element_Information.Surf_Shape_Elemental=Surf_Shape;
Element_Information.Surf_Dshape_Elemental=Surf_Dshape;
Element_Information.Surf_W_Elemental=Surf_W;
Element_Information.Par_Coord_Elemental=Par_Coord;

fprintf('%s\n\n',Line)

end

function[Mesh_Information,Nodes,IEN,Type,S_Nodes,S_Type,Reg_Nodes,Reg_Elems,BC_Nodes]=Abaqus_File_Processor_4(Name,S_Names,Reg_Name,BC_Name,Element_Information)
%%)

% Name='FEM_Matlab_25.inp';
% S_Names={'Sx','S1'};
% Reg_Name={'Sx','S1'};
% BC_Name={'Sx','S1'};

fprintf('\n\n------------------------------------------------------\n')
fprintf('Converting Abaqus to Matlab\n')
fprintf('------------------------------------------------------\n')
fprintf('Developer: Adnan Shahriar\n           Contact: adnan.shahriar@utsa.edu\n           Website: adnansh.com\n           Department of Mechanical Engineering, UTSA, Texas\n')
fprintf('------------------------------------------------------\n')

% fprintf(2,'Reading and Importing Abaqus Input File...\n')
fprintf('Reading and Importing Abaqus Input File...\n')

% Element assignment
E_Type={'C3D4','C3D4R','C3D10','C3D10R','C3D6','C3D6R','C3D15','C3D15R','C3D8','C3D8R','C3D20','C3D20R' ,'DC3D20', 'DC3D8'
           2      2       3        3       4      4       5        5      6       6       7        7        7         6};



% E_Type=[1  2 3  4  5 5  6  8;   % Element Type
%         3 10 6 15 18 8 20 27;   % Node Number
%         1  1 6 21 21 8 27 20];  % Gauss Points
% 
% % 8 and 20 Node Brick Element Surface Data    
% Surf6=[2     3     6     7  % g=1
%        1     4     5     8  % g=-1
%        3     4     7     8  % h=1
%        1     2     5     6  % h=-1
%        5     6     7     8  % r=1
%        1     2     3     4];% r=-1    
%    
% % 6 and 15 Node Wedge Element Data   
% Surf3=[9 9 9;9 9 9; 9 9 9; % x, y surface are not incorporated
%        4 5 6;       % r=1
%        1 2 3];      % r=-1
    

   
% Import    
    F_Name = Importfile(Name);

% Node Import
    fprintf('Importing Nodes...')
    
    for i=1:length(F_Name)
        if F_Name(i,1)=='*Node'
            N_s=i+1;
        elseif F_Name(i,1)=='*Element'
            N_e=i-1;
            break
        end
    end

    Nodes=zeros(N_e-N_s,4);
    count=1;
    for i=N_s:N_e
        Nodes(count,1)=str2num(F_Name(i,1));
        Nodes(count,2)=str2num(F_Name(i,2));
        Nodes(count,3)=str2num(F_Name(i,3));
        Nodes(count,4)=str2num(F_Name(i,4));
        count=count+1;
    end
    fprintf('Done\n')

%% Element Import
fprintf('Importing Elements...')

    E_d=0;          % E_d store all Element array data
                    % [ Element1 start      Element2 start ...
                    %   Element1 end        Element2 end   ...
                    %   Element1 Type       Element2 Type  ...]
    count=1;
    for i=N_e:length(F_Name)
        
        % Get E_d
        if F_Name(i,1)=='*Element'
            E_d(1,count)=i+1;
            if count>1
                E_d(2,count-1)=i-1; 
            end

            for eT=1:size(E_Type,2)
                ET=sprintf('type=%s',E_Type{1,eT});
                if F_Name(i,2)==ET
                    E_d(3,count)=E_Type{2,eT};
                end

            end

%             if F_Name(i,2)=='type=C3D8' || F_Name(i,2)=='type=C3D8R'
%                 E_d(3,count)=6;
%             elseif F_Name(i,2)=='type=C3D6' || F_Name(i,2)=='type=C3D6R'
%                 E_d(3,count)=3;
%             elseif F_Name(i,2)=='type=C3D20' || F_Name(i,2)=='type=C3D20R'
%                 E_d(3,count)=7;
%             elseif F_Name(i,2)=='type=C3D15' || F_Name(i,2)=='type=C3D15R'
%                 E_d(3,count)=4;
%             end

            count=count+1;
        elseif F_Name(i,1)=='*Nset' || F_Name(i,1)=='*End Part'
            E_d(2,count-1)=i-1;
            break
        end
    end

    if length(E_d)~=3
        fprintf(1,'Most likely element type from Abaqus is unknown')
        error
    end
    
%%

    % Generate 
    IEN=0;
    S=size(E_d);
    
    for eType=1:S(2)    % Go through each element types
       
        i=E_d(1,eType);
        nne=Element_Information.nnp_Elemental(E_d(3,eType));
        
        while i<=E_d(2,eType)
            
            enum=str2num(F_Name(i,1));
            
            n=1;
            j=2;
            while n<=nne
                IEN(n,enum)=str2num(F_Name(i,j));
                j=j+1;
                n=n+1;
                if j>16 && nne>16 %E_d(3,eType)>6 
                    j=1;
                    i=i+1;
                end
            end
            
            Type(1,enum)=E_d(3,eType);
            Type(2,enum)=Element_Information.nnp_Elemental(E_d(3,eType));
            Type(3,enum)=Element_Information.ngp_Elemental(E_d(3,eType));
            Type(4,enum)=1;
            Type(5,enum)=0;
            
            i=i+1;
        end
        
        
    end
    fprintf('Done\n')


    
%% Surface Import   
    % Get Surface Nodes
    %fprintf('\n--------------------------------------\n')
    fprintf('Importing surface sets...\n')
    %fprintf('--------------------------------------\n')

    S_Nodes=0;
    for SN=1:length(S_Names)
       
        % Check Start and End
        Sd=sprintf('nset=%s',S_Names{SN});
        St=0;
        En=0;
        for i=1:length(F_Name)
            if F_Name(i,2)==Sd
                St=i+1;
            elseif F_Name(i,1)=='*Elset' && St>0
                En=i-1;
                break
            end
        end
        
        if St==0 && En==0
            fprintf(2,'   Warning! Nodes of "%s" surface set can not be detected\n',S_Names{SN});
        else

                % Assign Nodes
                if F_Name(St-1,4)=='generate'
                    st=str2num(F_Name(St,1));
                    en=str2num(F_Name(St,2));
                    ns=[st:str2num(F_Name(St,3)):en];
                    S_Nodes(1:length(ns),SN)=ns;     
                else
                    count=1;
                    for i=St:En
                        for j=1:16
                            if F_Name(i,j)~=""
                                S_Nodes(count,SN)=str2num(F_Name(i,j));
                                count=count+1;
                            end
                        end
                    end
                end  
                
                fprintf('   Success! Nodes of "%s" surface set has been detected\n',S_Names{SN});
        end
    end
        

        % Get Surface Elements
        count=1;
        S_Type_B=0;
        for SN=1:length(S_Names)

            Sd=sprintf('elset=%s',S_Names{SN});

            % Check Start and End
            St=0;
            En=0;
            for i=1:length(F_Name)
                if F_Name(i,2)==Sd
                    St=i+1;
                elseif (F_Name(i,1)=='*Elset' || F_Name(i,1)=='*Nset' || F_Name(i,1)=='*End Assembly' || F_Name(i,1)=='*Surface') && St>0
                    En=i-1;
                    break
                end
            end

            % Breaking loop if no surface is detected
            if St==0 && En==0
                fprintf(2,'   Warning! Elements of "%s" surface set can not be detected\n',S_Names{SN});
            else

                % Assign Elements
                if F_Name(St-1,4)=='generate'
                    st=str2num(F_Name(St,1));
                    en=str2num(F_Name(St,2));
                    ns=[st:str2num(F_Name(St,3)):en];
                    S_Type_B(count:count+length(ns)-1,2)=ns; 
                    S_Type_B(count:count+length(ns)-1,1)=ns*0+SN;
                    count=count+length(ns);
                else

                    for i=St:En
                        for j=1:16
                            if F_Name(i,j)~=""
                                S_Type_B(count,2)=str2num(F_Name(i,j));
                                S_Type_B(count,1)=SN;
                                count=count+1;
                            end
                        end
                    end

                end 
                fprintf('   Success! Elements of "%s" surface set has been detected\n',S_Names{SN});
            end
        end
    
        
    %% Hardest Part! Surface type extraction!
    S=size(S_Type_B);
    
    c=1;
    for i=1:S(1)               % For each Element

        SN=S_Type_B(i,1);         % Extract Surface Number
        
            % Breaking loop if no surface is detected
            if S_Type_B==0
                break
            end

            e=S_Type_B(i,2);      % Extract the Element

            % Later this Stupid algorithm should be replaced!
            count=1;
            L=0;
            nn=Element_Information.nnp_Elemental(Type(1,e));
            for k=1:nn   % Extract the array no by checking each Nodes from List

                for l=1:length(S_Nodes)         % Run through all nodes and check!
                    if S_Nodes(l,SN)==IEN(k,e)
                        L(count,1)=k;
                        count=count+1;
                        break
                    end
                end

            end

            ss=size(Element_Information.Surf_ngp_Elemental,1);
                cs=1;
                for l=1:ss

                    SN_Limit=Element_Information.Surf_nnp_Elemental(l,Type(1,e));

                    if SN_Limit>0
                        SNS=Element_Information.Surf_Nodes_Elemental(1:SN_Limit,l,Type(1,e));
                        if length(intersect(SNS,L))==SN_Limit
                           S_Type_BU(c,3)=l;
                           S_Type_BU(c,1:2)=S_Type_B(i,1:2);
      
                           cs=cs+1;
                           c=c+1;
                        end
                    end
                    
                end

                if cs>3
                    fprintf(2,'Something wrong with Surface type extraction\n')
                end


    end

    if S_Type_B~=0
        S_Type_B=S_Type_BU;
    end

    S_Type=S_Type_B';


%% Boundary Extract
    BC_Nodes=0;
    %fprintf('\n--------------------------------------\n')
    fprintf('Importing Boundary sets...\n')
    %fprintf('--------------------------------------\n')
    for BN=1:length(BC_Name)
        
        Sd=sprintf('nset=%s',BC_Name{BN});
        St=0;
        En=0;
        for i=1:length(F_Name)
            if F_Name(i,2)==Sd
                St=i+1;
            elseif (F_Name(i,1)=='*Elset' || F_Name(i,1)=='*Nset') && St>0
                En=i-1;
                break
            end
        end
        
        if St==0 && En==0
            fprintf(2,'   Warning! "%s" boundary set can not be detected\n',BC_Name{BN});
        else
            
        
        
        if F_Name(St-1,4)=='generate' 
            st=str2num(F_Name(St,1));
            en=str2num(F_Name(St,2));
            ns=[st:str2num(F_Name(St,3)):en];
            BC_Nodes(1:length(ns),BN)=ns;     
        else        
            count=1;
            for i=St:En
                for j=1:16
                    if F_Name(i,j)~=""
                        BC_Nodes(count,BN)=str2num(F_Name(i,j));
                        count=count+1;
                    end
                end
            end
        end  
        
        
        fprintf('   Success! "%s" boundary set has been detected\n',BC_Name{BN});
        
        end
        
    end
    
 Type(end-1,:)=1;   
 Type(end,:)=0;   
    
%% Material Extract

% for MN=1:length(Mat_Names)
%     
%     Sd=sprintf('name=%s',Mat_Names{MN});
%     
%     % Get element data
%     for i=1:length(F_Name)
% 
%         if F_Name(i,2)==Sd
%             rho(MN,1)=str2num(F_Name(i+2,1));
%             MoE(MN,1)=str2num(F_Name(i+4,1));
%             PoRat(MN,1)=str2num(F_Name(i+4,2));
%             break
%         end
%         
%     end
%     
%     Sd2=sprintf('material=%s',Mat_Names{MN});
%     % Update the Type 
%     for i=1:length(F_Name)
% 
%         % Look for the Section part
%         if F_Name(i,3)==Sd2
%             
%             Elem_search=F_Name(i,2);
%             
%             % Look with the element name from Section part
%             for j=1:length(F_Name)
%                 
%                 % Got the part! Now Update Type
%                 if F_Name(j,1)=='*Elset' && F_Name(j,2)==Elem_search
%                     
%                     if F_Name(j,4)=='generate' || F_Name(j,3)=='generate'
%                        
%                         EList=[str2num(F_Name(j+1,1)):str2num(F_Name(j+1,3)):str2num(F_Name(j+1,2))];
%                         
%                         Type(4,EList)=MN;
%                         
%                     else
%                         
%                         R_Up=1;
%                         while R_Up~=0
%                             for k=1:16
%                                 if F_Name(j+1,k)==""
%                                     break
%                                 end
%                                 EList=str2num(F_Name(j+1,k));
%                                 Type(4,EList)=MN;
%                                 
%                             end
% 
%                             j=j+1;
%                             
%                             if F_Name(j+1,1)=="*Nset"
%                                 R_Up=0;
%                             end
%                             if F_Name(j+1,1)=="*Elset"
%                                 R_Up=0;
%                             end
%                             
%                         end
%                     end
%                     
%                     break
%                 end
%   
%             end
%                break
%         end
%         
%     end
%           
% end


%% Region Extract

%% Region Element Extract
    %fprintf('\n--------------------------------------\n')
    fprintf('Importing Region sets...\n')
    %fprintf('--------------------------------------\n')
    
    Reg_Elems=0;
    for RN=1:length(Reg_Name)
        
        Sd=sprintf('elset=%s',Reg_Name{RN});
        St=0;
        En=0;
        for i=1:length(F_Name)
            if F_Name(i,2)==Sd
                St=i+1;
            elseif sum(str2num(F_Name(i,1)))==0 && St>0
                En=i-1;
                break
            end
        end
        
        
       if St==0 && En==0  
       else
        
            if F_Name(St-1,4)=='generate' 
                st=str2num(F_Name(St,1));
                en=str2num(F_Name(St,2));
                ns=[st:str2num(F_Name(St,3)):en];
                Reg_Elems(1:length(ns),RN)=ns;     
            else        
                count=1;
                for i=St:En
                    for j=1:16
                        if F_Name(i,j)~=""
                            Reg_Elems(count,RN)=str2num(F_Name(i,j));
                            count=count+1;
                        end
                    end
                end
            end  
       end
        
    end
    

        

%% Region Node Extract
    c=0;
    Siz=size(Reg_Elems);
    for RN=1:length(Reg_Name)
        
        Sd=sprintf('nset=%s',Reg_Name{RN});
        Sd2=sprintf('elset=%s',Reg_Name{RN});
        St=0;
        En=0;
        for i=1:length(F_Name)
            if F_Name(i,2)==Sd %|| F_Name(i,2)==Sd2
                St=i+1;
            elseif (F_Name(i,1)=='*Elset' || F_Name(i,1)=='*Nset') && St>0
                En=i-1;
                break
            end
        end
        
        
        if St==0 && En==0
            fprintf(2,'   Warning! "%s" Region set can not be detected\n',Reg_Name{RN});
            c=c+1;
        else
 
            if St>0

                if F_Name(St-1,4)=='generate' 
                    st=str2num(F_Name(St,1));
                    en=str2num(F_Name(St,2));
                    ns=[st:str2num(F_Name(St,3)):en];
                    Reg_Nodes(1:length(ns),RN)=ns;     
                else        
                    count=1;
                    for i=St:En
                        for j=1:16
                            if F_Name(i,j)~=""
                                Reg_Nodes(count,RN)=str2num(F_Name(i,j));
                                count=count+1;
                            end
                        end
                    end
                end   

            elseif St==0

                Reg_Nodes=0;
                fprintf(2,'\n     WARNING! No Region Nodes were found! Probably, region has been created using Elements. Generate Externally if you want so much')
                fprintf('\n')            
            end
            
            
            fprintf('   Success! "%s" Region set has been detected\n',Reg_Name{RN});
        end
        
    end
    
    if c==length(Reg_Name)
       Reg_Nodes=0;
    end

    Mesh_Information.Nodes=Nodes;
    Mesh_Information.IEN=IEN;
    Mesh_Information.Type=Type;
    Mesh_Information.S_Nodes=S_Nodes;
    Mesh_Information.S_Type=S_Type;
    Mesh_Information.Reg_Nodes=Reg_Nodes;
    Mesh_Information.Reg_Elements=Reg_Elems;
    Mesh_Information.BC_Nodes=BC_Nodes;



    % fprintf('\n..Done\n')
    fprintf('------------------------------------------------------\n')

end    
        


%% Importing Function
function Ab2Mat = Importfile(filename)

startRow=1;
fid=fopen(filename);
g = textscan(fid,'%s','delimiter','\n');
fclose(fid);
endRow=length(g{1});


%IMPORTFILE Import numeric data from a text file as a matrix.
%   AB2MAT = IMPORTFILE(FILENAME) Reads data from text file FILENAME for
%   the default selection.
%
%   AB2MAT = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   Ab2Mat = importfile('Ab2Mat.txt', 1, 4493);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2020/05/30 17:35:22

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: text (%s)
%	column4: text (%s)
%   column5: text (%s)
%	column6: text (%s)
%   column7: text (%s)
%	column8: text (%s)
%   column9: text (%s)
%	column10: text (%s)
%   column11: text (%s)
%	column12: text (%s)
%   column13: text (%s)
%	column14: text (%s)
%   column15: text (%s)
%	column16: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
Ab2Mat = [dataArray{1:end-1}];

end


function[gp_X,gp_Y,gp_Z]=Gauss_Point_locations_3D(Element_Information,Mesh_Information)

    nnp_Elemental=Element_Information.nnp_Elemental;
    ngp_Elemental=Element_Information.ngp_Elemental;
    Type=Mesh_Information.Type;
    IEN=Mesh_Information.IEN;
    Nodes=Mesh_Information.Nodes;
    Shape=Element_Information.Shape_Elemental;
    
    SIEN=size(IEN);
    nel 	= SIEN(2);

%% Memory Preallocation
    gp_X=zeros(ngp_Elemental(Type(1,1)),nel);
    gp_Y=gp_X;
    gp_Z=gp_X;

%% Each element calculation  

for e=1:nel

    nn=nnp_Elemental(Type(1,e)); % No of Nodes
    xyz=[Nodes(IEN(1:nn,e),2) Nodes(IEN(1:nn,e),3) Nodes(IEN(1:nn,e),4)];
    
    gp_X(1:ngp_Elemental(Type(1,e)),e)=Shape(1:ngp_Elemental(Type(1,e)),1:nnp_Elemental(Type(1,e)),Type(1,e))*xyz(:,1);
    gp_Y(1:ngp_Elemental(Type(1,e)),e)=Shape(1:ngp_Elemental(Type(1,e)),1:nnp_Elemental(Type(1,e)),Type(1,e))*xyz(:,2);
    gp_Z(1:ngp_Elemental(Type(1,e)),e)=Shape(1:ngp_Elemental(Type(1,e)),1:nnp_Elemental(Type(1,e)),Type(1,e))*xyz(:,3);
    
        
end  

end
