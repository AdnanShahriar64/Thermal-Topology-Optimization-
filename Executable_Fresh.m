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

      Material_Information.ground_albedo = 0.4;
      Material_Information.epsilon_gr = 0.85;


      Material_Information.StefBC=5.67*10^-8;


%% Model Input File
% Read Abaqus Input file

      % Name='SPL_Section_Linear_2.inp';
      % Name='SPL_Thermal-Copy-Copy2.inp';
      % Name='SPL_Thermal-Copy-Copy2-Copy.inp';
      Name='SPL_Linear_3.inp';
      S_Names={'S4','S2','S3'};
      % S_Names={'S4'};
      Reg_Name={'Sx','S1'};
      BC_Name={'S3'};

      [Mesh_Information]=Abaqus_File_Processor_4(Name,S_Names,Reg_Name,BC_Name,Element_Information);



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
