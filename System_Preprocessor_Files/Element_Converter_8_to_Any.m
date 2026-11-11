function[Mesh_Information]=Element_Converter_8_to_Any(Mesh_Information,Element_Information,Type_In)
%%

Lib_Dir=fullfile(fileparts(fileparts(mfilename('fullpath'))),'Element_Library');  % <main folder>/Element_Library
list=dir(fullfile(Lib_Dir,'*.m'));
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
