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
