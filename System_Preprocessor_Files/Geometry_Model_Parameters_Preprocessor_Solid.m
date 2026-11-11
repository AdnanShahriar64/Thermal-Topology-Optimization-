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
