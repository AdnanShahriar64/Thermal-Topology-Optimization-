function[Nodes,INTP,W,Shape,Dshape,Surf_Nodes,Surf_INTP,Surf_W,Surf_Shape,Surf_Dshape,Parametric_Coordinates]=Integration_Points_Data_3D_Hexahedron_Node_8_Type_6

Nodes=[
-1.00000  -1.00000  -1.00000
 1.00000  -1.00000  -1.00000
 1.00000   1.00000  -1.00000
-1.00000   1.00000  -1.00000
-1.00000  -1.00000   1.00000
 1.00000  -1.00000   1.00000
 1.00000   1.00000   1.00000
-1.00000   1.00000   1.00000
];
INTP=[
    -1.00000  -1.00000  -1.00000
     1.00000  -1.00000  -1.00000
     1.00000   1.00000  -1.00000
    -1.00000   1.00000  -1.00000
    -1.00000  -1.00000   1.00000
     1.00000  -1.00000   1.00000
     1.00000   1.00000   1.00000
    -1.00000   1.00000   1.00000
];
W=[
     1.0000000000
     1.0000000000
     1.0000000000
     1.0000000000
     1.0000000000
     1.0000000000
     1.0000000000
     1.0000000000
];

ngp=8;
nnp=8;
Shape=zeros(ngp,nnp);
Dshape=zeros(3,nnp,ngp);

for i=1:ngp
     x=INTP(i,1);
     y=INTP(i,2);
     z=INTP(i,3);
     Shape(i,1:nnp)=Shape_Functions(x,y,z);

     [dNbydx,dNbydy,dNbydz]=Shape_Functions_Derivetives(x,y,z);
     Dshape(1,1:nnp,i)=dNbydx;
     Dshape(2,1:nnp,i)=dNbydy;
     Dshape(3,1:nnp,i)=dNbydz;
end

Surf_Nodes=zeros(1,6);
    Surf_Nodes_L=[
     2
     3
     6
     7
    ];
Surf_Nodes(1:4,1)=Surf_Nodes_L;
    Surf_Nodes_L=[
     1
     4
     5
     8
    ];
Surf_Nodes(1:4,2)=Surf_Nodes_L;
    Surf_Nodes_L=[
     3
     4
     7
     8
    ];
Surf_Nodes(1:4,3)=Surf_Nodes_L;
    Surf_Nodes_L=[
     1
     2
     5
     6
    ];
Surf_Nodes(1:4,4)=Surf_Nodes_L;
    Surf_Nodes_L=[
     5
     6
     7
     8
    ];
Surf_Nodes(1:4,5)=Surf_Nodes_L;
    Surf_Nodes_L=[
     1
     2
     3
     4
    ];
Surf_Nodes(1:4,6)=Surf_Nodes_L;
Surf_W=[
      1.000000      1.000000      1.000000      1.000000      1.000000      1.000000
      1.000000      1.000000      1.000000      1.000000      1.000000      1.000000
      1.000000      1.000000      1.000000      1.000000      1.000000      1.000000
      1.000000      1.000000      1.000000      1.000000      1.000000      1.000000
    ];
Surf_ngp=[
     4     4     4     4     4     4    ];
Surf_INTP=zeros(max(Surf_ngp),3,length(Surf_ngp));
INTP_Surf_L=[
      1.000000    -1.000000      -1.000000
      1.000000    1.000000      -1.000000
      1.000000    -1.000000      1.000000
      1.000000    1.000000      1.000000
    ];
Surf_INTP(1:4,1:3,1)=INTP_Surf_L;
INTP_Surf_L=[
      -1.000000    -1.000000      -1.000000
      -1.000000    1.000000      -1.000000
      -1.000000    -1.000000      1.000000
      -1.000000    1.000000      1.000000
    ];
Surf_INTP(1:4,1:3,2)=INTP_Surf_L;
INTP_Surf_L=[
      1.000000    1.000000      -1.000000
      -1.000000    1.000000      -1.000000
      1.000000    1.000000      1.000000
      -1.000000    1.000000      1.000000
    ];
Surf_INTP(1:4,1:3,3)=INTP_Surf_L;
INTP_Surf_L=[
      -1.000000    -1.000000      -1.000000
      1.000000    -1.000000      -1.000000
      -1.000000    -1.000000      1.000000
      1.000000    -1.000000      1.000000
    ];
Surf_INTP(1:4,1:3,4)=INTP_Surf_L;
INTP_Surf_L=[
      -1.000000    -1.000000      1.000000
      1.000000    -1.000000      1.000000
      1.000000    1.000000      1.000000
      -1.000000    1.000000      1.000000
    ];
Surf_INTP(1:4,1:3,5)=INTP_Surf_L;
INTP_Surf_L=[
      -1.000000    -1.000000      -1.000000
      1.000000    -1.000000      -1.000000
      1.000000    1.000000      -1.000000
      -1.000000    1.000000      -1.000000
    ];
Surf_INTP(1:4,1:3,6)=INTP_Surf_L;
Surf_Shape=zeros(max(Surf_ngp),8,length(Surf_ngp));
Surf_Dshape=zeros(3,8,max(Surf_ngp),length(Surf_ngp));

for ss=1:length(Surf_ngp)
    for gp=1:Surf_ngp(ss)
         x=Surf_INTP(gp,1,ss);
         y=Surf_INTP(gp,2,ss);
         z=Surf_INTP(gp,3,ss);
         Surf_Shape(gp,1:nnp,ss)=Shape_Functions(x,y,z);

     [dNbydx,dNbydy,dNbydz]=Shape_Functions_Derivetives(x,y,z);
     Surf_Dshape(1,1:nnp,gp,ss)=dNbydx;
     Surf_Dshape(2,1:nnp,gp,ss)=dNbydy;
     Surf_Dshape(3,1:nnp,gp,ss)=dNbydz;
    end
end

Parametric_Coordinates=zeros(3,2,length(Surf_ngp));
P_C_L=[
    0.000000    0.000000
    1.000000    0.000000
    0.000000    1.000000
];
Parametric_Coordinates(:,:,1)=P_C_L;
P_C_L=[
    0.000000    0.000000
    0.000000    1.000000
    1.000000    0.000000
];
Parametric_Coordinates(:,:,2)=P_C_L;
P_C_L=[
    0.000000    1.000000
    0.000000    0.000000
    1.000000    0.000000
];
Parametric_Coordinates(:,:,3)=P_C_L;
P_C_L=[
    1.000000    0.000000
    0.000000    0.000000
    0.000000    1.000000
];
Parametric_Coordinates(:,:,4)=P_C_L;
P_C_L=[
    1.000000    0.000000
    0.000000    1.000000
    0.000000    0.000000
];
Parametric_Coordinates(:,:,5)=P_C_L;
P_C_L=[
    0.000000    1.000000
    1.000000    0.000000
    0.000000    0.000000
];
Parametric_Coordinates(:,:,6)=P_C_L;
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
function[dNbydx,dNbydy,dNbydz]=Shape_Functions_Derivetives(x,y,z)
dNbydx=[
        y/8 + z/8 - (y*z)/8 - 1/8
        (y*z)/8 - z/8 - y/8 + 1/8
        y/8 - z/8 - (y*z)/8 + 1/8
        z/8 - y/8 + (y*z)/8 - 1/8
        y/8 - z/8 + (y*z)/8 - 1/8
        z/8 - y/8 - (y*z)/8 + 1/8
        y/8 + z/8 + (y*z)/8 + 1/8
        - y/8 - z/8 - (y*z)/8 - 1/8
];
dNbydy=[
        x/8 + z/8 - (x*z)/8 - 1/8
        z/8 - x/8 + (x*z)/8 - 1/8
        x/8 - z/8 - (x*z)/8 + 1/8
        (x*z)/8 - z/8 - x/8 + 1/8
        x/8 - z/8 + (x*z)/8 - 1/8
        - x/8 - z/8 - (x*z)/8 - 1/8
        x/8 + z/8 + (x*z)/8 + 1/8
        z/8 - x/8 - (x*z)/8 + 1/8
];
dNbydz=[
        x/8 + y/8 - (x*y)/8 - 1/8
        y/8 - x/8 + (x*y)/8 - 1/8
        - x/8 - y/8 - (x*y)/8 - 1/8
        x/8 - y/8 + (x*y)/8 - 1/8
        (x*y)/8 - y/8 - x/8 + 1/8
        x/8 - y/8 - (x*y)/8 + 1/8
        x/8 + y/8 + (x*y)/8 + 1/8
        y/8 - x/8 - (x*y)/8 + 1/8
];
end