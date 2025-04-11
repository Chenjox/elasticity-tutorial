

width = 48;
heightleft = 44;
heightright = 16;
DefineConstant[num_nodes=3]; // Definiere die Konstante, wenn sie nicht schon existiert.
//num_nodes= If (!Exists(num_nodes)) 
// 3
//EndIf

// Zuerst die Punkte
Point(1) = {0,0,0};
Point(2) = {width,heightleft,0};
Point(3) = {width,heightleft+heightright,0};
Point(4) = {0,heightleft,0};


// Dann die äußeren Linien
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(100) = {1,2,3,4};

Surface(100) = {100};

/*
Physical Curve("Left",4)   = {4};
Physical Curve("Right",2)  = {2};
Physical Curve("Bottom",1) = {1};
Physical Curve("Top",3)    = {3};

Physical Surface("Domain",1) = {100};
*/
// Da dune-gmsh4 keine Strings einlesen kann, wird folgendes gemacht:
Physical Curve(4)   = {4};
Physical Curve(2)  = {2};
Physical Curve(1) = {1};
Physical Curve(3)    = {3};

Physical Surface(1) = {100};

// // Viele Linien, und vorallem Knoten
Transfinite Line {1,3}  = num_nodes Using Progression 1;
Transfinite Line {2,4}  = num_nodes Using Progression 1;

Transfinite Surface {100};

Coherence;
// SetOrder 1;
Recombine Surface {100};
Mesh 2;
RenumberMeshNodes;
str_num_nodes = Sprintf("%03g",num_nodes);
// Save StrCat("filename.msh");
Save StrCat("cooksMembrane-",str_num_nodes,".msh");