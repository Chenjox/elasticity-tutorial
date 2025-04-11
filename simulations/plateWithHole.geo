
mm = 1e-3;
squareLength = 10 * mm;
circleRadius = 1 * mm;

DefineConstant[size_ratio=3];
meshSize = circleRadius / size_ratio;

// Square Coordinates
// Counter clockwise starting bottom left
Point(1) = { 0,0,0,meshSize };
Point(2) = { squareLength,0,0,meshSize };
Point(3) = { squareLength,squareLength,0,meshSize };
Point(4) = { 0,squareLength,0,meshSize };

// Additional Circle Coordinates
Point(5) = { squareLength - circleRadius, 0,0,meshSize};
Point(6) = { squareLength, circleRadius, 0,meshSize };

// Defining outer boundary segments
Line(1) = {1, 5};
Circle(2) = {5, 2, 6};
Line(3) = {6, 3};
Line(4) = {3, 4};
Line(5) = {4, 1};

// Defining the Boundary Curve
Curve Loop(6) = {1,2,3,4,5};
Plane Surface(7) = {6};

// This is the measure set, that needs to be meshed
Physical Surface(10) = {7};

Mesh.Algorithm = 6; // Frontal Delaunay

Mesh 2;
//RecombineMesh;

//Mesh.MshFileVersion = 2.0;
Save StrCat("plateWithHole.msh");

