// Gmsh project created on Sun Mar  3 15:46:55 2019
SetFactory("OpenCASCADE");

h = 2.0; // Element scale factor >> 1 is coarse
x = 100.0; 
y = 100.0;
saltthickness = 10.0; // Thickness of the salt unit
sandstonethickness = 20.0; // Thickness of each sandstone unit 


// Points
// Base of Domain

Point(1) = {0,0,-(2*sandstonethickness+saltthickness),h};
Point(2) = {x,0,-(2*sandstonethickness+saltthickness),h};
Point(3) = {x,y,-(2*sandstonethickness+saltthickness),h};
Point(4) = {0,y,-(2*sandstonethickness+saltthickness),h};

// Lines

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Line Loops: Used to help define surfaces, these connect lines
Line Loop(1) = {1, 2, 3, 4};

// Plane Surfaces: Define actual surfaces based on Line Loops
Plane Surface(1) = {1};

Physical Surface(1) = {1}; // Bottom Boundary

// Top of Bottom Sandstone Unit

Point(5) = {0,0,-(sandstonethickness+saltthickness),h};
Point(6) = {x,0,-(sandstonethickness+saltthickness),h};
Point(7) = {x,y,-(sandstonethickness+saltthickness),h};
Point(8) = {0,y,-(sandstonethickness+saltthickness),h};

// Lines

Line(5) = {1, 5};
Line(6) = {2, 6};
Line(7) = {3, 7};
Line(8) = {4, 8};

// Bottom of Anticline

Point(9) = {0.2*x,0,-(sandstonethickness),h};
Point(10) = {0.5*x,0,-(sandstonethickness-saltthickness),h};
Point(11) = {0.8*x,0,-(sandstonethickness),h};


Point(12) = {0.2*x,y,-(sandstonethickness),h};
Point(13) = {0.5*x,y,-(sandstonethickness-saltthickness),h};
Point(14) = {0.8*x,y,-(sandstonethickness),h};

//+
Bezier(9) = {5, 9, 10, 11, 6};
//+
Bezier(10) = {7, 14, 13, 12, 8};
//+
Line(11) = {8, 5};
//+
Line(12) = {6, 7};
//+
Curve Loop(2) = {10, -8, -3, 7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, 11, -5, -4};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {9, -6, -1, 5};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {6, 12, -7, -2};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {9, 12, 10, 11};
//+
Surface(6) = {6};
//+
Surface Loop(1) = {6, 4, 5, 2, 3, 1};
//+
Volume(1) = {1};
//+
Physical Surface("Salt Bottom Interface", 2) = {6};
//+
Physical Volume("Bottom Sandstone Volume", 6) = {1};
//+
Point(15) = {0,0,-(sandstonethickness),h};
//+
Point(16) = {0.2*x,0,-(sandstonethickness-saltthickness),h};
//+
Point(17) = {0.5*x,0,-(sandstonethickness-2*saltthickness),h};
//+
Point(18) = {0.8*x,0,-(sandstonethickness-saltthickness),h};
//+
Point(19) = {x,0,-(sandstonethickness),h};
//+
Bezier(13) = {15,16,17,18,19};
//+
Point(20) = {0,y,-(sandstonethickness),h};
//+
Point(21) = {0.2*x,y,-(sandstonethickness-saltthickness),h};
//+
Point(22) = {0.5*x,y,-(sandstonethickness-2*saltthickness),h};
//+
Point(23) = {0.8*x,y,-(sandstonethickness-saltthickness),h};
//+
Point(24) = {x,y,-(sandstonethickness),h};
//+
Bezier(14) = {20,21,22,23,24};
//+
Line(15) = {20, 15};
//+
Line(16) = {20, 8};
//+
Line(17) = {15, 5};
//+
Line(18) = {19, 24};
//+
Line(19) = {7, 24};
//+
Line(20) = {6, 19};
//+
Curve Loop(8) = {19, -14, 16, -10};
//+
Plane Surface(7) = {8};
//+
Curve Loop(9) = {11, -17, -15, 16};
//+
Plane Surface(8) = {9};
//+
Curve Loop(10) = {13, -20, -9, -17};
//+
Plane Surface(9) = {10};
//+
Curve Loop(11) = {20, 18, -19, -12};
//+
Plane Surface(10) = {11};
//+
Curve Loop(12) = {18, -14, 15, 13};
//+
Surface(11) = {12};
//+
Physical Surface("Top Sandstone Interface", 3) = {11};
//+
Surface Loop(2) = {11, 10, 9, 8, 7, 6};
//+
Curve Loop(14) = {10, -16, 14, -19};
//+
Surface(12) = {14};
//+
Curve Loop(16) = {13, -20, -9, -17};
//+
Surface(13) = {16};
//+
Surface Loop(3) = {11, 10, 8, 6, 7, 9};
//+
Volume(2) = {3};
//+
Physical Volume("Salt Volume", 7) = {2};

Point(25) = {0,0,0,h};
Point(26) = {x,0,0,h};
Point(27) = {x,y,0,h};
Point(28) = {0,y,0,h};
//+
Line(21) = {15, 25};
//+
Line(22) = {26, 19};
//+
Line(23) = {24, 27};
//+
Line(24) = {28, 20};
//+
Line(25) = {28, 27};
//+
Line(26) = {27, 26};
//+
Line(27) = {26, 25};
//+
Line(28) = {28, 25};
//+
Curve Loop(18) = {28, -27, -26, -25};
//+
Plane Surface(14) = {18};
//+
Physical Surface("Top Boundary", 4) = {14};
//+
Curve Loop(19) = {21, -27, 22, -13};
//+
Plane Surface(15) = {19};
//+
Curve Loop(20) = {28, -21, -15, -24};
//+
Plane Surface(16) = {20};
//+
Curve Loop(21) = {24, 14, 23, -25};
//+
Plane Surface(17) = {21};
//+
Curve Loop(22) = {23, 26, 22, 18};
//+
Plane Surface(18) = {22};
//+
Surface Loop(4) = {15, 16, 14, 18, 17, 11};
//+
Volume(3) = {4};
//+
Physical Volume("Top Sandstone Volume", 8) = {3};
//+
Curve Loop(23) = {21, -27, 22, -20, -6, -1, 5, -17};
//+
Curve Loop(24) = {21, -27, 22, -20, -6, -1, 5, -17};
//+
Plane Surface(19) = {24};
//+
Physical Surface("Side Boundary", 5) = {19, 5, 10, 18, 2, 7, 17, 16, 8, 3};
