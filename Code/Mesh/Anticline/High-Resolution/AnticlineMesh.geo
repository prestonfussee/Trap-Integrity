h = 2.0; // Element scale factor >> 1 is coarse
x = 100.0; 
y = 100.0;
saltthickness = 10.0; // Thickness of the salt unit
sandstonethickness = 20.0; // Thickness of each sandstone unit 


// Points
// Base of Domain

Point(1) = {0,0,-(2*sandstonethickness+saltthickness),h}; //
Point(2) = {x,0,-(2*sandstonethickness+saltthickness),h}; //
Point(3) = {x,y,-2*sandstonethickness-saltthickness, h};
Point(4) = {0,y,-(2*sandstonethickness+saltthickness),h};
Point(5) = {0,0,-(sandstonethickness+saltthickness),h};
Point(6) = {x,0,-(sandstonethickness+saltthickness),h};
Point(7) = {x,y,-(sandstonethickness+saltthickness),h};
Point(8) = {0,y,-(sandstonethickness+saltthickness),h};
Point(9) = {0.35*x,0,-(sandstonethickness),h};
//Point(10) = {0.5*x,0,-sandstonethickness+saltthickness, h};
Point(11) = {0.65*x,0,-(sandstonethickness),h};
Point(12) = {0.35*x,y,-(sandstonethickness),h};
//Point(13) = {0.5*x,y,-(sandstonethickness-saltthickness),h};
Point(14) = {0.65*x,y,-(sandstonethickness),h};
Point(15) = {0,0,-(sandstonethickness),h};
Point(16) = {0.35*x,0,-(sandstonethickness-saltthickness),h};
//Point(17) = {0.5*x,0,-(sandstonethickness-2*saltthickness),h};
Point(18) = {0.65*x,0,-(sandstonethickness-saltthickness),h};
Point(19) = {x,0,-(sandstonethickness),h}; 
Point(20) = {0,y,-sandstonethickness, h};
Point(21) = {0.35*x,y,-(sandstonethickness-saltthickness),h};
//Point(22) = {0.5*x,y,-(sandstonethickness-2*saltthickness),h};
Point(23) = {0.65*x,y,-(sandstonethickness-saltthickness),h};
Point(24) = {x,y,-(sandstonethickness),h};
Point(25) = {0,0,0,h};
Point(26) = {x,0,0,h};
Point(27) = {x,y,0,h};
Point(28) = {0,y,0,h};

// Lines

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {1, 5};
Line(6) = {2, 6};
Line(7) = {3, 7};
Line(8) = {4, 8};
Bezier(9) = {5,9,11,6};
Bezier(10) = {7, 14, 12, 8};
Line(11) = {8, 5};
Line(12) = {6, 7};
Bezier(13) = {15,16,18,19};
Bezier(14) = {20,21,23,24};
Line(15) = {20, 15};
Line(16) = {20, 8};
Line(17) = {15, 5};
Line(18) = {19, 24};
Line(19) = {7, 24};
Line(20) = {6, 19};
Line(21) = {15, 25};
Line(22) = {26, 19};
Line(23) = {24, 27};
Line(24) = {28, 20};
Line(25) = {28, 27};
Line(26) = {27, 26};
Line(27) = {26, 25};
Line(28) = {28, 25};

// Line Loops: Used to help define surfaces, these connect lines
Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {10, -8, -3, 7};
Line Loop(3) = {8, 11, -5, -4};
Line Loop(4) = {9, -6, -1, 5};
Line Loop(5) = {6, 12, -7, -2};
Line Loop(6) = {9, 12, 10, 11};
Line Loop(8) = {19, -14, 16, -10};
Line Loop(9) = {11, -17, -15, 16};
Line Loop(10) = {13, -20, -9, -17};
Line Loop(11) = {20, 18, -19, -12};
Line Loop(12) = {18, -14, 15, 13};
Line Loop(14) = {10, -16, 14, -19};
Line Loop(16) = {13, -20, -9, -17};
Line Loop(18) = {28, -27, -26, -25};
Line Loop(19) = {21, -27, 22, -13};
Line Loop(20) = {28, -21, -15, -24};
Line Loop(21) = {24, 14, 23, -25};
Line Loop(22) = {23, 26, 22, 18};
Line Loop(24) = {21, -27, 22, -20, -6, -1, 5, -17};
Line Loop(25) = {2, 7, 19, 23, 26, 22, -20, -6};
Line Loop(26) = {25, -23, -19, -7, 3, 8, -16, -24};
Line Loop(27) = {24, 16, -8, 4, 5, -17, 21, -28};


// Plane Surfaces: Define actual surfaces based on Line Loops
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Surface(6) = {6};
Plane Surface(7) = {8};
Plane Surface(8) = {9};
Plane Surface(9) = {10};
Plane Surface(10) = {11};
Surface(11) = {12};
Surface(12) = {14};
Surface(13) = {16};
Plane Surface(14) = {18};
Plane Surface(15) = {19};
Plane Surface(16) = {20};
Plane Surface(17) = {21};
Plane Surface(18) = {22};
Plane Surface(19) = {24};
Plane Surface(20) = {25};
Plane Surface(21) = {26};
Plane Surface(22) = {27};

// Surface Loops: Used to help define volumes
Surface Loop(1) = {6, 4, 5, 2, 3, 1}; // Bottom Sandstone Volume
Surface Loop(2) = {11, 10, 9, 8, 7, 6};
Surface Loop(3) = {15, 16, 14, 18, 17, 11};


// Volumes: Defined from Surface Loops
Volume(1) = {1};
Volume(2) = {2};
Volume(3) = {3};

// Physical Groups
// 1: Bottom Boundary
// 2: Salt Bottom Interface
// 3: Salt Top Interface
// 4: Top Boundary
// 5: Side Boundaries (all)
// 6: Bottom Sandstone Volume
// 7: Salt Volume
// 8: Top Sandstone Volume

Physical Surface(1) = {1}; // Bottom Boundary
Physical Surface(2) = {6}; // "Salt Bottom Interface" 
Physical Surface(3) = {11}; // "Top Sandstone Interface"
Physical Surface(4) = {14}; // "Top Boundary"
Physical Surface(5) = {19,20,21,22}; // "Side Boundary"
Physical Volume(6) = {1}; // "Bottom Sandstone Volume"
Physical Volume(7) = {2}; // "Salt Volume"
Physical Volume(8) = {3}; // "Top Sandstone Volume"


