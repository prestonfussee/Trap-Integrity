h = 6.0; // Element scale factor >> 1 is coarse
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

// Top of Bottom Sandstone Unit

Point(5) = {0,0,-(sandstonethickness+saltthickness),0.25*h};
Point(6) = {x,0,-(sandstonethickness+saltthickness),0.25*h};
Point(7) = {x,y,-(sandstonethickness+saltthickness),0.25*h};
Point(8) = {0,y,-(sandstonethickness+saltthickness),0.25*h};

// Top of Salt Unit
Point(9) = {0,0,-sandstonethickness,0.25*h};
Point(10) = {x,0,-sandstonethickness,0.25*h};
Point(11) = {x,y,-sandstonethickness,0.25*h};
Point(12) = {0,y,-sandstonethickness,0.25*h};

// Top of Top Sandstone Unit
Point(13) = {0,0,0,h};
Point(14) = {x,0,0,h};
Point(15) = {x,y,0,h};
Point(16) = {0,y,0,h};

// Lines

Line(1) = {1, 2};
Line(2) = {6, 2};
Line(3) = {5, 6};
Line(4) = {5, 1};
Line(5) = {6, 7};
Line(6) = {7, 3};
Line(7) = {2, 3};
Line(8) = {7, 8};
Line(9) = {4, 8};
Line(10) = {3, 4};
Line(11) = {8, 5};
Line(12) = {4, 1};
Line(13) = {9, 12};
Line(14) = {12, 11};
Line(15) = {11, 7};
Line(16) = {12, 8};
Line(17) = {11, 10};
Line(18) = {10, 6};
Line(19) = {9, 10};
Line(20) = {9, 5};
Line(21) = {13, 14};
Line(22) = {14, 10};
Line(23) = {14, 15};
Line(24) = {15, 11};
Line(25) = {16, 15};
Line(26) = {16, 12};
Line(27) = {13, 16};
Line(28) = {13, 9};

// Line Loops: Used to help define surfaces, these connect lines
Line Loop(1) = {10, 12, 1, 7};
Line Loop(2) = {8, 11, 3, 5};
Line Loop(3) = {8, -9, -10, -6};
Line Loop(4) = {7, -6, -5, 2};
Line Loop(5) = {3, 2, -1, -4};
Line Loop(6) = {11, 4, -12, 9};
Line Loop(7) = {19, -17, -14, -13};
Line Loop(8) = {16, 11, -20, 13};
Line Loop(9) = {14, 15, 8, -16};
Line Loop(10) = {15, -5, -18, -17};
Line Loop(11) = {18, -3, -20, 19};
Line Loop(12) = {25, 24, -14, -26};
Line Loop(13) = {13, -26, -27, 28};
Line Loop(14) = {28, 19, -22, -21};
Line Loop(15) = {22, -17, -24, -23};
Line Loop(16) = {25, -23, -21, 27};

// Plane Surfaces: Define actual surfaces based on Line Loops
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(13) = {13};
Plane Surface(14) = {14};
Plane Surface(15) = {15};
Plane Surface(16) = {16};

// Surface Loops: Define Volumes from Plane Surfaces
Surface Loop(1) = {5, 2, 3, 6, 1, 4};
Surface Loop(2) = {7, 11, 10, 9, 8, 2};
Surface Loop(3) = {14, 13, 12, 16, 15, 7};

// Volumes: Define Volumes from Surface Loops
Volume(1) = {1};
Volume(2) = {2};
Volume(3) = {3};


// Physical Surfaces: Used to label parts of the domain from Plane Surface data
// Note: To be compatible with FEniCS, labeling must be only be by integer values with index beginning at 0

// Bottom Boundary = 1
// Bottom Salt Interface = 2
// Top Salt Interface = 3
// Top Boundary = 4
// Side Boundary = 5
// Bottom Sandstone Volume = 6
// Salt Volume = 7
// Top Sandstone Volume = 8

Physical Surface(1) = {1}; // Bottom Boundary
Physical Surface(2) = {2}; // Bottom Salt Interface
Physical Surface(3) = {7}; // Top Salt Interface
Physical Surface(4) = {16}; // Top Boundary
Physical Surface(5) = {14, 11, 5, 4, 10, 15, 6, 13, 8, 3, 9, 12}; // Side Boundaries
Physical Volume(6) = {1}; // Bottom Sandstone Volume
Physical Volume(7) = {2}; // Salt Volume
Physical Volume(8) = {3}; // Top Sandstone Volume
//+
Field[1] = MathEval;
//+
Field[1].F = "0.1";
//+
Field[2] = Restrict;
//+
Field[2].FacesList = {2, 3};
//+
Field[2].IField = 0;
//+
Field[3] = MathEval;
//+
Field[3].F = "0.1";
//+
Field[4] = Restrict;
//+
Field[4].FacesList = {2, 3};
//+
Field[4].IField = 1;
//+
Field[4].RegionsList = {6};
