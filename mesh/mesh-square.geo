a = 20;
b = 7;
c = 3;

h = 0.1; // Mesh size.

// Create one point in the origin.
Point(1) = {0, 0, 0, h};

// Extrude the point along x to create one side. The Layers option indicates the
// number of mesh subdivisions along the extrusion.
Extrude {a, 0, 0} { Point{1}; Layers{a / h}; }

// Extrude that side along y to create the square.
Extrude {0, b, 0} { Line{1};  Layers{b / h}; Recombine; }

// Extrude that side along y to create the square.
Extrude {0, 0, c} { Surface{5};  Layers{c / h}; Recombine; }

// Define the tags.
// Physical Line(0) = {3};
// Physical Line(1) = {4};
// Physical Line(2) = {1};
// Physical Line(3) = {2};

// Physical Surface(10) = {5};

// Generate a 2D mesh.
Mesh 3;

// Save mesh to file.
str_h = Sprintf("%f", h);
Save StrCat("../build/mesh/mesh-square-h", str_h, ".msh");
