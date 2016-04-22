// Gmsh project created on Tue May 19 13:43:01 2015
lc = 0.1;  //0.01
li = 0.1;
lb = 0.4;

xc3 = 1.0;
yc3 = 5.0;
R3  = 0.125;

xc2 = 2.0;
yc2 = 5.0;
R2  = 0.125;

xc1 = 1.0;
yc1 = 0.3;
R1  = 0.125;

xc4 = 4.0;
yc4 = 5.0;
R4  = 0.125;

xc5 = 5.0;
yc5 = 5.0;
R5  = 0.125;

//BOX
Point(1) = {0, 0, 0, li};
Point(2) = {2, 0, 0, li};
Point(3) = {2, 6, 0, lb};
Point(4) = {0, 6, 0, lb};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(9) = {4, 1, 2, 3};

//fluid physicals
Physical Line(1) = {4};  //left
Physical Line(2) = {1};  //bottom
Physical Line(3) = {2};  //right
Physical Line(4) = {3};  //superior
Physical Point(5) = {4, 3, 2, 1}; //box's corners



//CIRCLES
//C1
Point(15) = {xc1, yc1, 0, lc};
Point(16) = {xc1 + R1, yc1, 0, lc};
Point(17) = {xc1, yc1 + R1, 0, lc};
Point(18) = {xc1 - R1, yc1, 0, lc};
Point(19) = {xc1, yc1 - R1, 0, lc};
Circle(15) = {16, 15, 17};
Circle(16) = {17, 15, 18};
Circle(17) = {18, 15, 19};
Circle(18) = {19, 15, 16};
Line Loop(10) = {16, 17, 18, 15};
Physical Line(101) = {16, 17, 18, 15};  //interface fs
Physical Point(101) = {17, 16, 19, 18}; //interfacial points
/*//C2
Point(25) = {xc2, yc2, 0, lc};
Point(26) = {xc2 + R2, yc2, 0, lc};
Point(27) = {xc2, yc2 + R2, 0, lc};
Point(28) = {xc2 - R2, yc2, 0, lc};
Point(29) = {xc2, yc2 - R2, 0, lc};
Circle(25) = {26, 25, 27};
Circle(26) = {27, 25, 28};
Circle(27) = {28, 25, 29};
Circle(28) = {29, 25, 26};
Line Loop(20) = {26, 27, 28, 25};
Physical Line(102) = {26, 27, 28, 25};  //interface fs
Physical Point(102) = {27, 26, 29, 28}; //interfacial points
//C3
Point(35) = {xc3, yc3, 0, lc};
Point(36) = {xc3 + R3, yc3, 0, lc};
Point(37) = {xc3, yc3 + R3, 0, lc};
Point(38) = {xc3 - R3, yc3, 0, lc};
Point(39) = {xc3, yc3 - R3, 0, lc};
Circle(35) = {36, 35, 37};
Circle(36) = {37, 35, 38};
Circle(37) = {38, 35, 39};
Circle(38) = {39, 35, 36};
Line Loop(30) = {36, 37, 38, 35};
Physical Line(103) = {36, 37, 38, 35};  //interface fs
Physical Point(103) = {37, 36, 39, 38}; //interfacial points
//C4
Point(45) = {xc4, yc4, 0, lc};
Point(46) = {xc4 + R4, yc4, 0, lc};
Point(47) = {xc4, yc4 + R4, 0, lc};
Point(48) = {xc4 - R4, yc4, 0, lc};
Point(49) = {xc4, yc4 - R4, 0, lc};
Circle(45) = {46, 45, 47};
Circle(46) = {47, 45, 48};
Circle(47) = {48, 45, 49};
Circle(48) = {49, 45, 46};
Line Loop(40) = {46, 47, 48, 45};
Physical Line(104) = {46, 47, 48, 45};  //interface fs
Physical Point(104) = {47, 46, 49, 48}; //interfacial points
//C5
Point(55) = {xc5, yc5, 0, lc};
Point(56) = {xc5 + R5, yc5, 0, lc};
Point(57) = {xc5, yc5 + R5, 0, lc};
Point(58) = {xc5 - R5, yc5, 0, lc};
Point(59) = {xc5, yc5 - R5, 0, lc};
Circle(55) = {56, 55, 57};
Circle(56) = {57, 55, 58};
Circle(57) = {58, 55, 59};
Circle(58) = {59, 55, 56};
Line Loop(50) = {56, 57, 58, 55};
Physical Line(105) = {56, 57, 58, 55};  //interface fs
Physical Point(105) = {57, 56, 59, 58}; //interfacial points
*/

Plane Surface(11) = {9, 10};
Physical Surface(10) = {11};       //fluid

//Plane Surface(12) = {10};
//Physical Surface(16) = {12};
//Transfinite Line {5, 6, 7, 8} = 20;// Using Bump 0.2;
//Transfinite Line {1, 2, 3, 4} = 20;
//Transfinite Surface {12};
//Transfinite Surface {11};
//Recombine Surface {12};
//Recombine Surface {11};
