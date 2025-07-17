function c2_boundary =initialize_boundaries()

%Boundaries are defined as a set of contours e.g. [x1 NaN x2 ; y1 NaN y2]

%%Kobuk Upper Lake
Lx = 3e4; Ly = 4e4; 
x=[-1 -1 1 1]*Lx; 
y=[-1 1 1 -1]*Ly;

c2_boundary = [x; y]; %defines box around the boundaries

end