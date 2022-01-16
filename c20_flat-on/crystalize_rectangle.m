
%space between alkanes in x direction (Ang)
spc_x = 2; 

%space between alkanes in y direction (Ang)
spc_y = 4.45; %flat-on
%spc_y = 3.80; %edge-on

%length of the alkane chain (Ang)
l = 25.51;

%number of atoms in one alkane
alk_atoms = 62;

%rhombic simulation cell for graphene
a = 2.46; %graphene lat. param. (Ang)
N = 108; %number of graphene unit cells
S1 = N*a*[1; 0]; %sim. cell vector   
S2 = N*a*[0.5; sqrt(3)/2]; %sim. cell vector

%rectangle area of the alkane crystal
T = (S1 + S2)/4; %left corner (origin)
D1 = (norm(S1) - norm(T)*norm(S1)/(norm(S1+S2)))*[1;0]; %cell vector #1 of the rectangle
D2 = norm(T)*norm(S1)/(norm(S1+S2))*tand(60)*[0;1]; %cell vector #2 of the rectangle

%load alkane structure
alkane = dlmread('c20_eicosane_flat-on_v2.lammps', ' ', 9, 0); %skip first 9 lines
%alkane = dlmread('c20_eicosane_edge-on_v2.lammps', ' ', 9, 0); %skip first 9 lines

%main loop to place alcaines
%================================
%unit cell vectors of the alkane crystalline phase
dx = (l+spc_x)*[1; 0];
dy = spc_y*[0; 1];

Nrows = floor(norm(D2)/norm(dy)); %number of alkane rows in sim. cell
Ncols = floor(norm(D1)/norm(dx)); %number of alkanes columns in sim. cell
tot_atoms = Nrows*Ncols*alk_atoms; %total number of atoms

full_structure = zeros(tot_atoms,5); %store alkane crystal (id type x y z)

n1 = 1;
n2 = alk_atoms;

for i = 1:Nrows
    for j = 1:Ncols

        %add one alkane to layer
        full_structure(n1:n2,:) = alkane(:,:); 

        %shift alkane to crystal position
        full_structure(n1:n2,3) = full_structure(n1:n2,3) + (j-1)*dx(1) + (i-1)*dy(1) + T(1);
        full_structure(n1:n2,4) = full_structure(n1:n2,4) + (i-1)*dy(2) + T(2);

        %shift atom id-s
        full_structure(n1:n2,1) = full_structure(n1:n2,1) + (j-1 + (i-1)*Ncols)*alk_atoms;

        n1 = n1+alk_atoms;
        n2 = n2+alk_atoms; 
        
    end
end

%a final shift of the alkane layer to fine tune position
% remx = rem(norm(S1),norm(dx));
% remy = rem(norm(S2),norm(dy));
% full_structure(:,3) = full_structure(:,3) + remx;
% full_structure(:,4) = full_structure(:,4) + remy;
%=================================

%output file to lammps data format
filename = 'alkane_crystal.lammps';
fid = fopen(filename, 'wt');

%header
fprintf(fid, '#\n');
fprintf(fid, '%d atoms\n', tot_atoms);
fprintf(fid, '5 atom types\n');
fprintf(fid, '0.0 %f xlo xhi\n',S1(1));
fprintf(fid, '0.0 %f ylo yhi\n',S2(2));
fprintf(fid, '0.0 1.78 zlo zhi\n');
fprintf(fid, '%f 0.0 0.0 xy xz yz\n\n',S2(1));
fprintf(fid, 'Atoms # atomic\n\n');

for q=1:size(full_structure,1)
    fprintf(fid, '%d %d %6.3f %6.3f %6.3f\n', full_structure(q,1), full_structure(q,2), full_structure(q,3), full_structure(q,4), full_structure(q,5));
end

fclose(fid);
