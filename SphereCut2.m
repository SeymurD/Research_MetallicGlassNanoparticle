%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Dean's Research Award Project
% Name: Seymur Dadashov
% CCID: seymur
% U of A ID: 1495584
%
%Professor: Dr. Hao Zhang
%
% Description: This program is designed to analyze a 'dump' file consisting
% of atom coordinates within a certain fixed volume cell. The program
% organizes atoms by their ID into a matrix. The matrix is then used to
% determine which atoms satisfy a sphere equation within the fixed volume
% cell such that the ratio of atom of type 1 to type 2 meets a certain
% criteria. Ex. Cu64Zr36 where the percent of atoms of Copper must be
% identicall 64% and Zirconium must be identically 36%. If a solution of
% such as sphere is not found with inputed radius, then the program will
% switch atom types (i.e. type 1 turns to type 2 or type 2 turn type 1 
% until the closest possible ratio to 64/36 is achieved. The program will
% then graph the atoms onto a scatter plot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

%Opening the target file and extracting the atom index, type and x,y,z
%coordinates into a cell matrix. The cell matrix is converted into a
%standard matrix for convenience.
fid = fopen('coords.dat', 'r');
startRow = 15;                      %Row where the data on the first atom starts.
endRow = 13514;

formatSpec = '%f%f%f%f%f%[^\n\r]';
textscan(fid, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fid, formatSpec, endRow-startRow+1);

dataArray{6} = [];
C = cell2mat(dataArray);

%The matrix is organized based on the atom ID number going from 1 to 13500.
[~, idx] = unique(C(:,1));
sD = C(idx, :);

%The logical array check for which indeces of the matrix the type of atoms
%is Type 1 i.e. Copper.
Cu_check = sD(:,2)==1;

%Creating individual arrays for x, y, and z coordinates for convenience in
%later use.
x = sD(:,3);
y = sD(:,4);
z = sD(:,5);

CuGraph = [];
ZrGraph = [];
for i=1:length(sD)
    if Cu_check(i) == 1
        CuGraph(i,:) = sD(i,3:5);
    else
        CuGraph(i,:) = [0 0 0];
    end
    if Cu_check(i) ~= 1
        ZrGraph(i,:) = sD(i,3:5);
    else
        ZrGraph(i,:) = [0 0 0];
    end
end

figure(1)
scatter3(CuGraph(:,1), CuGraph(:,2), CuGraph(:,3), 'blue');
hold on;
scatter3(ZrGraph(:,1), ZrGraph(:,2), ZrGraph(:,3), 'red');

%Initializing logical variables to keep the while loop below activated.
possibility = true; 
validInput = true;

while possibility == true
    %User is prompted to enter a max value for boundary of the fixed volume cell.
    bLimit = input('Enter the maximum value for the boundary of the model: \n');  
    %Initializing variables.
    newMin = 51;
    minShape = [];
    %Check for user input to make sure that the radius inputed does not
    %exceed boundary limits, i.e. sphere is not big enough to leave the
    %fixed volume cell.
    while validInput == 1
        newR = input('Enter radius (such that the diameter does not exceed boundary limits): \n');
        if 2*newR > bLimit
            fprintf('Boundary of cut atom exceeds dimensions, choose smaller radius \n');
        else
            break;
        end
    end
   
   %The vec_a,b,c are designed to shift a sphere throughout the fixed
   %volume cell. The three nested for loops are responsible for making sure
   %that the sphere equation is shifted for all possible values of
   %vec_a,b,c where initially the sphere equation is shifted for all values
   %of z for every y for every x. Therefore, as each new equation is formed
   %with a slightly different shift:
   %(x(i)-newR-a)^2 + (y(i)-newR-b)^2 + (z(i)-newR-c)^2 <= newR^2
   %new possibilities arise for the possible number of atoms found in each
   %sphere. In the equation newR shifts the sphere from origin completely
   %into the positive x,y,z quadrant. Then a, b, and c shift the sphere
   %within the positive quadrant which is the fixed volume cell. The loop
   %will complete every possible shift with the given radius such that is
   %does not leave the fixed volume cell. The minimum difference between
   %Type 1 and Type 2 that results in the closest combination to the
   %required ratio i.e. the closest number of Cu and Zr is taken to the
   %ratio Cu64Zr36. If an identical ratio exists, then the difference is 0.
    vec_a = 0:(bLimit-2*newR);
    vec_b = 0:(bLimit-2*newR);
    vec_c = 0:(bLimit-2*newR); 
    for m=1:length(vec_a)
        a = vec_a(1,m);
        for j=1:length(vec_b)
            b = vec_b(1,j);
            for k=1:length(vec_c)
                c = vec_c(1,k);
                for i=1:length(sD)
                    if (x(i)-newR-a)^2 + (y(i)-newR-b)^2 + (z(i)-newR-c)^2 <= newR^2
                        cutShape(i,:) = [x(i) y(i) z(i)]; %Storing data that fits inside sphere into coordinate matrix.
                    else
                        cutShape(i,:) = [0 0 0]; %If the values are outside sphere, assigned value of 0.
                    end
                end
                    
                    %Identifying which atoms are Cu and which are Zr.
                    cutIndex = cutShape(:,1)~=0;
                    countCu=[];
                    countZr=[];
                    for i=1:length(sD)
                        if cutIndex(i)~=0
                            if sD(i,2)==1
                                countCu(i) = 1;
                            else
                                countCu(i) = 0;
                            end

                            if sD(i,2)==2
                                countZr(i) = 1;
                            else
                                countZr(i) = 0;
                            end
                        end
                    end

                    AT = nnz(cutIndex(:,1));
                    CuT = nnz(countCu);
                    ZrT = nnz(countZr);

                    idealCu = (64/100)*AT;
                    idealZr = (36/100)*AT;

                    diffCu = idealCu - CuT;
                    diffZr = idealZr - ZrT;
                    
                    results(m, j, k) = diffCu; %If the difference between ideal and current is 0, the ideal ratio exists.

                    %The if statement collects all relevant information for
                    %the closest current sphere to the requested ratio.w
                    if abs(diffCu)<=50 && abs(diffCu)<newMin
                        newMin = max(diffCu,diffZr);
                        minAT=AT;
                        minCuT=CuT;
                        minZrT=ZrT;
                        minShape= cutShape;    
                    end
            end
        end
    end
    
    %Locates the minimum difference value found after all combinations are
    %ran.
    [minCu, minIDX] = min(abs(results(:)));
    [row,col,page]=ind2sub(size(results),minIDX);
    
    %Prints location of the combination used.
    disp(row);disp(col);disp(page);
    
    %Prints the minimum value.
    minCu = results(minIDX);
    disp(minCu);
    
    %Prints current total number of atoms, the number of Copper atoms, and
    %the number of Zirconium atoms.
    fprintf('AT: %5f, Cu: %5f, Zr: %5f \n', minAT, minCuT, minZrT);
   
    %Prints whether or not Copper atoms should be removed because they
    %exceed the ideal number of Copper atoms or if Zirconium atoms should
    %be removed if they exceed ideal number of Zirconium atoms.
    if minCu < 0
        fprintf('The minimum number of Copper atoms to be removed and replaced with Zirconium is: %5f. (ideal Cu atoms < current Cu atoms) \n', abs(minCu));
    else
        fprintf('The minimum number of Zirconium atoms to be removed and replaced with Copper is: %5f. (ideal Zr atoms < current Zr atoms)\n', abs(minCu));
    end
    
    %Asks user if they would like to try a different radius or continue
    %with finishing the calculations.
    tryAgain = input('Would you like to continue with this or select new radius? 1(Continue) or 0(New Radius) \n');
    if tryAgain == 1
        break;
    end
end

%Appending the column for ID and type of atom to the matrix, since it is
%only a coordinate matrix as of now.
temp = minShape;
minShape = [];
minShape = [sD(:,1) sD(:,2) temp];

for i=1:length(minShape)
    if minShape(i,3) == 0
        minShape(i,2) = 0;
    end
end

%Switching Cu to Zr atoms since there are more Cu than ideal Cu atoms.
if  minCu~=0 && abs(minCu)>=0.05 && minCu<0
    rCu = abs(round(minCu));
    mix = randperm(size(minShape,1));
    s_minShape = minShape(mix,:);
    counter = 0;
    for i=1:length(minShape)
        if s_minShape(i,2) == 1 && counter < rCu
            s_idx = s_minShape(i,1);
            minShape(s_idx, 2) = 2;
            counter = counter + 1;
        end
    end
end

%Switching Zr to Cu since there are more Zr than ideal Zr atoms.
if minCu ~= 0 && abs(minCu)>=0.05 && minCu>0
    rCu = abs(round(minCu));
    mix = randperm(size(minShape,1));
    s_minShape = minShape(mix,:);
    counter = 0;
    for i=1:length(minShape)
        if s_minShape(i,2)==2 && counter < rCu
            s_idx = s_minShape(i,1);
            minShape(s_idx, 2) = 1;
            counter = counter + 1;
        end
    end
end

%Recalculating the indeces of Cu and Zr atoms and then proceeding to count
%the total number of atoms, the adjusted ('a_') number of Copper atoms, and
%the adusted number of Zirconium atoms.
a_countCu=[];
a_countZr=[];
a_cutIndex = minShape(:,3)~=0;
for i=1:length(minShape)
    if a_cutIndex(i)~=0
        if minShape(i,2)==1
            a_countCu(i) = 1;
        else
            a_countCu(i) = 0;
        end

        if minShape(i,2)==2
            a_countZr(i) = 1;
        else
            a_countZr(i) = 0;
        end
    end
end

a_AT = nnz(minShape(:,3));
a_CuT = nnz(a_countCu);
a_ZrT = nnz(a_countZr);

a_idealCu = (64/100)*a_AT;
a_idealZr = (36/100)*a_AT;

a_diffCu = a_idealCu - a_CuT;
a_diffZr = a_idealZr - a_ZrT;

%Printing new values for total number of atoms, number of Copper atoms and
%the number of Zirconium atoms.
fprintf('AT: %5f, Cu: %5f, Zr: %5f DiffCu: %5f, DiffZr: %5f \n', a_AT, a_CuT, a_ZrT, a_diffCu, a_diffZr);                    

%Plotting a scatter plot of all the pots with the ideal or closest possible
%ratio to 64/36.

Cu_minCheck = minShape(:,2)==1;

Cu_minGraph = [];
Zr_minGraph = [];
for i=1:length(minShape)
    if Cu_minCheck(i) == 1
        Cu_minGraph(i,:) = minShape(i,3:5);
    else
        Cu_minGraph(i,:) = [0 0 0];
    end
    if Cu_minCheck(i) ~= 1
        Zr_minGraph(i,:) = minShape(i,3:5);
    else
        Zr_minGraph(i,:) = [0 0 0];
    end
end

figure(2)
scatter3(Cu_minGraph(:,1), Cu_minGraph(:,2), Cu_minGraph(:,3), 'red');
hold on;
scatter3(Zr_minGraph(:,1), Zr_minGraph(:,2), Zr_minGraph(:,3), 'blue');




