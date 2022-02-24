%Joan Akibode's FEA Project

%Coarse mesh with circulat obstacle

%UNI : jga2129


%importing the files containing the points coordinates 

coordinates = readmatrix("rounded_coarse_node_coordinate.txt");
connectivity = readmatrix("rounded_coarse_node_connectivity.txt");

v = coordinates(:,2:3);
f = connectivity(:,2:4);

%ploting of the mesh

figure

patch("Faces",f,"Vertices",v,"FaceColor",'white')
[t,s] = title("Coarse mesh with rounded obstacle","q = 0");
text(11, 0,'T = 0')
text(-11.5, 0,'qx = -1')
text(0, -6, 'q = 0')
xlabel("X")
ylabel("Y")

% retriving the size of the coordiantes and connectivity matrice sizes

sizef = size(f);
n = sizef(1); 
p = sizef(2);

sizem = size(v);
m = sizem(1);

%Creation of the matrices and variables

B = zeros(2*m,3);       % will contain the gradients of local shape functions
K = zeros(m);           % will contain the global conductance matrix
F = zeros(m,1);         % will countain the global boundary flux matrix
T = zeros(m,1);         % Will countain the emperatures
flux = zeros(n,2);      % Will countain the coordinates of the flux vectors
centroids = zeros(n,2); % Will countain the coordinates of the triangles centroîds

maxX = max(v(:,1)); % Maximum abscisse of the points (where T = 0 )
minX = min(v(:,1)); % Minimum abcsisse of the points (Where q = -1)


% construction of the global conductance Matrix K and the global boundary
% flux matrix F

syms y
  
for i = 1:n % We have to construct the local F and K of the n triangles

    % store index of each point
    
    numPoint = [f(i,1) f(i,2) f(i,3)]; 

    % coordinates of triangle's points

    x1 = v(numPoint(1),1);
    y1 = v(numPoint(1),2);
    x2 = v(numPoint(2),1);
    y2 = v(numPoint(2),2);
    x3 = v(numPoint(3),1);
    y3 = v(numPoint(3),2);

    % Calculation of the centroids that we will use to plot the flux

    centroids(i,:) = [(x1+x2+x3)/3, (y1+y2+y3)/3]; 
    
 
    % Surface of the triangle multiplied by 2 :

    Ae2 = (x2 * y3-x3*y2) - (x1*y3-x3*y1) + (x1*y2-x2*y1);
    
    
    %contruction of the gradient of local shape functions Bloc :

    Bloc = (1/Ae2)*[(y2-y3),(y3-y1),(y1-y2);(x3-x2),(x1-x3),(x2-x1)];
    
    % Storing the Bloc matrix by addind it at the end on B
    B(1+(i-1)*2 : 2+(i-1)*2,:) = Bloc;

    Kloc = (1/2)*Ae2 * (Bloc.')*Bloc;

    % Add the value of Kloc to the global K:

    for j = 1:3
        for k = 1:3

            K(numPoint(j),numPoint(k)) = K(numPoint(j),numPoint(k)) + Kloc(j,k);

        end
    end
    
    %contruction of F :

    %Is one of the segment of the triangle on the border where q = -1 ?

    nbPointsOnBorder = 0;           % will countain the number of point on the right border where q = -1
    yPointsOnBorder = zeros(2,1);   % will countain the y coordinates if those points

    for k = 1:3 % We go through the 3 points of the triangle

        if v(numPoint(k),1) <= minX % is the point on the left border ?
            
            % if yes we increment and take its y coordinate

            nbPointsOnBorder = nbPointsOnBorder + 1; 

            yPointsOnBorder(nbPointsOnBorder) = v(numPoint(k),2); 
        end
    end
    
    % if two points have y = -5 then one of the segment of the triangle is
    % on the border where q = -1

    if nbPointsOnBorder == 2


        yPointsOnBorder = sort(yPointsOnBorder);
           
        % Calculation of the primitives F1, F2 and F3 of N1, N2 and N3
        

        F1(y) = (1/Ae2)*( ((y2-y3)*(minX)+ (x2*y3-x3*y2))*y + (x3-x2)*(y*y)/2);
        F2(y) = (1/Ae2)*(((y3-y1)*(minX)+(x3*y1-x1*y3))*y + (x1-x3)* (y*y)/2);
        F3(y) = (1/Ae2)*(((y1-y2)*(minX)+ (x1*y2-x2*y1))*y + (x2-x1)*(y*y)/2);
        


        % F local

        Floc = [F1(yPointsOnBorder(2))-F1(yPointsOnBorder(1)), F2(yPointsOnBorder(2))-F2(yPointsOnBorder(1)), F3(yPointsOnBorder(2))-F3(yPointsOnBorder(1))];
        
        %F global : we increment the values in the local matrix to the
        %global matrix

        F(numPoint(1)) = F(numPoint(1)) + Floc(1);
        F(numPoint(2)) = F(numPoint(2)) + Floc(2);
        F(numPoint(3)) = F(numPoint(3)) + Floc(3);

    end     

end

% Calculation of the temperatures
% on the right border the temperature is equal to 0. To calculate the
% temperature of the other points we have to firs eliminate the lines and
% columns corresponding to the points where the temperature is 0. We will
% call this matrices and vectors bis. the we can do Kbis * Tbis = Fbis eq
% to Tbis = inv(Kbis)*Fbis. We can then inject the temperatures from Tbis
% to T in theiur right positions

vbis = coordinates(:,2:3);
Kbis = K;
Fbis = F;

Tbis = [transpose(1:1:m),zeros(m,1)];% the first column of this matrix is 
% the index of the temperatures so that when I delete one of the rows I
% still now what was its original index

%pulledBack =[];

 for i = m:-1:1 % We pull back the rows of K,F and T corresponding to where T = 0.
     % The for loop is descending so that deleting rows doesnt impact the 
     % index of the rows before the one being deleted
        if v(i,1) == maxX
            vbis(i,:) = []; % this one is to check if we eliminate the good ones
            Tbis(i,:) = [];
            Fbis(i) = [];
            Kbis(i,:) = [];
            Kbis(:,i) = [];
            %pulledBack(end+1) = i;
        end
 end 



%Tbis(:,2) = inv(Kbis)*Fbis; % Unefficent way to calculate the temperatures
%according to matlab

Tbis(:,2) = Kbis\Fbis; % Calculation of the temperatures
 

sizeBis = size(Tbis);

for i = 1:sizeBis(1) % Injecting the calculated temperatures in the original T vectors
    T(Tbis(i,1)) = T(Tbis(i,1)) + Tbis(i,2); % Tbis(i,1) is the original index of the temperature and Tbis(i,2) is the calculated temp
end

%Plotting of the temperature

figure
patch("Faces",f,"Vertices",v,"FaceVertexCData",T,"FaceColor","interp");
title("Temperature plot of the coarse mesh with rounded obstacle")
xlabel("X")
ylabel("Y")
colorbar


%Claculation of the fluxes
for i=1:n
    flux(i,:) = transpose(-1*B(1+(i-1)*2:2+(i-1)*2,:)*[T(f(i,1));T(f(i,2));T(f(i,3))]);
end

%Plotting of the fluxes

figure
quiver(centroids(:,1),centroids(:,2),flux(:,1),flux(:,2));
title("Flux plot of the coarse mesh with rounded obstacle")
xlabel("X")
ylabel("Y")

