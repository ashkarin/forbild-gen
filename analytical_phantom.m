%{

Yu, Z., Noo, F., Dennerlein, F., Wunderlich, A., Lauritsch, G., & 
Hornegger, J. (2012). Simulation tools for two-dimensional experiments in 
x-ray computed tomography using the FORBILD head phantom. 
Physics in medicine and biology, 57(13), N237.

%}

function [ phantom ] = analytical_phantom(oL, oR)
    if (oL~=1) 
        oL=0;
    end;
    if (oR~=1)
        oR=0;
    end;
    
    sha = 0.2*sqrt(3);
    y016b = -14.294530834372887;
    a16b = 0.443194085308632;
    b16b = 3.892760834372886;
    
    E = [-4.7 4.3 1.79989 1.79989 0 0.010 0; %1
          4.7 4.3 1.79989 1.79989 0 0.010 0; %2
        -1.08 -9 0.4 0.4 0 0.0025 0; %3
         1.08 -9 0.4 0.4 0 -0.0025 0; %4
         0 0 9.6 12 0 1.800 0; %5
         0 8.4 1.8 3.0 0 -1.050 0; %7
         1.9 5.4 0.41633 1.17425 -31.07698 0.750 0; %8
        -1.9 5.4 0.41633 1.17425 31.07698 0.750 0; %9
        -4.3 6.8 1.8 0.24 -30 0.750 0; %10
         4.3 6.8 1.8 0.24 30 0.750 0; %11
         0 -3.6 1.8 3.6 0 -0.005 0; %12
         6.39395 -6.39395 1.2 0.42 58.1 0.005 0; %13
         0 3.6 2 2 0 0.750 4; %14
         0 9.6 1.8 3.0 0 1.800 4; %15
         0 0 9.0 11.4 0 0.750 3; %16a
         0 y016b a16b b16b 0 0.750 1; %16b
         0 0 9.0 11.4 0 -0.750 oR; %6
         9.1 0 4.2 1.8 0 0.750 1];%R_ear

    %generate the air cavities in the right ear
    cavity1 = transpose(8.8:-0.4:5.6);
    cavity2 = zeros(9,1);
    cavity3_7 = ones(53,1)*[0.15 0.15 0 -1.800 0];
    for j = 1:3 
        kj = 8-2*floor(j/3); 
        dj = 0.2*mod(j,2);
        cavity1 = [cavity1; cavity1(1:kj)-dj; cavity1(1:kj)-dj];
        cavity2 = [cavity2; j*sha*ones(kj,1); -j*sha*ones(kj,1)];
    end
    E_cavity = [cavity1 cavity2 cavity3_7];

    %generate the left ear (resolution pattern)
    x0 = -7.0;
    y0 = -1.0;
    d0_xy = 0.04;

    d_xy = [0.0357, 0.0312, 0.0278, 0.0250];
    x00 = zeros(0,0);
    y00 = zeros(0,0);
    ab = 0.5*ones(5,1)*d_xy;
    ab = ab(:)*ones(1,4);
    leftear4_7 = [ab(:) ab(:) ones(80,1)*[0 0.750 0]];
    for i = 1:4
        y00 = [y00; transpose(y0+(0:4)*2*d_xy(i))];
        x00 = [x00; (x0+2*(i-1)*d0_xy)*ones(5,1)];
    end
    x00 = x00*ones(1,4);
    y00 = [y00;y00+12*d0_xy; y00+24*d0_xy; y00+36*d0_xy];
    
    leftear = [x00(:) y00 leftear4_7];
    C=[ 1.2 1.2 0.27884 0.27884 0.60687 0.60687 0.2 ...
    0.2 -2.605 -2.605 -10.71177 y016b+10.71177 8.88740 -0.21260;
    0 180 90 270 90 270 0 ...
    180 15 165 90 270 0 0 ];

    if (oL==0 & oR==0)
        phantom.E=E(1:17,:);
        phantom.C=C(:,1:12);
    elseif(oL==0 & oR==1)
        phantom.E=[E;E_cavity];
        phantom.C=C;
    elseif(oL==1 & oR==0)
        phantom.E=[leftear;E(1:17,:)];
        phantom.C=C(:,1:12);
    else
        phantom.E=[leftear;E;E_cavity];
        phantom.C=C;
    end
end