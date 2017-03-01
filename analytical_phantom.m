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
    
    function [value] = valif(cond, val1, val2)
        if (cond)
            value = val1;
        else
            value = val2;
        end
    end
    
    % these values are summed yp
    default = true;
    tiss_1 = valif(default, 0.010, 1);   
    tiss_2 = valif(default, 0.010, 1);
    tiss_3 = valif(default, 0.0025,0);
    tiss_4 = valif(default, -0.0025, 0);
    tiss_5 = valif(default, 1.800, 3);    % white phantom part
    tiss_6 = valif(default, -0.750, -2);  % internal phantom part
    tiss_7 = valif(default, -1.050, -1);   % black oval on internal
    tiss_8 = valif(default, 0.750, 2);
    tiss_9 = valif(default, 0.750, 2);
    tiss_10 = valif(default, 0.750, 2);
    tiss_11 = valif(default, 0.750, 2);
    tiss_12 = valif(default, -0.005, 0);
    tiss_13 = valif(default, 0.005, 0);
    tiss_14 = valif(default, 0.750, 2);
    tiss_15 = valif(default, 1.800, 3);
    tiss_16a = valif(default, 0.750, 2);
    tiss_16b = valif(default, 0.750, 2);
    tiss_R_ear = valif(default, 0.750, 2);
    tiss_left_ear = valif(default, 0.750, 2);
    tiss_right_ear = valif(default, -1.800, -3);
    
        E = [-4.7 4.3 1.79989 1.79989 0 tiss_1 0; %1
          4.7 4.3 1.79989 1.79989 0 tiss_2 0; %2
        -1.08 -9 0.4 0.4 0 tiss_3 0; %3
         1.08 -9 0.4 0.4 0 tiss_4 0; %4
         0 0 9.6 12 0 tiss_5 0; %5
         0 8.4 1.8 3.0 0 tiss_7 0; %7
         1.9 5.4 0.41633 1.17425 -31.07698 tiss_8 0; %8
        -1.9 5.4 0.41633 1.17425 31.07698 tiss_9 0; %9
        -4.3 6.8 1.8 0.24 -30 tiss_10 0; %10
         4.3 6.8 1.8 0.24 30 tiss_11 0; %11
         0 -3.6 1.8 3.6 0 tiss_12 0; %12
         6.39395 -6.39395 1.2 0.42 58.1 tiss_13 0; %13
         0 3.6 2 2 0 tiss_14 4; %14
         0 9.6 1.8 3.0 0 tiss_15 4; %15
         0 0 9.0 11.4 0 tiss_16a 3; %16a
         0 y016b a16b b16b 0 tiss_16b 1; %16b
         0 0 9.0 11.4 0 tiss_6 oR; %6
         9.1 0 4.2 1.8 0 tiss_R_ear 1];%R_ear

    %generate the air cavities in the right ear
    cavity1 = transpose(8.8:-0.4:5.6);
    cavity2 = zeros(9,1);
    cavity3_7 = ones(53,1)*[0.15 0.15 0 tiss_right_ear 0];
    for j = 1:3 
        kj = 8-2*floor(j/3); 
        dj = 0.2*mod(j,2);
        cavity1 = [cavity1; cavity1(1:kj)-dj; cavity1(1:kj)-dj];
        cavity2 = [cavity2; j*sha*ones(kj,1); -j*sha*ones(kj,1)];
    end
    E_cavity = [cavity1 cavity2 cavity3_7];

    %generate the left ear (resolution pattern)
    leScale = 3;
    x0 = -7.0;
    y0 = -1.0;
    d0_xy = 0.04 * leScale;   % spacing

    d_xy = [0.0357, 0.0312, 0.0278, 0.0250];
    x00 = zeros(0,0);
    y00 = zeros(0,0);
    ab = (0.5 * leScale)*ones(5,1)*d_xy;
    ab = ab(:)*ones(1,4);
    leftear4_7 = [ab(:) ab(:) ones(80,1)*[0 tiss_left_ear 0]];
    for i = 1:4
        y00 = [y00; transpose(y0+(0:4)*(2 * leScale)*d_xy(i))];
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