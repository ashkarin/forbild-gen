%{

Yu, Z., Noo, F., Dennerlein, F., Wunderlich, A., Lauritsch, G., & 
Hornegger, J. (2012). Simulation tools for two-dimensional experiments in 
x-ray computed tomography using the FORBILD head phantom. 
Physics in medicine and biology, 57(13), N237.

%}
function [ phantom ] = forbild_gen(W, H, oL, oR, physical, muW, muB)
    if nargin < 3
        oL = 1;
    end
    
    if nargin < 4
        oR = 1;
    end
    
    if nargin < 5
        physical = 0;
    end
    
    if nargin < 6
        muW = 1;
    end

    if nargin < 7
        muB = 1;
    end
    
    x_step = 25 / W; % 25cm should be sampled on W pixels
    y_step = 25 / H; % 25cm should be sampled on H pixels

    
    x = ((0:W-1) - W / 2) * x_step;
    y = ((0:H-1) - H / 2) * y_step;

    xcoord = ones(W, 1) * x;
    ycoord = transpose(y) * ones(1, H);

    meta_phantom = analytical_phantom(oL, oR);
    if physical
        meta_phantom = physical_value(meta_phantom, muW, muB);
    end

    phantom = discrete_phantom(xcoord, ycoord, meta_phantom);
end

