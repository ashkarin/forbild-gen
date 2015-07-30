% Transforms all attenuation coefficients within 
% the phantom into physical values based on
%{

Yu, Z., Noo, F., Dennerlein, F., Wunderlich, A., Lauritsch, G., & 
Hornegger, J. (2012). Simulation tools for two-dimensional experiments in 
x-ray computed tomography using the FORBILD head phantom. 
Physics in medicine and biology, 57(13), N237.

%}

% two materials, water and bone
% muW - attenuation coefficient of water
% muB - attenuation coefficient of bone
function [ physphantom ] = physical_value(phantom,muW,muB)
    physphantom.E=phantom.E;
    physphantom.C=phantom.C;
    nrows=size(phantom.E,1);
    shift=0;
    if (nrows > 71) 
        shift=80;
    end

    if (nrows >= 97)
        physphantom.E(1:80,6)=muB-1.05*muW;
    end

    if (nrows == 71 || nrows == 151)
        physphantom.E(18+shift,6)=muB-1.05*muW;
        physphantom.E((19:71)+shift,6)= -muB;
    end

    physphantom.E(5+shift,6)=muB;
    physphantom.E(17+shift,6)=1.05*muW-muB;
    physphantom.E(6+shift,6)=-1.05*muW;
    physphantom.E(14+shift,6)=muB;
    physphantom.E([7 8 9 10 13 15 16]+shift,6)=muB-1.05*muW;
    j=[1 2 3 4 11 12];
    physphantom.E(j+shift,6)=muW*phantom.E(j+shift,6);
end