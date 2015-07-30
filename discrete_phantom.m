%{

Yu, Z., Noo, F., Dennerlein, F., Wunderlich, A., Lauritsch, G., & 
Hornegger, J. (2012). Simulation tools for two-dimensional experiments in 
x-ray computed tomography using the FORBILD head phantom. 
Physics in medicine and biology, 57(13), N237.

%}

function [ image ] = discrete_phantom(xcoord,ycoord,phantom)
    image = zeros(size(xcoord));
    nclipinfo = 0;
    for k = 1:length(phantom.E(:,1))
        Vx0 = [transpose(xcoord(:))-phantom.E(k,1);
               transpose(ycoord(:))-phantom.E(k,2)];
          D = [1/phantom.E(k,3) 0; 0 1/phantom.E(k,4)];
          phi = phantom.E(k,5)*pi/180;
          Q = [cos(phi) sin(phi); -sin(phi) cos(phi)];
          f = phantom.E(k,6);
          nclip = phantom.E(k,7);
          equation1 = sum((D*Q*Vx0).^2);
          i = find(equation1<=1.0);

          if (nclip > 0)
              for j = 1:nclip
                  nclipinfo = nclipinfo+1;
                  d = phantom.C(1,nclipinfo);
                  psi = phantom.C(2,nclipinfo)*pi/180;
                  equation2 = ([cos(psi) sin(psi)]*Vx0);
                  i = i(find(equation2(i)<d));
              end
          end
          image(i) = image(i)+f;
    end
end