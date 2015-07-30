function [ sino ] = line_integrals(scoord,theta,phantom)
    sinth = sin(transpose(theta(:)));
    costh = cos(transpose(theta(:)));
    meps = 1e-10;
    nclipinfo = 0;
    mask = zeros(size(theta));

    for k=1:size(phantom.C,2)
        psi = phantom.C(2,k)*pi/180;
        tmp = -sinth*cos(psi)+costh*sin(psi);
        kk = find(abs(tmp)<meps); mask(kk) = meps;
    end

    theta= theta+mask;
    sino = zeros(1,length(scoord(:)));
    sinth = sin(transpose(theta(:)));
    costh = cos(transpose(theta(:)));
    sx = transpose(scoord(:)).*costh;
    sy = transpose(scoord(:)).*sinth;

    for k=1:length(phantom.E(:,1))
        x0 = phantom.E(k,1);
        y0 = phantom.E(k,2);
        a = phantom.E(k,3);
        b = phantom.E(k,4);
        phi = phantom.E(k,5)*pi/180;
        f = phantom.E(k,6);
        nclip = phantom.E(k,7);
        s0 = [sx-x0;sy-y0];
        DQ = [cos(phi)/a sin(phi)/a; -sin(phi)/b cos(phi)/b];
    DQthp = DQ*[-sinth;costh];
    DQs0 = DQ*s0;
    A = sum(DQthp.^2);
    B = 2*sum(DQthp.*DQs0);
    C = sum(DQs0.^2)-1;
    equation = B.^2-4*A.*C;
    i = find(equation>0);
    tp = 0.5*(-B(i)+sqrt(equation(i)))./A(i);
    tq = 0.5*(-B(i)-sqrt(equation(i)))./A(i);

    if (nclip>0)
        for j = 1:nclip
            nclipinfo = nclipinfo+1;
            d = phantom.C(1,nclipinfo);
            psi= phantom.C(2,nclipinfo)*pi/180;
            xp = sx(i)-tp.*sinth(i);
            yp = sy(i)+tp.*costh(i);
            xq = sx(i)-tq.*sinth(i);
            yq = sy(i)+tq.*costh(i);
            tz = d-cos(psi)*s0(1,i)-sin(psi)*s0(2,i);
            tz = tz./(-sinth(i)*cos(psi)+costh(i)*sin(psi));
            equation2 = ((xp-x0)*cos(psi)+(yp-y0)*sin(psi));
            equation3 = ((xq-x0)*cos(psi)+(yq-y0)*sin(psi));
            m1 = find(equation3>=d); tq(m1) = tz(m1);
            m2 = find(equation2>=d); tp(m2) = tz(m2);
        end
    end
    sinok = f*abs(tp-tq);
    sino(i) = sino(i)+sinok;
    end
    sino = reshape(sino,size(theta));
end 