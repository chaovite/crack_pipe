function Y = coswindow(y,t,T)

    if length(T)~=4,disp('length of T must be 4'),return,end
    if length(t)~=length(y),disp('length of t and y must match'),return,end

    I1 = find(t<T(1));
    I2 = find(T(1)<=t&t<=T(2));
    I3 = find(T(3)<=t&t<=T(4));
    I4 = find(T(4)<t);

    w2 = 0.5-0.5*cos(pi*(t(I2)-T(1))/(T(2)-T(1)));
    w3 = 0.5+0.5*cos(pi*(t(I3)-T(3))/(T(4)-T(3)));
    
    Y = y;
    Y(I1) = 0;
    Y(I2) = w2.*Y(I2);
    Y(I3) = w3.*Y(I3);
    Y(I4) = 0;
