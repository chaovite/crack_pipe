function xout = taper(x, ratio)
% taper the input x using cosine function
% use built in function tukeywin.
% see also tukeywin
if size(x, 1) > 1
    xout = x.*tukeywin(length(x),ratio);
else
    xout = x.*(tukeywin(length(x),ratio))';
end

