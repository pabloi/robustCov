function eh=drawEllipse2D(M)
%Draws a 2D ellipse given by %x'*M*X=1 on current axes. M bust be a
%postivie semi-definite matrix for it to work properly,
if nargin<1 || isempty(M)
    M=eye(2);
end
th=0:.1:2*pi;
x=sin(th);
y=cos(th);
a=sqrt(1./sum([x;y].*(M*[x;y])));
eh=plot(a.*x,a.*y);
end