function Yb = Geom_bezier(XP,YP,XQ)
%
y=zeros(length(XQ),1);
%
p = [XP;YP]';y(:)=XQ(:);
n = size(p,1);
m = length(y);
T = zeros(n,n);
X(:,1) = p(:,1);
Y(:,1) = p(:,2);

b = XP(end);a=XP(1);
for j = 1:m
    for i = 2:n
        Y(i:n,i) = (b-y(j))/(b-a)*Y(i-1:n-1,i-1) + (y(j)-a)/(b-a)*Y(i:n,i-1);
    end
    Yb(j) = Y(n,n);
end

end
