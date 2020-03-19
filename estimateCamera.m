function [R,t] = estimateCamera2(xy,XYZ)
    %M = K [R t]
	u = xy(:,1);
	v = xy(:,2);
	X = XYZ(:,1);
	Y = XYZ(:,2);
	Z = XYZ(:,3);
	A=[];

	for i = 1:1:size(xy,1)
	    r1=[X(i) Y(i) Z(i) 1 0 0 0 0 -u(i)*X(i) -u(i)*Y(i) -u(i)*Z(i) -u(i)];
	    r2=[0 0 0 0 X(i) Y(i) Z(i) 1 -v(i)*X(i) -v(i)*Y(i) -v(i)*Z(i) -v(i)];
	    A = [A; r1; r2];
	end

	[~,~,V] = svd(A);
	M = reshape(V(:,end),4,3)';
    R = M(:,1:end-1);
    R = Gram_Schimidt(R); 
    t = M(:,end)
end
function [Q] = Gram_Schimidt(A)
    [m,n] = size(A); Q = zeros(m,n); R = zeros(n,n);
    for i=1:n
        v = A(:,i);
        for j=1:i-1
            R(j,i)=Q(:,j)'*A(:,i);
            v=v-R(j,i)*Q(:,j);
        end
        R(i,i)=norm(v);
        Q(:,i)=v/R(i,i);
    end
end