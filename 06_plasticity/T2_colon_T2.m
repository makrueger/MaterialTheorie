function c = T2_colon_T2(A,B)
%
% double contraction 
% c = A(i,j)*B(i,j)
%
[n,m] = size(A);

c = 0;
for i=1:n
    for j=1:n
        c = c + A(i,j)*B(i,j);
    end;
end;
