function spd = spdiag( dvec ,k)
%
% sparse diagonal matrix
%
% spd = spdiag( dvec,k )
%
m = length(dvec);
n = abs(k)+m;
if (k == 0),
   index = 1:n;
   spd = sparse( index, index, dvec(:),n,n);
else 
  if (k > 0),
   % upper triangular
   index = 1:m;
   spd = sparse( index(:), index(:) + abs(k), dvec(:),n,n);
  else 
   % k < 0, lower triangular
   index =  1:m;
   spd = sparse( index(:)+abs(k), index(:), dvec(:),n,n);
  end;
end;


