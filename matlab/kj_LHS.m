function [LHS] = kj_LHS (rIn,E)

% -curl curl E


nPhi = 13;
kz = 0;
f = 13e6;

w = 2*pi*f;
c = 2.99792458e8;
eps_0 = 8.854187817e?12;

[M,N] = size(r);

for j=1+1:N-1
    
    r = rIn(j);
    
    % See the $KINETICJ/mathematica/kj_curl_curl_finite_difference file for the
    % derivation of this code chunk.
        
    curlcurl_r=i.*r.^(-2).*((-1).*i.*(nPhi.^2+kz.^2.*r.^2).*Er(j)+nPhi.*Et( ...
        j)+(1/2).*h.^(-1).*r.*((-1).*nPhi.*Et((-1)+j)+nPhi.*Et(1+j)+ ...
        kz.*r.*((-1).*Ez((-1)+j)+Ez(1+j))));
    
    
    curlcurl_t=r.^(-2).*((-1).*i.*nPhi.*Er(j)+Et(j)+(-1/2).*h.^(-2).*r.*( ...
        2.*h.^2.*i.^2.*kz.^2.*r.*Et(j)+h.*((-1).*Et((-1)+j)+Et(1+j)) ...
        +2.*r.*(Et((-1)+j)+(-2).*Et(j)+Et(1+j))+h.*i.*nPhi.*(Er((-1) ...
        +j)+(-1).*Er(1+j)+(-2).*h.*i.*kz.*Ez(j))));
    
    
    curlcurl_z=(1/2).*h.^(-2).*(h.*r.^(-2).*(i.*(kz.*r.*((-1).*r.*Er((-1)+ ...
        j)+r.*Er(1+j)+2.*h.*(Er(j)+i.*nPhi.*Et(j)))+(-2).*h.*i.* ...
        nPhi.^2.*Ez(j))+r.*(Ez((-1)+j)+(-1).*Ez(1+j)))+(-2).*(Ez(( ...
        -1)+j)+(-2).*Ez(j)+Ez(1+j)));
    
    
    % +w^2/c^2 * ( E + i/(w*eps_0)*Jp )
    
    
    
end

end