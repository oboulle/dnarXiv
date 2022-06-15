function [m0] = init_m(param,Y)
% Build the initial matrix of messages from V.N. to C.N.

% *** Symbol part ***
m0s = kron([0, log(param.px(1)./param.px(2:end))],ones(param.N,1));

% *** Channel part ***
% For P(Y|X=0) 

mx0y = log(kron(param.pz(Y+1,1),ones(1,2^param.q)));

% Indices j-k
%jmk = gf(kron(Y,ones(1,2^param.q)),param.q) - gf(kron(ones(param.N,1),[0:2^param.q-1]),param.q); 
% For P(Y|X=k)
mxky = zeros(length(Y),2^param.q);
for k=1:2^param.q
	mxky(:,k) = log(param.pz(Y+1,k));
end

% *** Initial messages ***
% Rem : One vector of messages per line
m0 = mx0y - mxky + m0s;

end
