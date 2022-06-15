function [param,cst] = c_param(q,px,pz,N,it,H)
% Build the class of parameters for the LDPC algorithm
%
% q : Number of input symbols, power of 2 (q : 2^q symbols)
% px : Probability distribution for the input symbols 
% pz : Probability distribution for the channel
% N : Sequence length
% it : Number of iterations for the LDPC decoder
% H : Coding matrix
% add also param.M : number of Check nodes 
% param.r : n-th root of the unit
%
% Also return a class of constant expressions (for faster computation)
% cst.W : W matrices for decoding (matrix coefficient)
% cst.R : R matrices for decoding (check node value)
% cst.Rij : for the computation of the fft
% cst.iRij : for the inverse computation of the fft 

% Determine M
M = min(size(H,1), size(H,2));

% Build W and R
W = gen_w(q);
R = gen_r(q);

% Build Rij and iRij
icj = gf( kron([0:(2^q-1)]',ones(1,2^q)) , q ) .* gf( kron(ones(2^q,1),[0:(2^q-1)])   , q  );
Rij = exp(i*2*pi/2).^double(icj.x);
iRij = 1/(2^q)*exp(-i*2*pi/2).^double(icj.x);

% Build the structure param
param = struct('q', q,'px', px, 'pz', pz, 'N', N, 'it',it, 'H', H, 'M', M);

% Build the structure cst
cst = struct('W', W, 'R', R, 'Rij', Rij, 'iRij', iRij);

end



function [W] = gen_w(q)
	% Build the matrices W required for the Check Node processing (matrix coefficients)
	% As we are in the Galois field G(2^q), there are 2^q possible matrices of size 2^q*2^q each
	% Return W, 3D matrix, with one matrix per dimension

	% The matrix
	W = zeros(2^q,2^q,2^q);

	% For each possible value of the GF
	for j=0:(2^q-1)
		% Matrix with j(barre) * k
		Wk = kron([0:2^q-1]',ones(1,2^q)) * gf(j,q);

		% Matrix with n
		Wn = gf(kron(ones(2^q,1),[0:2^q-1]),q);

		% Substraction
		list = find((Wk-Wn)==0);

		% Build W
		A = zeros(size(W(:,:,j+1)));
		A(list) = 1;
		W(:,:,j+1) = A';
	end
end

function [R] = gen_r(q)
	% Build the matrices R required for the Check Node processing (check-node value)
	% As we are in the Galois field G(2^q), there are 2^q possible matrices of size 2^q*2^q each
	% Return R, 3D matrix, with one matrix per dimension

	% The matrix
	R = zeros(2^q,2^q,2^q);

	% For each possible value of the GF
	for j=0:(2^q-1)
		% Matrix with j(barre)-k
		Rk = gf(j,q) - gf(kron([0:2^q-1]',ones(1,2^q)),q);
	
		% Matrix with n
		Rn = gf(kron(ones(2^q,1),[0:2^q-1]),q);

		% Addition
		list = find((Rk+Rn)==0);

        	% Build R
        	A = zeros(size(R(:,:,j+1)));
        	A(list) = 1;
        	R(:,:,j+1) = A;
	end
end
