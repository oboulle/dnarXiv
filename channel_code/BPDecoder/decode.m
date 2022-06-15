function [Xh,Xhval,sx] = decode(U,m0,param,cst)
% Realizes the decoding of U from the initial messages m0 and
% param, cst : see c_param.m
% Nb-binary case, Circular Gaussian noise

% N and M
N = size(m0,1);
M = length(U);

% *** Initial matrix of messages from m0 ***
mvc = zeros(N,M,2^param.q);
% Non-zero coefficients in the coding matrix = link between a V.N. and a C.N.
[line,col] = find(param.H); % 
% Places the initial messages in the initial matrix
for j=1:M
        mvc(line(find(col==j)),j,:) = m0(line(find(col==j)),:);
end

% *** Initialization of the decoding algorithm ***
%Initial check of the syndrom
sx = 11;
% Initial number of iterations
j = 0;
% Contains the estimate of X
Xh = zeros(N,1);


% *** Iterations of the algorithm *** 
while((j<param.it) && (sx>0))

	% Messages from C.N. to V.N.
	mcv = check_node(param,cst,mvc,U);
	
	% Messages from V.N. to C.N.
	mvc = variable_node(param,mcv,m0);
	
	% Estimated version of X
	Xhval = reshape(sum(mcv,2),N,2^param.q) + m0;
	[v,p] = min(Xhval');
	Xh = (p-1)';

	% Syndrom
	% sx = length(find( U - rem(param.G'*Xh,2^param.q) ))
	sx = gf(U,param.q) - gf(full(param.H)',param.q) * gf(Xh,param.q);
	sx = length(find(sx.x));

	% Iterations
	j = j+1;
end


end


function [mcv] = check_node(param,cst,mvc,U)
% Response messages of the check nodes from the messages of the variable nodes

% N and M
N = size(mvc,1);
M = length(U);

% Initialize the response matrix
mcv = zeros(size(mvc));

% For each check node...
for j=1:M
	% Find the non-zero components
	pos = find(param.H(:,j));
	% Associated coefficients
	coeff = param.H(pos,j);

	% For each non-zero component : response message
	for k=1:length(pos)
		% Right vector of components (current extracted)
		posk = pos;
		posk(k) = 0;
		coeffk = coeff(find(posk));
		posk = posk(find(posk));

		% Compute the modified coefficients H(barre) and z(barre)
		Hb = -gf(coeffk,param.q)./gf(coeff(k),param.q);
		zb = gf(U(j),param.q)/gf(coeff(k),param.q);

		% Extract the concerned messages
		mc = reshape(mvc(posk,j,:),length(posk),2^param.q);
		Ma = Hb.x;
		% Multiply the messages by the associated W matrices
		for l=1:length(posk)
			mc(l,:) = (cst.W(:,:,Ma(l)+1) * mc(l,:)')'; 
		end

		% Compute the product of the Fourier transforms then the inverse Fourier transform
		mcff = iftmgf( prod(ftmgf(mc,cst),1),cst );
		% Apply the R matrix (in p-space)
		pcf = (cst.R(:,:,zb.x+1) * (exp(-mcff)/sum(exp(-mcff)))')';
		% Come back to m-space
		mcf = [0, log(pcf(1)./pcf(2:end))];

		% Position the resulting message in the message matrix
		mcv(pos(k),j,:) = real(reshape(mcf,1,1,2^param.q));
	end
end

end

function [mvc,Xh] = variable_node(param,mcv,m0)
% Response messages of the variable nodes from the messages of the check nodes

% N and M
N = size(m0,1);
M = size(m0,2);

% Initialize the response matrix
mvc = zeros(size(mcv));

% For each variable node...
for k=1:N
	% Extract the positions
	pos = find(param.H(k,:));

	% For each active component
	for j=1:length(pos)
		% Extract the right components (exclude the current one)
		posj = pos;
		posj(j) = 0;
		posj = posj(find(posj));

		% Compute the messages
		mvc(k,pos(j),:) = reshape(m0(k,:),1,1,2^param.q) + sum(mcv(k,posj,:),2);
	end

end

end

function [F] = ftmgf(m,cst)
% Compute the Fourier transforms of a matrix of messages m (one per line) in LLR domain
% The Fourier transform is the one associated to GF(2^param.q)

% Constant term
c = kron(sum(exp(-m),2),ones(1,size(m,2)));

% Fourier Transform
F = (cst.Rij * exp(-m)')' ./c;

end

function [m] = iftmgf(F,cst)
% Compute the inverse Fourier Transform (associated to GF(2^param.q) of F (matrix, one vector per line)
% Return m in the LLR domain

% Compute the probabilities in vectorized form
p = (cst.iRij * F')';

% Compute the messages
m = zeros(size(p));
m(:,2:end) = kron(log(p(:,1)),ones(1,size(F,2)-1)) - log(p(:,2:end));

end
