#!/usr/local/bin/octave --silent

pkg load communications

addpath(fileparts(mfilename('fullpath')));

function Decoder(input_path, output_path, H_Matrix_path)
 
    % Test the *Non-binary* LDPC algorithm
    %***page_output_immediately(true)
    %***pkg load communications

    % *** Parameters ***
    % Number of input symbols, power of 2 (q : 2^q symbols)
    q = 2;
    % Probability distribution for the input symbols
    px = [0.25 0.25 0.25 0.25];
    % Probability distribution for the channel
    %pz = [0.91 0.03 0.03 0.03];
    
            %p(Y=A|X=A)
            %Y
%     pz = [0.98342075,0.00338629,0.0091428,0.00405016; %X
%           0.00446663,0.98274879,0.00370053,0.00908405; 
%           0.00860151,0.00315,0.98431995,0.00392854;
%           0.00498667,0.0100006,0.00433847,0.98067426];
      %pz=pz';

           %p(X=A|Y=A)
           %X
     pz= [0.903544 0.023989 0.044734 0.027734; %Y
           0.016697 0.917887 0.015575 0.049841;
           0.045891 0.018859 0.913311 0.021939;
           0.019837 0.043914 0.018730 0.917519;
           ] ;

      
    % Number of iterations
    it = 20;
    
    % Load the coding matrix %use eval to allow matrixName argument use
    %H = matrix()';
    [path,name,ext]=fileparts(H_Matrix_path);
    
    if(strcmp(ext,'.txt'))
        H=dlmread(H_Matrix_path);
    else
        %eval(['H = load(''' H_Matrix_path ''').H']); 
        eval(['H = load(''' H_Matrix_path ''')']); 
        H=H.H;
    end

    H=H'; %do it in prev instruction if possible

    % Sequence length
    N = length(H);

    % Build the class of parameters
    [param,cst] = c_param(q,px,pz,N,it,H);

   
    %==== read Reference file X to compare it with decoded seq OR DO it externally 
    %=====   [X] =  DO COMPARISON OUTSIDE TO SAVE TIME

    %==== Y :read consensus sequence in quaternary format ('A'->0, 'C'->1,...)
    fid_input = fopen(input_path); % open the input file
    if fid_input == -1
        error('Author:Function:OpenFile', 'consensus file not found\n');
    end
    
    fid_output = fopen(output_path,'w');
    
    while ~feof(fid_input) % feof(fid) is true when the file ends
        
        consensus_name = fgetl(fid_input); % get sequence name
    	consensus = fgetl(fid_input); % get sequence
    	
    	% if consensus sequence longer than expected -> cut excess
    	if length(consensus) > N
    		consensus = consensus(1:N);
    	end
    	
    	% if consensus sequence shorter than expected -> fill with A
        while length(consensus) < N
            consensus = [consensus ("A")];
        end
        
        quad_consensus = DNA_TO_QUAD(consensus); %convert bases into digits
        
        % *** Realizes the decoding of U ***
        U = zeros(param.M,1);
        %error(num2str(param.M));
        % Initial matrix of messages from V.N. to C.N.
        % Rem : init_m is the function to adapt to the channel
        % Messages : log P(X|Y=0)/P(X|Y=k)

        m0 = init_m(param,str2num(quad_consensus));
        % Decoding
        %U=1:N
        [Xh,Xhval] = decode(U,m0,param,cst);
        decoded_sequence = QUAD_TO_DNA(num2str(Xh(1:end)));
        %====Write Xh (estimated sequence ) on a file B.
        %Xh ->file
        fprintf(fid_output, '%s\n', consensus_name);
        fprintf(fid_output, '%s\n', decoded_sequence);
    end
	fclose(fid_input);
    fclose(fid_output);
end


function [sequence] = DNA_TO_QUAD(sequence)
    sequence = strrep(sequence,'A','0 ');
	sequence = strrep(sequence,'C','1 ');
	sequence = strrep(sequence,'G','2 ');
	sequence = strrep(sequence,'T','3 ');

end


function [sequence] = QUAD_TO_DNA(sequence)
    sequence( sequence == "0" ) = "A";
    sequence( sequence == "1" ) = "C";
    sequence( sequence == "2" ) = "G";
    sequence( sequence == "3" ) = "T";
end


arg_list = argv();
if (length(arg_list) != 3)
    error("usage : DNA_BP_LDPC_decoder.m input_path output_path H_Matrix_path\n");
end

Decoder(arg_list{1}, arg_list{2}, arg_list{3});
