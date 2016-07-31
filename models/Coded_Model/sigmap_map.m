% -------------------------------------------------------------------------
% Author: Rodrigo Calderon Rico
%         Phd. Candidate,
%         University of Guadalajara,
%         Guadalajara, Mexico.
% -------------------------------------------------------------------------

function [symbols,ldM,maxout] = sigmap_map(bit,mod_typ,MAP)

% =========================================================================
% Funktion: Signal-Mapping
%           This function has been created for bit-encoding. The variable BIT 
%           is a vector or matrix with incomming bit. Depending on the 
%           modulation-type mod_typ a complex symbol-matrix (or vector) named 
%           SYMBOLS will be created. Additional the number of bit per symbol 
%           ldM will be given back, too.
%
% call:   [symbols,ldM,maxout] = sigmap2(bit,mod_typ);
%
% Input:  bit     : Bitvector or -matrix with size (high x dummy)
%				mod_typ : Modulation-Technique
%		            		1		BPSK        Mapping: 0->+1 und 1->-1
%				            2		QPSK
%				            3		8-PSK
%				            4		16-QAM
%				            5		16-PSK
%				            6		64-QAM
%				            8		256-QAM
%           MAP     : Bit-Symbol-Mapping
%                       0     Gray-Mapping mit 0->+1 und 1->-1 (default)
%                       1     Gray-Mapping mit 0->-1 und 1->+1
%
% Output:  symbols : complex symbolvector or matrix with size (high/M x dummy)
%				ldM     : number of bit per symbol
%				maxout  : maximal output (REAL or IMAG);
%
%
% =========================================================================


% Anzahl der Eingabeparameter ---------------------------------------------
if nargin<3
   MAP = 0;
end

% Variablendimension festellen --------------------------------------------
[rows,cols] = size(bit);	   	% get matrix size
if (cols == 1)				    % make row vector (only if line vector)
   bit = bit.';
   [rows,cols] = size(bit);
end;


switch mod_typ 
   
case 1 	 %  ###########  BPSK  #############

   ldM   	= 1; 
   maxout   = 1;
   if ~MAP
      symbols 	= (1 - 2 * bit);
   elseif MAP==1
      symbols 	= (2 * bit - 1);
   end
   
case 2   %  ###########  QPSK  #############
   
   ldM 	   = 2;
   new_cols = floor (cols / ldM);
   w_qpsk 	= 1/sqrt(2);
   maxout   = w_qpsk;
   if ~MAP
      symbols(:,1:new_cols) = w_qpsk * ...
         (1-2*bit(:,ldM*[1:new_cols]-1) + j*(1 -2*bit(:,ldM*[1:new_cols])));
      
   elseif MAP==1
      symbols(:,1:new_cols) = w_qpsk * ...
         (2*bit(:,ldM*[1:new_cols]-1)-1 +i*(2*bit(:,ldM*[1:new_cols])-1));
   end
   
case 3   %  ###########  8-PSK  ############
   
   ldM 	   = 3;
   new_cols = floor (cols / ldM);
   a 	      = sin(pi/8);
   b 	      = cos(pi/8);
   maxout   = b;
   if ~MAP
      symbols(:,1:new_cols) = ((b-(b-a)*bit(:,3*[1:new_cols])).*(1-2*bit(:,3*[1:new_cols]-1)) ...
         +j*(a+(b-a)*bit(:,3*[1:new_cols])).*(1-2*bit(:,3*[1:new_cols]-2)));
   elseif MAP==1
      symbols(:,1:new_cols) = ((b-(b-a)*bit(:,3*[1:new_cols])).*(2*bit(:,3*[1:new_cols]-1)-1) ...
         +j*(a+(b-a)*bit(:,3*[1:new_cols])).*(2*bit(:,3*[1:new_cols]-2)-1));
   end
   
case 4   %  ###########  16-QAM  ############
   
   ldM 	   = 4;
   new_cols = floor (cols / ldM);
   
   if ~MAP
      bit = ~bit;
   elseif MAP==1
      bit = bit;
   end

   bit3    = 3-2*bit(:,ldM*[1:new_cols]  );
   bit2    = 2*  bit(:,ldM*[1:new_cols]-1)-1;
   bit1    = 3-2*bit(:,ldM*[1:new_cols]-2);
   bit0    = 2*  bit(:,ldM*[1:new_cols]-3)-1;
   norm    = 1/sqrt(10);
   maxout  = 3*norm;
   symbols = norm * ((bit0.*bit1) + j*(bit2.*bit3));
   
case 5   %  ###########  16-PSK  ###########
   
   ldM 	   = 4;
   new_cols = floor (cols / ldM);
   
   dc = cos(pi/16);
   cc = cos(3*pi/16);
   bc = cos(5*pi/16);
   ac = cos(7*pi/16);
   
   ds = sin(pi/16);
   cs = sin(3*pi/16);
   bs = sin(5*pi/16);
   as = sin(7*pi/16);
   
   if ~MAP
      bit1 = 1-2*bit(:,ldM*[1:new_cols]-3);
      bit0 = 1-2*bit(:,ldM*[1:new_cols]-2);
      bit2 =     bit(:,ldM*[1:new_cols]-1);
      bit3 =     bit(:,ldM*[1:new_cols]  );
      symbols = bit0 .* ((ac-dc) * bit2 + (cc-dc) *bit3 + (dc+bc-ac-cc)*bit2.*bit3 + dc) ...
         + j * (bit1 .* ((as-ds) * bit2 + (cs-ds) *bit3 + (ds+bs-as-cs)*bit2.*bit3 + ds));
      
   elseif MAP==1
      bit0 = 2*bit(:,ldM*[1:new_cols]-3)-1;
      bit1 = 2*bit(:,ldM*[1:new_cols]-2)-1;
      bit2 =   bit(:,ldM*[1:new_cols]-1);
      bit3 =   bit(:,ldM*[1:new_cols]  );
      symbols = bit0 .* ((ac-dc) * bit2 + (cc-dc) *bit3 + (dc+bc-ac-cc)*bit2.*bit3 + dc) ...
         + j * (bit1 .* ((as-ds) * bit2 + (cs-ds) *bit3 + (ds+bs-as-cs)*bit2.*bit3 + ds));
   end
   maxout = dc;
   
case 6   %  ###########  64-QAM  ############
   
   ldM  = 6;
   if ~MAP
      bit = ~bit;
   elseif MAP==1
      bit = bit;
   end
   new_cols = floor (cols / ldM);
   
   bit5 = 3-2*bit(:,ldM*[1:new_cols]  );
   bit4 =   2*bit(:,ldM*[1:new_cols]-1)-1;
   bit3 =   2*bit(:,ldM*[1:new_cols]-2)-1;
   bit2 = 3-2*bit(:,ldM*[1:new_cols]-3);
   bit1 =   2*bit(:,ldM*[1:new_cols]-4)-1;
   bit0 =   2*bit(:,ldM*[1:new_cols]-5)-1;
   norm = 1/sqrt(42);
   offs = 4;
   symbols = norm * ((bit0.*(offs-bit1.*bit2)) + j*(bit3.*(offs-bit4.*bit5)));
   maxout  = 7*norm;
   
case 8   %  ########### 256-QAM  ############
   
   ldM  = 8;   
   if ~MAP
      bit = ~bit;
   elseif MAP==1
      bit = bit;
   end

   new_cols = floor (cols / ldM);
   
   bit7 = 3-2*bit(:,ldM*[1:new_cols]  );
   bit6 =   2*bit(:,ldM*[1:new_cols]-1)-1;
   bit5 =   2*bit(:,ldM*[1:new_cols]-2)-1;
   bit4 =   2*bit(:,ldM*[1:new_cols]-3)-1;
   bit3 = 3-2*bit(:,ldM*[1:new_cols]-4);
   bit2 =   2*bit(:,ldM*[1:new_cols]-5)-1;
   bit1 =   2*bit(:,ldM*[1:new_cols]-6)-1;
   bit0 =   2*bit(:,ldM*[1:new_cols]-7)-1;
   norm = 1/sqrt(170);
   offs = 4;
   offs2 = 8;
   symbols = norm * ( bit0.*(offs2-(bit1.*(offs-bit2.*bit3))) + i* ( bit4.*(offs2-(bit5.*(offs-bit6.*bit7))) ) );
   maxout = 15*norm;
   
end;


% EOF


