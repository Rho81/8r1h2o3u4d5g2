function [softbit,hardbit] = sigdemap(symbols,mod_typ,MAP,Zeilen,Spalten)

% ==============================================================================
% Funktion: Signal-De-Mapping
%           This function has been created for bit-decoding. The variable 
%           symbols is a vector or matrix of outcomming symbols. Depending on 
%           the modulation-type mod_typ a binary bit-matrix (or vector) named 
%           BIT will be created. Additional a soft-information will be generated.
%
% Aufruf:   [softbit,hardbit] = sigdemap2(symbols,mod_typ,Zeilen,Spalten);
%
% Eingabe:  symbols : complex symbolvector or matrix with size (high/M x dummy)
%				mod_typ : Modulation-Technique
%		            		1		BPSK        Mapping: 0->+1 und 1->-1
%				            2		QPSK
%				            3		8-PSK
%				            4		16-QAM
%				            5		16-PSK
%				            6		64-QAM
%				            8		256-QAM
% Ausgabe:  softbit : Softbitvector or -matrix with size (high x dummy)
%           hardbit : Hardbitvecotr of -matrix with size (high x dummy)
%                     Mapping 0->+1 und 1->-1
%
% Anmerkung:   
%
% M-Files:  

if nargin <3
   MAP = 0;
end

if nargin<4
   RESHAPE = 0;
else
   RESHAPE = 1;

   if nargin <5
      Spalten = [];
   end
end


[rows,cols] = size(symbols);	   	% get matrix size


switch mod_typ 
   
case 1 	 %  ###########  BPSK  #############
   
   ldM   	= 1;
   a1 	   = 1;
   softbit  = real(symbols);
   %	softbit = sqrt(2)*a1*softbit;
   %	bit	= (sign(softbit)+1) /2;
   
   
case 2   %  ###########  QPSK  #############
   
   ldM 	 = 2;
   new_cols = (cols * ldM);
   iw_qpsk  = sqrt(2);
   a1	      =1/iw_qpsk;
   softbit(:,1:ldM:(new_cols-1))  = iw_qpsk*real(symbols);
   softbit(:,2:ldM:(new_cols  ))  = iw_qpsk*imag(symbols);
   
   softbit = softbit / iw_qpsk;
   %	bit	= (sign(softbit)+1) /2;
   
case 3   %  ###########  8-PSK  ############
   
   ldM 	   = 3;
%   a1       = sin(pi/8);
   new_cols = floor (cols * ldM);
   softbit(:,1:ldM:new_cols-2) = imag(symbols);
   softbit(:,2:ldM:new_cols-1) = real(symbols);
   
   if ~MAP
      softbit(:,3:ldM:new_cols) = (-abs(imag(symbols)) + abs(real(symbols)))/sqrt(2);
   elseif MAP==1
      softbit(:,3:ldM:new_cols) = (abs(imag(symbols)) - abs(real(symbols)))/sqrt(2);
   end
   
   
   
case 4   %  ###########  16-QAM  ############
   
   ldM 	 = 4;
   new_cols = floor (cols * ldM);
   inorm    = sqrt(10);
   symbols	 = symbols * inorm;
   
   softbit(:,1:ldM:new_cols-3) = real(symbols);
   softbit(:,2:ldM:new_cols-2) = -(abs(real(symbols))-2);
   softbit(:,3:ldM:new_cols-1) = imag(symbols);
   softbit(:,4:ldM:new_cols  ) = -(abs(imag(symbols))-2);
   
   softbit = softbit / inorm;
   %	bit = round((sign(softbit)+1) /2);
   
   
case 5   %  ###########  16-PSK  ############
   
   ldM 	 = 4;
   new_cols = floor (cols * ldM);
   a1 = sin(pi/16);
   a = sin(pi/8);
   b = cos(pi/8);
   
   if ~MAP
      softbit(:,1:ldM:new_cols-3) = imag(symbols);
      softbit(:,2:ldM:new_cols-2) = real(symbols);
      softbit(:,3:ldM:new_cols-1) = (-abs(imag(symbols)) + abs(real(symbols)))/sqrt(2);
      softbit(:,4:ldM:new_cols  ) = ...
         - ((1+sign(softbit(:,3:ldM:new_cols-1)))/2) .* (b*abs(imag(symbols)) - a*abs(real(symbols))) ...  % u3=+1
         - ((1-sign(softbit(:,3:ldM:new_cols-1)))/2) .* (b*abs(real(symbols)) - a*abs(imag(symbols))) ;    % u3=-1
   elseif MAP==1
      softbit(:,1:ldM:new_cols-3) = real(symbols);
      softbit(:,2:ldM:new_cols-2) = imag(symbols);
      softbit(:,3:ldM:new_cols-1) = (abs(imag(symbols)) - abs(real(symbols)))/sqrt(2);
      softbit(:,4:ldM:new_cols  ) = ((1-sign(softbit(:,3:ldM:new_cols-1)))/2) .* (b*abs(imag(symbols)) - a*abs(real(symbols))) + ((1+sign(softbit(:,3:ldM:new_cols-1)))/2) .* (b*abs(real(symbols)) - a*abs(imag(symbols))) ;
   end
   
case 6   %  ###########  64-QAM  ############
   
   ldM 	 = 6;
   new_cols = floor (cols * ldM);
   inorm    = sqrt(42);
   a1	 = 1/inorm;
   symbols	 = symbols * inorm;
   
   softbit(:,1:ldM:new_cols-5) = real(symbols);
   softbit(:,2:ldM:new_cols-4) = -(abs(real(symbols))-4);
   softbit(:,3:ldM:new_cols-3) = -(abs(abs(real(symbols))-4)-2);
   softbit(:,4:ldM:new_cols-2) = imag(symbols);
   softbit(:,5:ldM:new_cols-1) = -(abs(imag(symbols))-4);
   softbit(:,6:ldM:new_cols  ) = -(abs(abs(imag(symbols))-4)-2);
   
   softbit = a1*softbit;
   %	bit = round((sign(softbit)+1) /2);
   
   
case 8   %  ########### 256-QAM  ############
   
   ldM 	 = 8;
   new_cols = floor (cols * ldM);
   inorm    = sqrt(170);
   a1	 = 1/inorm;
   symbols	 = symbols * inorm;
   
   softbit(:,1:ldM:new_cols-7) = real(symbols);
   softbit(:,2:ldM:new_cols-6) = -(abs(real(symbols))-8);
   softbit(:,3:ldM:new_cols-5) = -(abs(abs(real(symbols))-8)-4);
   softbit(:,4:ldM:new_cols-4) = -(abs(abs(abs(real(symbols))-8)-4)-2);
   softbit(:,5:ldM:new_cols-3) = imag(symbols);
   softbit(:,6:ldM:new_cols-2) = -(abs(imag(symbols))-8);
   softbit(:,7:ldM:new_cols-1) = -(abs(abs(imag(symbols))-8)-4);
   softbit(:,8:ldM:new_cols  ) = -(abs(abs(abs(imag(symbols))-8)-4)-2);
   
   softbit = a1*softbit;
   %	bit = round((sign(softbit)+1) /2);
   
end;

if RESHAPE
   softbit = reshape(softbit,Zeilen,Spalten);
else
   softbit = softbit;
end

if ~MAP
   hardbit = (1-sign(softbit))/2;
elseif MAP==1
   hardbit = (1+sign(softbit))/2;
end
% EOF

