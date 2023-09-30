function [M,R,V,xs,M_Diag,V_Diag,d_Diag,D_max] = three_moment(L,I,E,w,P,x)
% [M,R,V,xs,M_Diag,V_Diag,d_Diag,D_max] = three_moment(L,I,E,w,P,x)
% solve the three-moment equations for a continuous beam of N spans.
%
% INPUT: 
%   L is a vector of length N containing the lengths of each span.
%   I is a vector of length N containing the moments of inertia of each span.
%   E is a scalar constant for the modulus of elasticity.  
%   w is a vector of length N containing the uniform loads on each span.
%   P is a vector of point loads magnitudes
%   x is a vector of point load locations, from the start of the first span
%
% OUTPUT:
%   M      is a vector of length N+1 containing the moments at each support
%   R      is a vector of length N+1 containing the reactions at each support
%   V      is a matrix of size 2xN+1 containing end-shears of each span
%   xs     is a vector of the x-axis for the shear, moment, and displ plots
%   M_Diag is a vector of the moment diagram
%   V_Diag is a vector of the shear diagram
%   d_Diag is a vector of the displacement diagram
%   D_max  is a vector of the max absolute displacement of each span
%
% This program assumes that none of the supports are moment resisting, 
% that there are no hinges in the beam, that all the spans are made of 
% the same material, and that each spans is prismatic.  

% H.P. Gavin, Civil and Environmental Engineering, Duke University, 3/24/09

 N = length(L);				% number  of spans

 nP = length(x);			% number of concentrated point loads
 span = zeros(1,nP);			% spans containing the point loads
 sumL = cumsum(L);
 xL   = zeros(1,nP);			% distance from point load to left rctn
 xR   = zeros(1,nP);			% distance from point load to right rctn
 for i=1:length(x)
	span(i) = min( find( x(i) < sumL ) );
 end

 for k=1:nP				% loop over all concentrated point loads
	if span(k) == 1			% the point load is in the first span
		xL(k) = x(k);
	else		
		xL(k) = x(k) - sumL(span(k)-1);
	end
	xR(k) = sumL(span(k)) - x(k);
 end

 F = zeros(N+1,N+1);			% initialize the flexibility matrix 

 for j=2:N				% create the flexibility matrix (8)-(10)

	F(j,j-1) = L(j-1) / I(j-1);

	F(j,j)   = 2 * (  L(j-1) / I(j-1)  +  L(j) / I(j)  );

	F(j,j+1) = L(j) / I(j);
 end

 F(1,1)     = 1.0;
 F(N+1,N+1) = 1.0;

 d = zeros(N+1,1);
 for j=2:N				% create the right-hand-side vector 

	l = j-1;			% j-1  is the number of the left  span 
	r = j;  			% j    is the number of the right span

	d(j) = -w(l)*L(l)^3 / (4*I(l)) - w(r)*L(r)^3 / (4*I(r));

	for k=1:nP			% loop over all concentrated point loads
		if span(k) == l		% the point load is in the left span
 			d(j) = d(j) - P(k)*xL(k)/(L(l)*I(l))*(L(l)^2-xL(k)^2);
		end
		if span(k) == r		% the point load is in the right span
 			d(j) = d(j) - P(k)*xR(k)/(L(r)*I(r))*(L(r)^2-xR(k)^2);
		end
	end
 end

 M = ( inv(F) * d )';			% compute the internal moments  (7)

 R = zeros(1,N+1);			% build the vector of reaction forces
 for j=1:N+1				% j is the reaction number
	
	l = j-1;			% j-1  is the number of the left  span 
	r = j;  			% j    is the number of the right span

	if j == 1
		R(j) = w(r)*L(r)/2 - M(j)/L(r) + M(j+1)/L(r); 
	end
	if j == N+1
		R(j) = w(l)*L(l)/2 - M(j)/L(l) + M(j-1)/L(l); 
	end
	if j > 1 && j < N+1
		R(j) = w(l)*L(l)/2 + w(r)*L(r)/2 ... 
                      - M(j)/L(l) - M(j)/L(r) + M(j-1)/L(l) + M(j+1)/L(r); 
	end

	for k=1:nP			% loop over all concentrated point loads
		if span(k) == l		% the point load is in the left span
			R(j) = R(j) + P(k)*xL(k)/L(l);
		end
		if span(k) == r		% the point load is in the right span
			R(j) = R(j) + P(k)*xR(k)/L(r);
		end
	end

 end

 slope = zeros(1,N);
 for j=1:N				% compute the slopes  (15)

	r = j;  		% j is the span to the right of reaction j

	slope(j) =  w(r)*L(r)^3 / 24  +  M(j+1)*L(r) / 6  +  M(j)*L(r) / 3;  

	for k=1:nP			% loop over all concentrated point loads
		if span(k) == r		% the point load is in the right span
 			slope(j) = slope(j)+P(k)*xR(k)/L(r)*(L(r)^2-xR(k)^2)/6;
		end
	end

	slope(j) = -slope(j) / ( E*I(r) );

 end

 if (  abs ( sum(R) - sum (w .* L) - sum(P) ) < 1e-9  )
	disp (' yes! ')		% equilibrium check ... should be close to zero
 end

% ----------  shear, moment, slope, and deflection data and plots  -----------

 for j=1:N	% x-axis data for shear, moment, slope, and deflection diagrams
	xs(:,j) = [ 0:L(j)/157:L(j) ]' ;
 end

 for j=1:N				% j is the span number
	Vo = ( M(j) - M(j+1) ) / L(j) - w(j)*L(j)/2;		% shear at left
	V_Diag(:,j) = Vo + w(j)*xs(:,j);
	for k=1:nP			% loop over all concentrated point loads
		if span(k) == j		% the point load is in span to the right
			i1 = find(xs(:,j)<xL(k));
			i2 = find(xs(:,j)>xL(k));
			V_Diag(i1,j) = V_Diag(i1,j) - P(k)*xR(k)/L(j);
			V_Diag(i2,j) = V_Diag(i2,j) + P(k)*(1-xR(k)/L(j));
		end
	end
	M_Diag(:,j) = M(j) + cumtrapz( -V_Diag(:,j) ) * xs(2,j);
	s_Diag(:,j) = cumtrapz( M_Diag(:,j) ) * xs(2,j) / (E*I(j)) + slope(j) ;
	d_Diag(:,j) = cumtrapz( s_Diag(:,j) ) * xs(2,j) ;
 end

% ------------- display key results to the screen --------------
 fprintf('------------------------------------------------------------------\n');
 fprintf('                 Moment              Shear             Deflection\n');
 fprintf('   Maximum      %12.5e        %12.5e      %12.5e\n', ...
             max(max(M_Diag)), max(max(-V_Diag)), max(max(d_Diag)) );
 fprintf('   Minimum      %12.5e        %12.5e      %12.5e\n', ...
             min(min(M_Diag)), min(min(-V_Diag)), min(min(d_Diag)) );
 fprintf('------------------------------------------------------------------\n');
 
 for j=1:N				% j is the span number
        V(1,j) = -V_Diag(1,j);		% shear force at left  end of span
        V(2,j) = -V_Diag(158,j);	% shear force at right end of span
 end

 for j=2:N		% x-axis data for shear and moment diagram plots
	xs(:,j) = xs(:,j) + sumL(j-1);
 end

% Plotting

 xs      =  xs(:);
 M_Diag  =  M_Diag(:);
 V_Diag  = -V_Diag(:);
 s_Diag  =  s_Diag(:);
 d_Diag  =  d_Diag(:);
 z = zeros(1,length(xs));

 D_max  = max(abs(d_Diag));

 figure(1)
  clf
  subplot(411)	 	
   plot ( xs, z, '-k', xs, V_Diag, '-b', 'LineWidth', 2 )	
    ylabel('Internal Shear')

  subplot(412)  		
   plot ( xs, z, '-k', xs, M_Diag, '-b', 'LineWidth', 2 )	
    ylabel('Internal Moment')
 
  subplot(413)	 
   plot ( xs, z, '-k', xs, s_Diag, '-b', 'LineWidth', 2 )
    ylabel('Slope')
 
  subplot(414)  	
   plot ( xs, z, '-k', xs, d_Diag, '-b', 'LineWidth', 2 )
    ylabel('Deflection')

% ----------------------------------------------------------- three_moment.m
