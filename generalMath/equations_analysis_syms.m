%% analysing some eqations using symbolic toolbox
%
%
%
%% eqations
syms d1 d2 d3 alpha beta theta 
syms dDist1 dDist2 dDist3 dAlpha dBeta dTheta 

% cosine law in a triangle in eucleadian 2d
% d3 = d1^2 + d2^2 - 2*d1*d2*cos(theta);

% re-presenting the cosine law as a function thea can be compared with a zero
F = (d1+dDist1)^2 + (d2+dDist2)^2 - 2*(d1+dDist1)*(d2+dDist2)*cos(theta+dTheta) - (d3+dDist3);
diff(F,dTheta)
%diff(F,dTheta) = sin(dTheta + theta)*(d2 + dDist2)*(2*d1 + 2*dDist1)

accucacy_dist = 0.001;
int(F,dTheta)