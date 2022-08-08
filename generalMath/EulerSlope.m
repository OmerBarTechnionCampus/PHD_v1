function [mE,msgs] = EulerSlope(varargin)
% this function calculates the euler slope / azimuth using the other
% edges slopes / azimuths
% based upon: https://en.wikipedia.org/wiki/Euler_line#Slope
% Omer Bar 2019 July, Version 1
mE = NaN;
msgs = '';
try
    switch nargin
        case 3
            m1 = varargin{1};
            m2 = varargin{2};
            m3 = varargin{3};
        case 6
            de1 = varargin{1};
            dn1 = varargin{2};
            de2 = varargin{3};
            dn2 = varargin{4};
            de3 = varargin{5};
            dn3 = varargin{6};
            
            m1 = Azimuth(de1,dn1); % function Azimuth is in "d:\GoogleDrive\_Coding\Matlab\Geodesy\"
            m2 = Azimuth(de2,dn2);
            m3 = Azimuth(de3,dn3);
        otherwise
            mE = NaN;
            msgs = 'input must be (slope1,slope2,slop3) Or (de1,dn1,de2,dn2,de3,dn3)';
    end% switch
    mE = -(m1*m2 + m1*m3 + m2*m3 + 3)/(m1+m2+m3+3*m1*m2*m3);
catch excp
    mE = NaN;
    msgs = excp.message;
end %try
end %function EulerSlope