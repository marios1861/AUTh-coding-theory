% This work is based on the paper:
% T. Richardson and R. Urbanke. "Efficient Encoding of Low-Density
% Parity-Check Codes". IEEE Trans. on Information theory, Feb. 2001


clear variables;
close all;
clc

% *******************************************************
% Taken from IEEE 802.11n: HT LDPC matrix definitions
% You can change this according to your needs :)
Z = 27;
rotmatrix = ...
    [0 -1 -1 -1 0 0 -1 -1 0 -1 -1 0 1 0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;
    22 0 -1 -1 17 -1 0 0 12 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1 -1;
    6 -1 0 -1 10 -1 -1 -1 24 -1 0 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1;
    2 -1 -1 0 20 -1 -1 -1 25 0 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1;
    23 -1 -1 -1 3 -1 -1 -1 0 -1 9 11 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1;
    24 -1 23 1 17 -1 3 -1 10 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1;
    25 -1 -1 -1 8 -1 -1 -1 7 18 -1 -1 0 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1;
    13 24 -1 -1 0 -1 8 -1 6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1;
    7 20 -1 16 22 10 -1 -1 23 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1;
    11 -1 -1 -1 19 -1 -1 -1 13 -1 3 17 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1;
    25 -1 8 -1 23 18 -1 14 9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0;
    3 -1 -1 -1 16 -1 -1 2 25 5 -1 -1 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0];

H = zeros(size(rotmatrix)*Z);
Zh = diag(ones(1,Z),0);

% Convert into binary matrix
for r=1:size(rotmatrix,1)
    for c=1:size(rotmatrix,2)
        rotidx = rotmatrix(r,c);
        if (rotidx > -1)            
            Zt = circshift(Zh,[0 rotidx]);
        else
            Zt = zeros(Z);
        end
        limR = (r-1)*Z+1:r*Z;
        limC = (c-1)*Z+1:c*Z;
        H(limR,limC) = Zt;
    end
end



% *******************************************************
% Don't change much from this line below please

% Define sizes of message and code
mLen = size(H,1);
cLen = size(H,2);

% Create a random message
msg = round(rand(1,mLen));


% Rename variables to match those used by the author
n = cLen;
m = cLen - mLen;
s = msg;
% Let c = [s p1 p2] be the codeword


% *******************************************************
% Check table III on the paper for a general description of 
% this algorithm:
Hrow1 = H(:,end);

% Find the 'gap' length
for i=1:cLen
    if Hrow1(i) == 1
        g = i;
        break;
    end
end
g = mLen - g;

% Width and last index of submatrices A, B
wA = n-m;
wB = g;
eA = wA;
eB = wA + wB;

% Extract the submatrices A, B, C, D, E and T
A = H(1:m-g,1:eA);
B = H(1:m-g,eA+1:eB);
T = H(1:m-g,eB+1:end);
C = H(m-g+1:end,1:eA);
D = H(m-g+1:end,eA+1:eB);
E = H(m-g+1:end,eB+1:end);


% Calculate p1 and p2
invT = (inv(T)); % or abs(inv(T)) ?
ET1 = -(E*invT);

Iup = diag(ones(1,size(ET1,2)),0);
Idn = diag(ones(1,size(ET1,1)),0);
X = [Iup zeros(size(Iup,1),size(Idn,2)); ET1 Idn];
Y = X*H;
spy(Y); % If exists a null space @ lower right corner, then we're ok

phi = ET1*B + D;
xtra = ET1*A + C;
p1 = mod(phi*xtra*(s'),2)';
p2 = mod(invT*(A*(s') + B*(p1')),2)';
c = [s p1 p2];

% Checking c*H'=0;
zero = mod(c*H',2);
if sum(zero)== 0
    % disp('OK!')
else
    disp('Error')
end
