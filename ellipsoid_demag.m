function H_demag = ellipsoid_demag(a, b, M)
% Calculate the demagnetising field of a single domain magnetic ellipsoid
% of rotation with magnetisation = M along
% long axis with radius a, other axes radius b, ellipsoid long axis is
% assumed to be along the z axis.

% To give H_demag in units of A/m M should be in A/m.

% Calculate de-mag tensor in co ordinates s.t. z-axis is along a
q = a./b;
N(3) = 1./(q.^2 - 1).* (q./sqrt(q.^2 - 1)) .* log(q + sqrt(q.^2 - 1) - 1); % forumla from Magnetism, Craik
N(2) = (1 - N(3))/2;
N(1) = N(2);

H_demag = -N.*M;