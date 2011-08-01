function H_demag = ellipsoid_demag(a, b, M)
% Calculate the demagnetising field of a single domain magnetic ellipsoid
% of rotation with magnetisation = M along
% long axis with radius a, other axes radius b, ellipsoid long axis is
% assumed to be along the z axis.

% Calculate de-mag "tensor" (in co ordinates s.t. z-axis is along a)
q = a./b;
N(3) = 1./(q.^2 - 1).* (q./sqrt(q.^2 - 1)) .* log(q + sqrt(q.^2 - 1) - 1); % forumla from Magnetism, Craik
N(2) = (1 - N(3))/2;
N(1) = N(2);

% Componentwise multiply N and M to get the demag field
H_demag = -[N(1)*M(1) N(2)*M(2) N(3)*M(3)];
% [Coded directly because .* hates vectors with different orientations;
% leads to lots of crashes. ]