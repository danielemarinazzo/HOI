function H = hoi_ent_g_COV(C,biascorrection,psiterms,dterm)
% ENT_G Entropy of a Gaussian variable in bits
%   H = ent_g(x) returns the entropy of a (possibly
%   multidimensional) Gaussian variable x with bias correction.
%   Rows of x correspond to samples, columns to dimensions/variables.
%   (Samples first axis)
%
%   biascorrect : true / false option (default true) which specifies
%   whether bias correction should be applied to the estimated entropy.


if ~ismatrix(C)
    error('ent_g: input arrays should be 2d')
end

Nvarx = size(C,1);

chC = chol(C);

% entropy in nats
HX = sum(log(diag(chC))) + 0.5*Nvarx*(log(2*pi)+1);

% bias correction
if biascorrection
    HX = (HX - Nvarx*dterm - sum(psiterms));
end

% convert to bits
ln2 = log(2);
H = HX / ln2;
