function I = hoi_ent_ggg_COV(C,modelorder,biascorrection,psiterms,dterm)
% ENT_G Entropy of a Gaussian variable in bits
%   H = ent_g(x) returns the entropy of a (possibly
%   multidimensional) Gaussian variable x with bias correction.
%   Rows of x correspond to samples, columns to dimensions/variables.
%   (Samples first axis)
%
%   biascorrect : true / false option (default true) which specifies
%   whether bias correction should be applied to the estimated entropy.


Nvarx = 1;
Nvary = size(C,1)-modelorder-Nvarx;
Nvarz = modelorder;

Nvarxyz = Nvarx+Nvary+Nvarz;
Nvarxz = Nvarx+Nvarz;
Nvaryz = Nvary + Nvarz;

xidx = 1;
zidx = 2:modelorder+1;
yidx = zidx(end)+1:size(C,1);

Cz = C(zidx,zidx);
Cxz = C([xidx zidx], [xidx zidx]);
Cyz = C([yidx zidx], [yidx zidx]);
Cxyz = C([xidx yidx zidx], [xidx yidx zidx]);

chCz = chol(Cz);
chCxz = chol(Cxz);
chCyz = chol(Cyz);
chCxyz = chol(Cxyz);

% entropies in nats
% normalisations cancel for cmi
HZ = sum(log(diag(chCz))); % + 0.5*Nvarz*log(2*pi*exp(1));
HXZ = sum(log(diag(chCxz))); % + 0.5*(Nvarx+Nvarz)*log(2*pi*exp(1));
HYZ = sum(log(diag(chCyz))); % + 0.5*(Nvary+Nvarz)*log(2*pi*exp(1));
HXYZ = sum(log(diag(chCxyz))); % + 0.5*(Nvarx+Nvary+Nvarz)*log(2*pi*exp(1));

if biascorrection
    HZ = (HZ - Nvarz*dterm - sum(psiterms(1:Nvarz)));
    HXZ = (HXZ - Nvarxz*dterm - sum(psiterms(1:Nvarxz)));
    HYZ = (HYZ - Nvaryz*dterm - sum(psiterms(1:Nvaryz)));
    HXYZ = (HXYZ - Nvarxyz*dterm - sum(psiterms(1:Nvarxyz)));
end

% convert to bits
ln2 = log(2);
I = (HXZ + HYZ - HXYZ - HZ) / ln2;

end