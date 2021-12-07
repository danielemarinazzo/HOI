function [l, r] = RangeIntersection(varargin)
% [l, r] = RangeIntersection(l1, r1, l2, r2, ... )
%
% Purpose: Intersection of n sets of interval-unions
%
% Given N sets of input closed interval-unions in braket form:
% I1 := [l1(i),r1(i)], i = 1,2...,M1
% ...
% In := [ln(j),rn(j)], j = 1,2...,Mn (mathematical notation)
% The set K = intersect(union{I1),...union{In)) can be written as a canonical
% partition by intervals Kk; i.e., union{Kk) = K, where Kk are P
% intervals and {Kk} are disjoint to each other (their intersections
% are empty).
% This function returns
% Kk = [l(k),r(k)], k=1,2,...P, in the ascending sorted order.
%
% EXAMPLE USAGE:
% >> [lower upper] = RangeIntersection([1 5], [3 9], 2, 8)
% lower = 2 5
% upper = 3 8
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Original: 19-April-2013

nset = length(varargin)/2;

if nset == 0
    l = [];
    r = [];
elseif nset == 1
    [l, r] = MergeBrackets(varargin{1:2});
else
    [l, r] = deal(varargin{1:2});
    for k = 2:nset
        [l, r] = RI2(l, r, varargin{2*k+(-1:0)});
    end
end

end % RangeIntersection

%%
function [l, r] = RI2(l1, r1, l2, r2)

l = bsxfun(@max, l1(:), l2(:)');
r = bsxfun(@min, r1(:), r2(:)');
intersect = l < r;
[l, r] = MergeBrackets(l(intersect)', r(intersect)');

end % RI2
