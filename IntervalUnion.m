function [left, right] = IntervalUnion(left, right)
% function [lower upper] = IntervalUnion(left, right)
%
% Purpose: Interval merging
%
% Given N input closed intervals in braket form:
%   Ii := [left(i),right(i)], i = 1,2...,N (mathematical notation)
% The set union{Ii) can be written as a canonical partition by
%   intervals Jk; i.e., union{Ii) = union(Jk), where Jk are M intervals
%   (with M<=N, so the partition is minimum cardinal), and {Jk} are
%   disjoint to each other (their intersections are empty). This function
%   returns Jk = [lower(k),upper(k)], k=1,2,...M, in the ascending sorted
%   order.
%
% EXAMPLE USAGE:
%   >> [lower upper] = IntervalUnion([0 1 2 3 4],[1.5 1.6 3.5 3 5])
%   	lower =   0    2  4
%       upper = 1.6  3.5  5
%
% Algorithm complexity: O(N*log(N))
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Original: 19-April-2013

[left, right] = MergeBrackets(left, right);

end