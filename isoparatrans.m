function [x,y,chi_eval] = isoparatrans(eps,eta,x_nodes,y_nodes)
%UNTITLED2 Maps a point (eps,eta) in quadrilateral reference element 
%( [-1,1]x[-1,1] ) to a point in an arbitary quadrilateral element with 
% corners (nodes) defined by x_nodes, y_nodes.

% probably not perfect by any means, just threw it together to get a feel
% for isoparametric transformations.

% chi(j) is a function equal to 1 at node j and 0 at all other nodes (in
% reference element).
chi = @(eps, eta) [(eps-1)*(eta-1)*0.25, -(eps+1)*(eta-1)*0.25, (eps+1)*(eta+1)*0.25, -(eps-1)*(eta+1)*0.25];
chi_eval = chi(eps,eta);

x = chi_eval*x_nodes;
y = chi_eval*y_nodes;

