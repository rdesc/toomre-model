function [a] = nbodyaccn(m, r)
% Computes the n-body accleration
%
% Input arguments
%
%     m:   (N x 1 array) Particle masses
%     r:   (N x 3 array) Particle position vectors
%  
% Return value
% 
%      a:  (N x 3 array) acclerations for all N particles in x, y, z
%      directions

s = size(m);
N = s(1);

% get the number of particles that have mass
M = length(m(m>0));

% init acceleration array
a = zeros(N, 3);

% iterate over all particles
for i = 1:N
    a_i = zeros(1,3);
    % iterate over all other particles except for i
    for j = 1:M
        if j == i
            continue
        end
        a_i = a_i + m(j) * (r(j,:) - r(i,:)) / ((r(j,1) - r(i,1)).^2 + ...
                                                (r(j,2) - r(i,2)).^2 + ...  
                                                (r(j,3) - r(i,3)).^2).^(3/2);
    end
    a(i,:) = a_i;
end    