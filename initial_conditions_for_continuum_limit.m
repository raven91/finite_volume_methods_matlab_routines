% mu = [0.5, pi];
% sigma = [0.2, 0; 0, 0.2];
% n = 100;
% l = 100;
% [x, phi] = meshgrid(0 : 1/n : 1, 0 : 2*pi/l : 2*pi);
% x = x(:);
% phi = phi(:);
% u = [x, phi];
% 
% f = zeros(size(u, 1), 1);
% for i = 1 : size(u, 1)
%     f(i, 1) = sqrt(det(inv(sigma))) / (2 * pi) * exp(-0.5 * (u(i, :) - mu) * inv(sigma) * transpose(u(i, :) - mu));
% end

%%
% figure;
% surf(reshape(x, [n + 1, l + 1]), reshape(phi, [n + 1, l + 1]), reshape(f, [n + 1, l + 1]), 'EdgeColor', 'none');
% axis equal;
% 
% phi = linspace(0, 1, l);
% mu = 0.5;
% sigma = 0.4;
% f = 1 / sqrt(2 * pi * sigma * sigma) * exp(-0.5 * (phi - mu) .* (phi - mu) / sigma / sigma);
% figure;
% plot(phi, f, '-');

%%
% phi = linspace(0, 1, l);
% f = (exp(-10 * phi .^ 2) / sqrt(0.1 * pi)) .^ 4;
% figure;
% plot(phi, f, '-');

%%
% phi = linspace(-5, 5, l);
% f = (exp(-0.5 * (phi - 2) .^ 2) + exp(-0.5 * (phi + 2) .^ 2)) / (2 * sqrt(2 * pi));
% figure;
% plot(phi, f, '-');

%%
l = 100;
% dl = 2 * pi / l;
dl = 1 / l;
x = - 0 * (l - 1) : dl : dl * (1 * l - 1);
kappa = 4;
mu = 0.5;
% von Mises
% y = (0.0 + exp(kappa * cos(x - mu)) / (2 * pi * besseli(0, kappa))) / (0.0 + 1);%0~2pi
y = (0.0 + exp(kappa * cos(x*(2*pi) - mu)) / (2 * pi * besseli(0, kappa))) / (0.0 + 1) * (2*pi);%0~1
% normal
% kappa = 0.8;
% y = (0.1 + exp(-((x - mu) .^ 2) / (2 * kappa ^ 2)) / sqrt(2 * pi * kappa ^ 2)) / (0.1 + 1);
% Cauchy
% a = 1;
% y = 1 / (2 * pi) * (1 + a * cos(x - mu));
% x_0 = mu;
% gamma = 0.4;
% y = 1 / (pi * gamma) * ((gamma ^ 2) ./ ((x - x_0) .^ 2 + gamma ^ 2));
% wrapped Cauchy
% y = 1 / (2 * pi) * sinh(gamma) ./ (cosh(gamma) - cos(x - x_0));
% wrapped Cauchy
% rho_0 = 0.68;
% y = 1 / (2 * pi) * (1 - rho_0 ^ 2) ./ (1 + rho_0 ^ 2 - 2 * rho_0 * cos(x - mu));%0~2pi
% y = 1 / (2 * pi) * (1 - rho_0 ^ 2) ./ (1 + rho_0 ^ 2 - 2 * rho_0 * cos(2 * pi * (x - mu))) * (2 * pi);%0~1

% subplot(3, 1, 1);
% plot(x, y, '-*');
% sum(y) * dl

% sine-skewed modification
lambda = 0.5;
k = 1;
% y = y .* (1 + lambda * sin(k * (x - mu)));%0~2*pi
% y = y .* exp(x + cos(1 * (x - mu)));%0~2*pi
% y = y .* (1 + 1 * besseli(1, x - mu));
% y = y .* (1 + lambda * sin(2 * pi * k * (x - mu)));%0~1
sum(y) * dl

% subplot(3, 1, 2);
% plot(x, (1 + lambda * sin(k * (x - mu))), '-*');

% figure;
% hold on;
% subplot(3, 1, 3);
plot(x, y, '-*');