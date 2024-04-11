function z = drawTruncatedNormal(N,mu,sig,xlo,xhi)
%DRAWTRUNCATEDNORMAL Draws from a truncated normal distribution.
%
% Uses the inversion method to generate random numbers from a truncated
% normal distribution.
%
% Parameters:
%   - N: Size of the resulting array of random numbers
%       Note: If N is a scalar, then the result will be a NxN-matrix
%   - mu: Scalar, mean of underlying normal distribution
%   - sig: Scalar, standard deviation of underlying normal distribution
%   - xlo: Scalar, low truncation point, if any
%   - xhi: Scalar, high truncation point, if any
%
% Returns a matrix with size(N) of random truncated normal numbers.

    % Defaults
    if (nargin<2)||isempty(mu)
      mu=0;
    end

    if (nargin<3)||isempty(sig)
      sig=0;
    end

    if (nargin<4)||isempty(xlo)
      xlo=-inf;
      plo=0;
    else
      plo = 0.5 * erfc(-(xlo-mu)/sig ./ sqrt(2));
    end

    if (nargin<5)||isempty(xhi)
      xhi=inf;
      phi=1;
    else
      phi = 0.5 * erfc(-(xhi-mu)/sig ./ sqrt(2));
    end

    % Test if trunation points are reversed
    if xlo > xhi
      error('MC', 'Must have xlo <= xhi if both provided');
    end

    % Generate uniform [0,1] random deviates
    r = rand(N);

    % Scale to [plo,phi]
    r = plo + (phi - plo) * r;

    % Invert through standard normal
    z = -sqrt(2) .* erfcinv(2*r);

    % Apply shift and scale
    z = mu + z * sig;
end