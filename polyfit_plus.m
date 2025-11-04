function [p,S,mu] = polyfit_plus(x,y,n,p0,relto,rezi)
%polyfit_plus Fit polynomial to data.
%
%polyfit_plus is based on the native Matlab function polyfit, but has some
%extra functionality - especially for fits with less data points than order
%n plus one.
%
%   p = polyfit_plus(x,y,n) finds the coefficients of a polynomial p(x) of
%   degree n that fits the data y best in a least-squares sense. p is a
%   row vector of length n+1 containing the polynomial coefficients in
%   descending powers, p(1)*x^n + p(2)*x^(n-1) +...+ p(n)*x + p(n+1).
%   If n is a scalar with n+1 < length(x(:)), the first n+1-length(x(:))
%   elements of p will be set to zero. If n is a vector with 2 elements,
%   n = n(1) will set the degree of the polynomial, and n_fit = n(2) sets
%   the order of the polynomial fit. n(2) <= n(1) and n(2) <= length(x(:))-1
%   are required. If n is omitted, n = length(x(:))-1 is assumed.
%
%   p = polyfit_plus(x,y,n,p0) works like p = polyfit_plus(x,y,n), but a
%   default polynomial p0 is used to determine the first n - n_fit elements
%   of p, with n_fit = min(n(2), length(x(:)) - 1). If p0 is omitted or has
%   less than n(1)+1 elements, the missing elements of p0 are set to 0.
%   The fit is done in the following manner:
%     1. calculate error y1 of initial polynomial p0 at the data points:
%        y1 = y - p0(1)*x^n + p0(2)*x^(n-1) +...+ p0(n)*x + p0(n+1).
%     2. fit the error with an intermediate polynomial p1
%        p1 = polyfit_plus(x,y1,[n,n_fit])
%     3. The result is obtained by p = p0 + p1.
%   Thus, for all i <= n-n_fit:
%   p(i) = p0(i).
%
%   p = polyfit_plus(x,y,n,p0,relto) works like p = polyfit_plus(x,y,n,p0), 
%   but an additional parameter relto gives the order of the polynomial
%   coefficient, with respect to which the other coefficients of p0 are
%   assumed to be relative. I.e. with j = n+1-relto for all i <= n-n_fit:
%   p(i) = p0(i) / p0(j) * p(j).
%   If p0(j) == 0, an error will be launched.
%   If relto is empty or < 0 or > n+1, p(i) = p0(i) will be assumed as if
%   relto were omitted. 
%
%   If p0 is a matrix with 2 rows, each column is assumed to define a
%   linear regression line for the corresponding coefficient.
%   I.e. with j = n+1-relto for all i <= n-n_fit:
%   p(i) = p0(1,i) * p(j) + p0(2,i)
%   If relto is empty or < 0 or > n+1, p(i) = p0(i) will be assumed as if
%   relto were omitted. 
%
%   p = polyfit_plus(x,y,n,p0,relto,rezi) works like p = polyfit_plus(x,y,n,p0,relto), 
%   but for rezi = 1 the transformation y = 1./y is done prior to the fit.
%   Afterwards, polyval is done with n+1 data points (x2,y2) in the interval given
%   by x, and the actual result p is fitted by polyfit(x2,1./y2,n). If p0
%   ~= 0, it is transformed to its reciprocal before.
%
%   [p,S] = polyfit_plus(...) returns the polynomial coefficients p and a
%   structure S for use with polyval to obtain error estimates for
%   predictions.  S contains fields for the triangular factor (R) from a QR
%   decomposition of the Vandermonde matrix of X, the degrees of freedom
%   (df), and the norm of the residuals (normr).  If the data Y are random,
%   an estimate of the covariance matrix of p is (Rinv*Rinv')*normr^2/df,
%   where Rinv is the inverse of R.
%
%   [p,S,mu] = polyfit_plus(...) finds the coefficients of a polynomial in
%   x^ = (x-mu(1))/mu(2) where mu(1) = mean(x) and mu(2) = std(x). This
%   centering and scaling transformation improves the numerical properties
%   of both the polynomial and the fitting algorithm.


%% manage sizes of x, y and scale x, if neccessary

if ~isequal(size(x),size(y))
    error(message('MATLAB:polyfit:XYSizeMismatch'))
end

x = x(:);
y = y(:);

ind = isfinite(x) & isfinite(y);
x = x(ind);
y = y(ind);

if nargout > 2
   mu = [mean(x); std(x)];
   x = (x - mu(1))/mu(2);
end

%% determine polynomial degrees n and n_fit (n: for output, n_fit: for fit to data)

if nargin < 3 || isempty(n)
    n = length(x) - 1;
end
n = round(n);
if length(n) == 1
    n_fit = length(x) - 1;
else
    n_fit = min(n(2), length(x) - 1);
    n     = n(1);
end
n     = max(0,n);
n_fit = max(0,min(n,n_fit));

%% defaults for parameters p0,relto and rezi

% just use native Matlab function, if no other features are required:
if nargin < 4
    [p,S] =  polyfit(x,y,n_fit);
    p = [zeros(1,n-n_fit), p];
    return;
end

% ammend p0 with leading zeros, if needed.
[row,col] = size(p0);
if col == 1 && row > 1
    p0 = p0';
end
p0rev = fliplr(p0);
p0rev(:,n+2) = 0;
p0 = fliplr(p0rev(:,1:n+1));

if nargin < 5
    relto = [];
end

if nargin < 6
    rezi = 0;
end

%% reciprocal prefix: prepare data and p0 in case of rezi == 1

if rezi
    if any(p0)
        p0 = recipoly(p0,x);
    end
    y = 1./y;
end

%% calculate polynomial fit with defaults from p0

if n == n_fit
    p = polyfit(x,y,n);
elseif isempty(relto)
% absolute version of default polynomial p0
    p0 = p0(end,:);
%   1. calculate error y1 of initial polynomial p0 at the data points:
    y1 = y - polyval(p0,x); 
%   2. fit the error with an intermediate polynomial p1
    p1 = polyfit_plus(x,y1,[n,n_fit]);
%   3. The result is obtained by
    p = p0 + p1;
elseif size(p0,1) == 1
% relative version of default polynomial p0
    relto = max(0, min(n, round(relto(1))));
    J = n+1-relto;
    if p0(J) == 0
        error('p0(n+1-relto) must not be zero!');
    end
       
    % Construct Vandermonde matrix.
    V(:,n+1) = ones(length(x), 1);
%     for j = n:-1:1
    for j = n:-1:J
       V(:,j) = x.*V(:,j+1);
    end
        
    % Ammend V and y in a way to add p(i)- p(J)*p0(i)/p0(J) = 0 to the matrix
    % equations.
    i = 1:(n-n_fit);
    V((end+1):(end+n-n_fit),:) = eye(n-n_fit,n+1);
    V((length(x)+1):end,J) = -p0(i) / p0(J);
    y(end+n-n_fit) = 0;

    % Solve least squares problem.
    [Q,R] = qr(V,0);
    ws = warning('off','all'); 
    p = R\(Q'*y);    % Same as p = V\y;
    warning(ws);
    if size(R,2) > size(R,1)
         warning(message('MATLAB:polyfit:PolyNotUnique'))
    elseif warnIfLargeConditionNumber(R)
         warning('MATLAB:polyfit:RepeatedPointsOrRescale', ...
                ['Polynomial is badly conditioned. Add points with distinct X\n' ...
                 '         values, reduce the degree of the polynomial, or try centering\n' ...
                 '         and scaling as described in HELP POLYFIT.']);
    end
    
else
% relative linear regression version of default polynomial p0
    relto = max(0, min(n, round(relto(1))));
    J = n+1-relto;
       
    % Construct Vandermonde matrix.
    V(:,n+1) = ones(length(x), 1);
%     for j = n:-1:1
    for j = n:-1:J
       V(:,j) = x.*V(:,j+1);
    end
        
    % Ammend V and y in a way to add p(i)- p(J)*p0(1,i) = p0(2,i) to the matrix
    % equations.
    i = 1:(n-n_fit);
    V((end+1):(end+n-n_fit),:) = eye(n-n_fit,n+1);
    V((length(x)+1):end,J) = -p0(1,i);
    y(end+n-n_fit) = p0(2,i);

    % Solve least squares problem.
    [Q,R] = qr(V,0);
    ws = warning('off','all'); 
    p = R\(Q'*y);    % Same as p = V\y;
    warning(ws);
    if size(R,2) > size(R,1)
         warning(message('MATLAB:polyfit:PolyNotUnique'))
    elseif warnIfLargeConditionNumber(R)
         warning('MATLAB:polyfit:RepeatedPointsOrRescale', ...
                ['Polynomial is badly conditioned. Add points with distinct X\n' ...
                 '         values, reduce the degree of the polynomial, or try centering\n' ...
                 '         and scaling as described in HELP POLYFIT.']);
    end
    
end

%% calculate additional output parameters, if needed
if nargout > 1
    r = y - V*p(:);
    % S is a structure containing three elements: the triangular factor from a
    % QR decomposition of the Vandermonde matrix, the degrees of freedom and
    % the norm of the residuals.
    S.R = R;
    S.df = max(0,length(x) - (n+1));
    S.normr = norm(r);
end

p = p(:)';          % Polynomial coefficients are row vectors by convention.

%% reciprocal postfix: calculate result as reciprocal of p in case of rezi == 1
if rezi
    p = recipoly(p,x);
end

function p_rezi = recipoly(p,x)
% Calculate the reciprocal p_rezi of polynomal p with datapoints in
% interval given by x.

if all(p == 0)
    error('Reciprocal of all zero polynomial not possible!');
end
minx = min(x);
maxx = max(x);
n = length(p)-1;
number_of_points = n;
x2 = [];
while length(x2) < n + 1
    number_of_points = number_of_points + 1;
    x2 = linspace(minx,maxx,number_of_points);
    y2 = polyval(p,x2);
    x2(y2 == 0) = [];
end
y2(y2 == 0) = [];
p_rezi = polyfit_plus(x2,1./y2,n);

function flag = warnIfLargeConditionNumber(R)
if isa(R, 'double')
    flag = (condest(R) > 1e+10);
else
    flag = (condest(R) > 1e+05);
end

