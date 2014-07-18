function Obj = NLSObjHW4(b_0, data)
% PURPOSE: compute (the average of) the sum of squared residuals, 
%          its first derivatives wrt. theta, hessian, and the asymptotic
%          covariance matrix 
%--------------------------------------------------------------------
% USAGE: Obj = NLSobjH4(beta, data)
% where: b_0   = parameters  (1 x 3) [1 b1 b2] is the truth
%        data  = data matrix (n x 4) [y 1 x1 x2]
%--------------------------------------------------------------------
% RETURNS: Obj
%--------------------------------------------------------------------
% Reference:    Gregory, Allan W. and Veall, Michael R. (1985) "Formulating 
%               Wald Tests of Nonlinear Restrictions". Econometrica Vol. 53,
%               No. 6(Nov., 1985)
% ---------------------------------------------------------------------
% Written by Robert Ackerman, UNC Chapel Hill.
% November 24, 2013.
% data nx4 [y 1 x_1 x2)
% beta 1 x 3 [1 b1 b2]
% --------------------------------------------------------------------
Obj = (data(:,1) - data(:,2:end)*b_0')'*(data(:,1) - data(:,2:end)*b_0');


