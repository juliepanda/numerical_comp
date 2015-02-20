import numpy as np 
import scipy

def bisection (a , b, xtol, maxit):
	# could improve rounding as a param
	"reliable but slow 1D zero-finding: returns x where f(x)=0, xtol defines precision"

	xtol = np.float32(10 ** -xtol);
	counter = 0;

	fa = f(a);
	fb = f(b);


	if fa == np.float32(0.0):
		print 'k =', counter,', a =', str(np.around(a, decimals=16)), ', b =', str(np.around(b, decimals=16)),', |ak-bk| =', str(np.around(np.absolute(b-a), decimals=16)) , ', f(a) =',str(np.around(fa, decimals=16)),', f(b) =',str(np.around(fb, decimals=16))
		return a
	if fb == np.float32(0.0):
		print 'k =', counter,', a =', str(np.around(a, decimals=16)), ', b =', str(np.around(b, decimals=16)),', |ak-bk| =', str(np.around(np.absolute(b-a), decimals=16)) , ', f(a) =',str(np.around(fa, decimals=16)),', f(b) =',str(np.around(fb, decimals=16))
		return b

	if sign(fa) == sign(fb):
		print('Error: f(a) and f(b) have same sign.');
		return;

	while counter != maxit:
		fa = f(a);
		fb = f(b);

		print 'k =', counter,', a =', str(np.around(a, decimals=16)), ', b =', str(np.around(b, decimals=16)),', |ak-bk| =', str(np.around(np.absolute(b-a), decimals=16)) , ', f(a) =',str(np.around(fa, decimals=16)),', f(b) =',str(np.around(fb, decimals=16))
		m = (a + b)/2.0;
		fm = f(m);
		counter = counter + 1;
		g = (b-a)/2.0
		if g < xtol and g > -xtol:
			return m
		else:
			if sign(fm)!=sign(fa):
				b = m;
			else:
				a = m;
	return "ran " + str(maxit) + ' times'

def newton (x0, ftol, maxit):
	"needs single point,fast, but f' must exists"
	# xk_1 = xk - f(xk)/f'(xk)
	counter = 0;
	ftol = np.float32(10 ** -ftol);
	fx = float(f(x0));

	if np.absolute(fx) < ftol:
		return x0;
	print 'k =', counter,', x =',str(x0),', f(x) =',fx;
	while counter < maxit and np.absolute(fx) > ftol:
		counter = counter + 1;
		fx = float(f(x0));
		fpx = float(myfunc_derive(x0));

		x1 = float(x0 - (fx/fpx));
		fx1 = float(f(x1));
		print 'k =', counter,', x =',str(x1),', f(x) =',fx1;
		if np.absolute(fx) < ftol:
			# found zero
			return x1;
		else:
			x0 = x1;
	return "ran " + str(maxit) + " times";

def secant (x, x_1, ftol, maxit):
	"safer method than newton, faster method than bisection, has problem with extrapolation when interval is large"
	# x1 = x - ((f(x)(x-x_1))/(f(x)-f(x_1))
	ftol = np.float32(10 ** -ftol);
	counter = 0;

	if np.absolute(f(x)) > ftol:
		print 'k =', counter,', x =',str(x),', f(x) =',f(x);
	if np.absolute(f(x_1)) > ftol:
		print 'k =', counter,', x =',str(x_1),', f(x) =',f(x_1);

	if sign(f(x)) == sign(f(x_1)):
		print('Error: f(x) and f(x_1) have same sign.');
		return;

	while counter < maxit and np.absolute(f(x)) > ftol:
		counter = counter + 1;
		x_temp = x - ((f(x)*(x - x_1))/(f(x)-f(x_1)));

		print 'k =', counter,', x =',str(x),', f(x) =',f(x_temp);
		if np.absolute(f(x_temp)) <= ftol:
			x_1 = x;
			x = x_temp;
	return "ran " + str(maxit) + " times";


def wheeler(x0, x1, ftol, maxit):
	"modified secant method, no problem with extrapolation"
	ftol = np.float32(10 ** -ftol)
	if f(x0) == 0.0 or f(x1) == 0.0 or f(x0)*f(x1)>0.0: return
 	mu = 1
	xlst = [0]*(maxit+2)
	xlst[0] = x0
	xlst[1] = x1
	k = 0
	k_star = 0
	while np.absolute(f(xlst[k])) > ftol and k < maxit:
		k = k + 1
		xlst[k+1] = xlst[k] - f(xlst[k])*(xlst[k] - xlst[k_star])/(f(xlst[k]) - mu*f(xlst[k_star]))
		print 'k =', k,', x =',str(xlst[k+1]),', f(x) =',f(xlst[k+1]);

		if f(xlst[k+1])*f(xlst[k]) < 0:
			mu = 1; k_star = k
		else:
			mu = mu/2

	return xlst[k], f(xlst[k])

def sign (val):
	"helper function that returns the sign of a float value"
	if val < 0:
		return '-';
	elif val > 0:
		return '+';
	else:
		return '0';

def f (x):
	return x**3 - 2.5*x - 4.011

