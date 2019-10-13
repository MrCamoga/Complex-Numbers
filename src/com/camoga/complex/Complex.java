package com.camoga.complex;

import static java.lang.Math.atan2;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Complex {

//	public static final Complex ONE = new Complex(1,0);
//	public static final Complex ZERO = new Complex(0,0);
	public static final Complex I = new Complex(0,1);
	public static final Complex EULERGAMMA = new Complex(0.5772156649015328606065,0);
	public static final Complex PI = new Complex(Math.PI, 0);
	public static final Complex INFINITY = new Complex(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
	
	private double r, i;
	
	public Complex() {
		this(0,0);
	}
	
	public Complex(double r, double i) {
		this.r = r;
		this.i = i;
	}
	
	public static Complex add(Complex ...zs) {
		Complex result = valueOf(0);
		for(Complex z : zs) {
			result.r += z.r;
			result.i += z.i;
		}
		return result;
	}
	
	public static Complex sub(Complex w, Complex z) {
		return new Complex(w.r-z.r,w.i-z.i);
	}
	
	public static Complex mul(Complex ...zs) {
		Complex r = valueOf(1);
		for(Complex z : zs) {
			r = new Complex(r.r*z.r-r.i*z.i, r.r*z.i+z.r*r.i);
		}
		return r;
	}
	
	public static Complex div(Complex w, Complex z) {
		return div(mul(w, conjugate(z)), modSq(z));
	}
	
	public static Complex div(Complex z, double d) {
		return new Complex(z.r/d, z.i/d);
	}
	
	/**
	 * a-bi
	 * @param z
	 * @return
	 */
	public static Complex conjugate(Complex z) {
		return new Complex(z.r,-z.i);
	}
	
	/**
	 * a^2+b^2
	 * 
	 * @param z
	 * @return
	 */
	public static double modSq(Complex z) {
		return z.r*z.r+z.i*z.i;
	}
	
	/**
	 * principal value of ln(z)
	 * @param z
	 * @return
	 */
	public static Complex ln(Complex z) {
		return new Complex(Math.log(mod(z)),argument(z));
	}
	
	/**
	 * e^ix = cos(x) + isin(x)
	 * @param x argument
	 * @return
	 */
	public static Complex euler(double x) {
		return new Complex(Math.cos(x), Math.sin(x));
	}
	
	/**
	 * e^iz
	 * @param z complex number
	 * @return
	 */
	public static Complex euler(Complex z) {
		return exp(mul(z, I));
	}
	
	/**
	 * e^z
	 * @param z
	 * @return
	 */
	public static Complex exp(Complex z) {
		return mul(valueOf(Math.exp(z.r)), euler(z.i));
	}	
	
	/**
	 * cos(z)
	 * @param z
	 * @return
	 */
	public static Complex cos(Complex z) {
		return new Complex(Math.cos(z.r)*Math.cosh(z.i), -Math.sin(z.r)*Math.sinh(z.i));
	}
	
	/**
	 * sin(z)
	 * @param z
	 * @return
	 */
	public static Complex sin(Complex z) {
		return new Complex(Math.sin(z.r)*Math.cosh(z.i), Math.cos(z.r)*Math.sinh(z.i));
	}
	
	public static Complex tan(Complex z) {
		return div(sin(z), cos(z));
	}
	
	public static Complex cosh(Complex z) {
		return cos(mul(I, z));
	}

	public static Complex sinh(Complex z) {
		return mul(new Complex(0,-1), sin(mul(I, z)));
	}
	
	public static Complex tanh(Complex z) {
		return mul(new Complex(0,-1), tan(mul(I, z)));
	}
	
	public static Complex acos(Complex z) {
		return sub(valueOf(Math.PI/2), asin(z));
	}
	
	public static Complex asin(Complex z) {
		return mul(negate(PI), ln(add(mul(z,PI), sqrt(sub(valueOf(1), pow(z,2))))));
	}
	
	public static Complex atan(Complex z) {
		return mul(valueOf(0.5), mul(PI, sub(ln(sub(valueOf(1),mul(PI, z))), ln(add(valueOf(1), mul(PI, z))))));
	}
	
	public static Complex asec(Complex z) {
		return acos(reciprocal(z));
	}
	
	public static Complex acsc(Complex z) {
		return asin(reciprocal(z));
	}
	
	public static Complex acot(Complex z) {
		return atan(reciprocal(z));
	}
	
	public static Complex atanh(Complex z) {
		return mul(valueOf(0.5), sub(ln(add(valueOf(1), z)), ln(sub(valueOf(1), z))));
	}
	
	public static Complex acosh(Complex z) {
		return ln(add(z, mul(sqrt(sub(z, valueOf(1))), sqrt(add(z, valueOf(1))))));
	}
	
	public static Complex asinh(Complex z) {
		return ln(add(z, sqrt(add(valueOf(1), pow(z, valueOf(2))))));
	}
	
	public static Complex acoth(Complex z) {
		return atanh(reciprocal(z));
	}
	
	public static Complex asech(Complex z) {
		return acosh(reciprocal(z));
	}
	
	public static Complex acsch(Complex z) {
		return asinh(reciprocal(z));
	}
	
	/**
	 * logb(z)
	 * @param z
	 * @param base
	 * @return
	 */
	public static Complex log(Complex z, Complex base) {
		return div(ln(z), ln(base));
	}
	
	/**
	 * z^r
	 * @param z
	 * @param real
	 * @return
	 */
	public static Complex pow(Complex z, double real) {
		return pow(z, valueOf(real));
	}
	
	/**
	 * z^w
	 * @param z
	 * @param w
	 * @return
	 */
	public static Complex pow(Complex z, Complex w) {
		double modz = mod(z);
		if(modz == 0) {
			if(w.i == 0 && w.r == 0) return valueOf(1);
			else if(w.r == -1 && w.i == 0) return INFINITY;
			else return valueOf(0);
		}
		double argz = argument(z);
		// (r*e^ia)^(c+di)= r^c*e^(iac) * e^(ln(r)*di)*e^(-ad)
		double real = Math.pow(modz, w.r)*Math.exp(-argz*w.i);
		Complex c = euler(Math.log(modz)*w.i+argz*w.r);
		return mul(c, valueOf(real));
	}
	
	/**
	 * sqrt(z)
	 * @param z
	 * @return
	 */
	public static Complex sqrt(Complex z) {
		return pow(z, 0.5);
	}
	
	/**
	 * Riemann Zeta function
	 * @param s
	 * @param it not recommended less than 100
	 * @return
	 */
	public static Complex zeta(Complex s, int it) {
		Complex result = new Complex();
		if(s.r > 0) {
			for(int i = 1; i < it; i++) {
				result = add(result, mul(valueOf(i%2==0 ? -1:1), pow(valueOf(i), mul(s, valueOf(-1)))));
			}
			Complex lastIt = add(result, mul(valueOf((it)%2==0 ? -1:1), pow(valueOf(it+1), mul(s, valueOf(-1)))));
			
			Complex factor = reciprocal(sub(valueOf(1), pow(valueOf(2), sub(valueOf(1), s))));
			result = div(add(result, lastIt), 2);
			
			result = mul(result, factor);
		} else {
			result = mul(pow(valueOf(2*Math.PI), s),reciprocal(PI),sin(mul(s, valueOf(Math.PI/2))), gamma(sub(valueOf(1), s)), zeta(sub(valueOf(1), s), it));
		}
		return result;
	}
	
	/**
	 * Gamma function
	 * @param z
	 * @return gamma(z)
	 */	
	public static Complex gamma(Complex z) {
		if(z.real() < 0.5) return div(PI, mul(sin(mul(PI,z)), gamma(sub(valueOf(1), z))));
		
		double[] q = new double[]{75122.6331530, 80916.6278952, 36308.2951477, 8687.24529705, 1168.92649479, 83.8676043424, 2.50662827511};
		
		Complex zpow = valueOf(1);
		Complex sum = valueOf(0);
		Complex prod = valueOf(1);
		
		for(int i = 0; i < q.length; i++) {
			sum = add(sum, mul(zpow, valueOf(q[i])));
			zpow = mul(zpow, z);
			
			prod = mul(prod, add(z,valueOf(i)));
		}
		
		Complex result = div(mul(sum, pow(add(z, valueOf(5.5)), add(z, valueOf(0.5))), reciprocal(exp(add(z,valueOf(5.5))))), prod);
		
		return result;
	}
	
	/**
	 * Derivative of the gamma function
	 * @param z
	 * @param it
	 * @return
	 */
	public static Complex gammaDerivative(Complex z, int it) {
		return mul(gamma(z), digamma(z,it));
	}
	
	/**
	 * LogGamma Function
	 * This function is similar to ln(gamma(z)) but with a single branch cut along the negative real axis
	 * @param z
	 * @param it
	 * @return
	 */
	public static Complex logGamma(Complex z, int it) {
		Complex sum = new Complex();
		
		for(int k = 1; k < it; k++) {
			Complex d = div(z, k);
			sum = add(sum, d, negate(ln(add(valueOf(1),d))));
		}
		
		return add(negate(ln(z)), mul(negate(EULERGAMMA), z), sum);
	}
	
	/**
	 * Digamma function
	 * @param z
	 * @param it
	 * @return
	 */
	public static Complex digamma(Complex z, int it) {
		Complex sum = new Complex();
		
		for(int k = 0; k < it; k++) {
			sum = add(sum, valueOf(1/(double)(k+1)), negate(reciprocal(add(valueOf(k),z))));
		}
		
		return add(negate(EULERGAMMA), sum);
	}
	
	public static Complex trigamma(Complex z, int it) {
		return polygamma(z, 1, it);
	}
	
	public static Complex polygamma(Complex z, int n, int it) {
		if(n == 0) return digamma(z, it);
		Complex sum = new Complex();
		Complex sign = n%2==0 ? valueOf(-1):valueOf(1);
		
		for(int k = 0; k < it; k++) {
			sum = add(sum, reciprocal(pow(add(z,valueOf(k)), n+1)));
		}
		
		return mul(sign, factorial(n), sum);
	}
	
	/**
	 * gamma(x)gamma(y)/gamma(x+y)
	 * @param x
	 * @param y
	 * @return
	 */
	public static Complex beta(Complex x, Complex y) {
		return div(mul(gamma(x), gamma(y)), gamma(add(x,y)));
	}
	
	/**
	 * 1/z
	 * @param z
	 * @return
	 */
	public static Complex reciprocal(Complex z) {
		return div(valueOf(1), z);
	}
	
	public static Complex factorial(int n) {
		double result = 1;
		for(int i = 1; i <= n; i++) {
			result *= i;
		}
		return valueOf(result);
	}
	
	/**
	 * -z
	 * @param z
	 * @return
	 */
	public static Complex negate(Complex z) {
		return new Complex(-z.r, -z.i);
	}
	
	/**
	 * new Complex(r,0)
	 * @param r
	 * @return
	 */
	public static Complex valueOf(double r) {
		return new Complex(r,0);
	}
	
	/**
	 * atan2(b/a)
	 * @param z
	 * @return
	 */
	public static double argument(Complex z) {
		return atan2(z.i, z.r);
	}
	
	public String toString() {
		String result = "";
		if(i==0) result += r;
		else if(r==0) result += (i==1 ? "":i)+"i";
		else result += r + (i >= 0 ? (i == 0 ? (""):(" + " + (i==1 ? "":i)+"i")):(" - " + (i==-1 ? "":-i)+"i"));
		return result;
	}

	/**
	 * sqrt(a^2+b^2)
	 * @param z
	 * @return
	 */
	public static double mod(Complex z) {
		return Math.sqrt(modSq(z));
	}
	
	/**
	 * Re(z)
	 * @return
	 */
	public double real() {
		return r;
	}
	
	/**
	 * Im(z)
	 * @return
	 */
	public double imaginary() {
		return i;
	}

	/**
	 * Returns a new Complex initialized to the value given by the String
	 * @param s
	 * @return
	 */
	public static Complex parseComplex(String s) { // basura
		s = s.replace(" ", "");
		
		String real = "";
		String imaginary = "";
		
		for(String a : s.split("(([+-])?\\d?i)")) real += a;
	
		Matcher matcher = Pattern.compile("(([+-])?\\d?i)").matcher(s);
		if(matcher.find()) imaginary = matcher.group(1);
		
		imaginary = imaginary.replace("i", "");
		if(imaginary.equals("")||imaginary.equals("+")) imaginary = "1";
		if(imaginary.equals("-")) imaginary = "-1";
		if(real.equals("")) real = "0";
		
		return new Complex(Double.parseDouble(real),Double.parseDouble(imaginary));
	}
	/**
	 * z^(1/3)
	 * @param z
	 * @return
	 */
	public static Complex cbrt(Complex z) {
		return pow(z, 1/3D);
	}
}