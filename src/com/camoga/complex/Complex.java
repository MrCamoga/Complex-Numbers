package com.camoga.complex;

import static java.lang.Math.PI;
import static java.lang.Math.atan2;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Complex {

	public static final Complex ONE = new Complex(1,0);
	public static final Complex ZERO = new Complex(0,0);
	public static final Complex I = new Complex(0,1);
	public static final Complex INFINITY = new Complex(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
	
	private double r, i;
	
	public Complex() {
		this(0,0);
	}
	
	public Complex(double r, double i) {
		this.r = r;
		this.i = i;
	}
	
	public static Complex add(Complex w, Complex z) {
		return new Complex(w.r+z.r,w.i+z.i);
	}
	
	public static Complex sub(Complex w, Complex z) {
		return new Complex(w.r-z.r,w.i-z.i);
	}
	
	public static Complex mul(Complex w, Complex z) {
		return new Complex(w.r*z.r-w.i*z.i,w.r*z.i+z.r*w.i);
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
	 * (a^2+b^2)
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
	 * e^iphi = cos(phi) + isin(phi)
	 * @param phi argument
	 * @return
	 */
	public static Complex euler(double phi) {
		return new Complex(Math.cos(phi), Math.sin(phi));
	}
	
	/**
	 * e^iz
	 * @param z complex number
	 * @return
	 */
	public static Complex euler(Complex z) {
		Complex iz = mul(z, I);
		return exp(iz);
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
		return sub(valueOf(PI/2), asin(z));
	}
	
	public static Complex asin(Complex z) {
		return mul(new Complex(0,-1), ln(add(mul(z,I), sqrt(sub(ONE, pow(z,2))))));
	}
	
	public static Complex atan(Complex z) {
		return mul(valueOf(0.5), mul(I, sub(ln(sub(ONE,mul(I, z))), ln(add(ONE, mul(I, z))))));
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
		return mul(valueOf(0.5), sub(ln(add(ONE, z)), ln(sub(ONE, z))));
	}
	
	public static Complex acosh(Complex z) {
		return ln(add(z, mul(sqrt(sub(z, ONE)), sqrt(add(z, ONE)))));
	}
	
	public static Complex asinh(Complex z) {
		return ln(add(z, sqrt(add(ONE, pow(z, valueOf(2))))));
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
	 * logb (z)
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
			if(w.i == 0 && w.r == 0) return ONE;
			else if(w.r == -1 && w.i == 0) return INFINITY;
			else return ZERO;
		}
		double argz = argument(z);
		// (r*e^ia)^(c+di)= r^c*e^(iac) * e^(ln(r)*di)*e^(-ad)
		//real product
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
	 * sum from n=1 to <b>it</b> of 1/n^z
	 * @param z
	 * @param it
	 * @return
	 */
	public static Complex riemannzetareal(Complex z, int it) {
//		if(z.r <= 1) return valueOf(0);
		Complex result = new Complex();
		if(z.r > 1) {
			for(int i = 1; i < it; i++) {
				result = add(result, pow(valueOf(i), mul(z, valueOf(-1))));
			}
		} else if(z.r > 0 && z.r < 1) {
			for(int i = 1; i < it; i++) {
				result = add(result, mul(valueOf(i%2==0 ? -1:1), pow(valueOf(i), mul(z, valueOf(-1)))));
			}
			Complex lastIt = add(result, mul(valueOf((it)%2==0 ? -1:1), pow(valueOf(it+1), mul(z, valueOf(-1)))));
			
			Complex factor = reciprocal(sub(ONE, pow(valueOf(2), sub(ONE, z))));
			result = div(add(result, lastIt), 2);
			
			result = mul(result, factor);
		}
		return result;
	}
	
	/**
	 * 
	 * @param n
	 * @return gamma(n)
	 */
	public static Complex gamma(Complex n) {
		double dt = 0.001;
		n = sub(n, ONE);
		
		Complex dexp = exp(valueOf(-dt));
		Complex exp = ONE;
		
		Complex sum = ZERO;
		
		for(double t = 0; t < 30; t+=dt) {
			sum = add(sum, mul(exp,pow(valueOf(t), n)));
			exp = mul(exp, dexp);
//			System.out.println(t);
		}
		sum = mul(sum, valueOf(dt));

		return sum;		
	}
	
	/**
	 * 1/z
	 * @param z
	 * @return
	 */
	public static Complex reciprocal(Complex z) {
		return div(ONE, z);
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
		else result += (r) + (i >= 0 ? (i == 0 ? (""):(" + " + (i==1 ? "":i)+"i")):(" - " + (i==-1 ? "":-i)+"i"));
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
	
	public double real() {
		return r;
	}
	
	public double imaginary() {
		return i;
	}

	/**
	 * Returns a new double initialized to the value given by the String
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