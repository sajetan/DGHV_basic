//Tejas Narayana

clear;

//implementation of DGHV: https://eprint.iacr.org/2009/616.pdf
function DGHVKeyGen(eta, gamma, N)
	repeat p := Random(Round(2^(eta-1)), Round(2^(eta))); until ((p mod N) ne 0) and  ((p mod 2) ne 0); //considering for N>2
	repeat q0 := Random(Round((2^gamma) div p)); until (q0 mod 2) ne 0;
	x0 := p*q0;
	return p, x0;
end function;

function DGHVEncrypt(m, N, p, rho, gamma)
	r := Random(-2^rho, 2^rho) ;
	q := Random(Round((2^gamma) div p)) ;
	c := p*q + N*r + m ;
	return c,(r);
end function;

//decryption function, choosing c_prime by c modulo p and subtracting by p such that [c]p
function DGHVDecrypt(c, N, p)
	c_prime := c mod p;
	if c_prime ge (p div 2) then
		c_prime := c_prime - p;
	end if;
 	
	m := c_prime mod N;
	return m;
end function;

function DGHVAdd(c1, c2, x0)
	c := (c1 + c2) mod x0;
	return c;
end function;

function DGHVMult(c1, c2, x0)
	c := (c1 * c2) mod x0;
	return c;
end function;


//similar to the decryption function, extracting noise term by Nr+m-m div N and returning the absolute value
function DGHVNoiseTerm(c, N, p)
	c_prime := c mod p;
	if c_prime ge (p div 2) then
                c_prime := c_prime - p;
        end if;
	r := (c_prime - (c_prime mod N)) div N;
	return (r);
end function;


