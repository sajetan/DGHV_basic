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

//ACD problem: https://eprint.iacr.org/2016/215.pdf
//implementation of OLA based attack as given in Section4
function LatAttackDGHV(eta,rho, cs)
        R := 2^rho;
        t := #cs;
        B := Matrix(Integers(), t, t+1, []); //deriving the basis matrix
        for i in [1..#cs] do
            B[i][1] := cs[i];
            B[i][i+1] := R;
        end for;
        Lat,Coef := BKZ(LatticeWithBasis(B),30);
        Coef := Transpose(Coef);
        smallMatrix := (Submatrix(Coef, 1,1,t,t-1)); //taking t-1 equations like mentioned in the paper
        q_vectors := Basis(Kernel(smallMatrix)); //finding the kernel of the matrix that maps to the lattice
        q0 := q_vectors[1][1];
        r0 := cs[1] mod q0;
        p := (cs[1] - r0) div q0;
        //if r0 is even p is even (we assumed p is odd), so adding 1 to match the parity
        if (p mod 2 eq 0) then
            p:=p+1;
        end if;
        return p;
end function;



//implementation of SDA based attack as given in Section3
function LatAttackDGHV2(eta,rho, cs)
        R := 2^(rho+1);
        t := #cs;
        B := Matrix(Integers(), t+1, t+1, [0: i in [1..(t+1)^2]]); //deriving the basis matrix
        for i in [1..#cs] do
            B[1][i+1] := cs[i];
            B[i+1][i+1] := -cs[1];
        end for;
        B[1][1]:=R;
        Lat := BKZ(LatticeWithBasis(B),2);
        q0 := Basis(Lat)[1][1] div R;
        r0 := cs[1] mod q0;
        p := (cs[1] - r0) div q0;
        //if r0 is even p is even (we assumed p is odd), so adding 1 to match the parity
        if (p mod 2 eq 0) then
            p := p+1;
        end if;
        return p;
end function;

