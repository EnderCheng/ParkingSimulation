package ca.uwaterloo.cheng.utils;

import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.security.NoSuchProviderException;
import java.security.SecureRandom;
import java.security.Security;

import org.bouncycastle.jce.provider.BouncyCastleProvider;

import it.unisa.dia.gas.jpbc.Element;
import it.unisa.dia.gas.jpbc.Pairing;
import it.unisa.dia.gas.jpbc.PairingParameters;
import it.unisa.dia.gas.plaf.jpbc.pairing.PairingFactory;

public class Tools {
	private static Pairing pairing = PairingFactory.getPairing("pairingparams.properties");
	private static PairingParameters pp = PairingFactory.getPairingParameters("pairingparams.properties");
	private static BigInteger p = pp.getBigInteger("r");
	
	public static Element hash_g(byte[] data)
    {
    	Security.addProvider(new BouncyCastleProvider());
    	MessageDigest mda;
    	byte [] digest = null;
		try {
			mda = MessageDigest.getInstance("SHA-256", "BC");
			digest = mda.digest(data);
		} catch (NoSuchAlgorithmException | NoSuchProviderException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println("The Hash function cannot be found: SHA-256");
			System.exit(-1);
		}
		Element g = pairing.getG1().newElementFromHash(digest, 0, digest.length).getImmutable();
		return g;
    }
	
	public static BigInteger hash_p(byte[] data)
	{
    	Security.addProvider(new BouncyCastleProvider());
    	MessageDigest mda;
    	byte [] digest = null;
		try {
			mda = MessageDigest.getInstance("SHA-256", "BC");
			digest = mda.digest(data);
		} catch (NoSuchAlgorithmException | NoSuchProviderException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println("The Hash function cannot be found: SHA-256");
			System.exit(-1);
		}
		
		BigInteger g = new BigInteger(digest).mod(p);
		return g;
		
	}
	
	public static BigInteger hash_pp(BigInteger x)
	{
    	Security.addProvider(new BouncyCastleProvider());
    	MessageDigest mda;
    	byte [] digest = null;
		try {
			mda = MessageDigest.getInstance("SHA-256", "BC");
			digest = mda.digest(x.toByteArray());
		} catch (NoSuchAlgorithmException | NoSuchProviderException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println("The Hash function cannot be found: SHA-256");
			System.exit(-1);
		}
		
		BigInteger g = new BigInteger(digest).mod(p);
		return g;
	}
	
	public static BigInteger randomZp(BigInteger p) {
		BigInteger r;

		do {
			r = new BigInteger(p.bitLength(), new SecureRandom());
		} while (r.compareTo(BigInteger.ZERO) <= 0 || r.compareTo(p) >= 0);

		return r;
	}
	
    // return a random integer in Z_n
    public static BigInteger randomZN(BigInteger n)
    {
        BigInteger r;
        
        do
        {
            r = new BigInteger(n.bitLength(), new SecureRandom());
        }
        while (r.compareTo(BigInteger.ZERO) <= 0 || r.compareTo(n) >= 0);
        
        return r;
    }
    
    // return a random integer in Z*_n
    public static BigInteger randomZStarN(BigInteger n)
    {
        BigInteger r;
        
        do
        {
            r = new BigInteger(n.bitLength(), new SecureRandom());
        }
        while (r.compareTo(n) >= 0 || r.gcd(n).intValue() != 1);
        
        return r;
    }
    
    //extended Eclid algorithm
    public static BigInteger[] gcd(BigInteger p, BigInteger q, BigInteger module) {
        if (q.mod(module) == BigInteger.ZERO)
           return new BigInteger[] { p.mod(module), BigInteger.ONE, BigInteger.ZERO };

        BigInteger[] vals = gcd(q, p.mod(q),module);
        BigInteger d = vals[0].mod(module);
        BigInteger a = vals[2].mod(module);
        BigInteger b = vals[1].subtract((p.divide(q).mod(module)).multiply(vals[2]).mod(module)).mod(module);
        return new BigInteger[] { d, a, b };
     }
}
