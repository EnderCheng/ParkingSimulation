package ca.uwaterloo.cheng.utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.Properties;

import it.unisa.dia.gas.jpbc.Element;
import it.unisa.dia.gas.jpbc.Pairing;
import it.unisa.dia.gas.jpbc.PairingParameters;
import it.unisa.dia.gas.plaf.jpbc.pairing.PairingFactory;
import it.unisa.dia.gas.plaf.jpbc.pairing.a.TypeACurveGenerator;
import it.unisa.dia.gas.plaf.jpbc.util.io.Base64;

/**
 * @author Cheng.Huang
 * @ClassName: ParameterGenerator
 * @Description: The class is for generating the parameters for setup phase
 * @date 09/20/2017
 */
public class ParametersGenerator {

    public static void main(String[] args) {
    	System.out.println("Pairing Parameters Generating...");
    	
        int rBits = 160;
        int qBits = 512;
        TypeACurveGenerator pg = new TypeACurveGenerator(rBits,qBits);
        PairingParameters pp = pg.generate();
        
        String pairing_file_name = "pairingparams.properties";
        File pairing_file = new File (pairing_file_name);
        try(PrintWriter out = new PrintWriter(pairing_file)) {
			out.write(pp.toString());
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.out.println("The file is not found: pairingparams.properties");
			System.exit(-1);
		}
        System.out.println("Success");
        /************************************************************/
        
        System.out.println("Other Public Parameters Generating...");
        BigInteger p =pp.getBigInteger("r");
        Pairing pairing = PairingFactory.getPairing("pairingparams.properties");
        Element g = pairing.getG1().newRandomElement().getImmutable();
        Element gt = pairing.pairing(g, g);
        
        BigInteger x = BigInteger.ONE;
        do{
        	x = new BigInteger(rBits, new SecureRandom());
        }while(x.compareTo(p)!=-1);
        Element X = g.pow(x);
        
        BigInteger y = BigInteger.ONE;
        do{
        	y = new BigInteger(rBits, new SecureRandom());
        }while(y.compareTo(p)!=-1);
        Element Y = g.pow(y);
        
        BigInteger z = BigInteger.ONE;
        do{
        	z = new BigInteger(rBits, new SecureRandom());
        }while(z.compareTo(p)!=-1);
        Element Z = g.pow(z);
        
        BigInteger mu = BigInteger.ONE;
        do{
        	mu = new BigInteger(rBits, new SecureRandom());
        }while(mu.compareTo(p)!=-1);
        
        BigInteger a = BigInteger.ONE;
        do{
        	a = new BigInteger(rBits, new SecureRandom());
        }while(a.compareTo(p)!=-1);
        
        Element A = g.pow(a);
        
        String params_file_name = "params.properties";
        String params_file_name_private = "params.properties_private";
        File params_file = new File (params_file_name);
        File params_file_private = new File (params_file_name_private);
        Properties pts = new Properties();
        pts.setProperty("p", p.toString());
        pts.setProperty("g", Base64.encodeBytes(g.toBytes()));
        pts.setProperty("gt", Base64.encodeBytes(gt.toBytes()));
        pts.setProperty("X", Base64.encodeBytes(X.toBytes()));
        pts.setProperty("Y", Base64.encodeBytes(Y.toBytes()));
        pts.setProperty("Z", Base64.encodeBytes(Z.toBytes()));
        pts.setProperty("mu", mu.toString());
        pts.setProperty("A", Base64.encodeBytes(A.toBytes()));
        try(PrintWriter out = new PrintWriter(params_file)) {
        	pts.store(out, null);
        } catch (IOException e) {
			e.printStackTrace();
			System.out.println("The file cannot be saved: params.properties");
			System.exit(-1);
		}
        pts.setProperty("a", a.toString());
        try(PrintWriter out = new PrintWriter(params_file_private)) {
        	pts.store(out, null);
        } catch (IOException e) {
			e.printStackTrace();
			System.out.println("The file cannot be saved: params.properties_private");
			System.exit(-1);
		}
        
//        Properties prop = new Properties();
//        try(FileReader reader = new FileReader(params_file_name)) {
//        	prop.load(reader);
//		} catch (IOException e) {
//			e.printStackTrace();
//			System.out.println("The file cannot be loaded: params.properties");
//			System.exit(-1);
//		}
        
        System.out.println("Success");
    }
}
