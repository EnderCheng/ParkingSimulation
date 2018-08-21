package ca.uwaterloo.cheng;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Files;
import java.util.Properties;

import org.bouncycastle.util.encoders.Base64;

import ca.uwaterloo.cheng.utils.BloomFilter;
import ca.uwaterloo.cheng.utils.Tools;
import it.unisa.dia.gas.jpbc.Element;
import it.unisa.dia.gas.jpbc.Pairing;
import it.unisa.dia.gas.plaf.jpbc.pairing.PairingFactory;

public class Simulation {

	private static Pairing pairing = PairingFactory.getPairing("pairingparams.properties");
//	private static PairingParameters pp = PairingFactory.getPairingParameters("pairingparams.properties");
	private static BloomFilter<String> Omega;
	private static BigInteger p, mu, d, r, v, a;
	private static Element g, X, Y, Z, A, W;

	public static void main(String[] args) {
		int NUM = 100; // repeat times
		int[] nlist = { 10, 20, 50, 1000 };
		long commsize = 0;
		long userstore = 0;
		long serverstore = 0;
		boolean result = false;
		int expectedNumberOfElements = 10000000;
		double falsepositiverate = 0.0001;
		String params_file_name = "params.properties_private";
		Properties prop = new Properties();
		try (FileReader reader = new FileReader(params_file_name)) {
			prop.load(reader);
		} catch (IOException e) {
			e.printStackTrace();
			System.out.println("The file cannot be loaded: params.properties_private");
			System.exit(-1);
		}
		p = new BigInteger(prop.getProperty("p"));
		g = pairing.getG1().newElementFromBytes(Base64.decode(prop.getProperty("g"))).getImmutable();
//		/gt = pairing.getGT().newElementFromBytes(Base64.decode(prop.getProperty("gt"))).getImmutable();
		X = pairing.getG1().newElementFromBytes(Base64.decode(prop.getProperty("X"))).getImmutable();
		Y = pairing.getG1().newElementFromBytes(Base64.decode(prop.getProperty("Y"))).getImmutable();
		Z = pairing.getG1().newElementFromBytes(Base64.decode(prop.getProperty("Z"))).getImmutable();
		mu = new BigInteger(prop.getProperty("mu"));
		A = pairing.getG1().newElementFromBytes(Base64.decode(prop.getProperty("A"))).getImmutable();
		a = new BigInteger(prop.getProperty("a"));
		Omega = new BloomFilter<>(falsepositiverate, expectedNumberOfElements);
//		Xi = new BloomFilter<>(falsepositiverate, expectedNumberOfElements);
//		Psi = new BloomFilter<>(falsepositiverate, expectedNumberOfElements);

		// offline
		Element eYg = pairing.pairing(Y, g).getImmutable();
		Element eZA = pairing.pairing(Z, A).getImmutable();
		Element eZgmu = pairing.pairing(Z, g.pow(mu)).getImmutable();
		Element eZg = pairing.pairing(Z, g).getImmutable();
		Element eXg = pairing.pairing(X, g).getImmutable();

		// user subscribe
		d = Tools.randomZp(p);
		r = Tools.randomZp(p);
		Element M = Y.pow(d).mul(Z.pow(r));
		BigInteger temphash = Tools.hash_pp(d);

		BigInteger alpha = Tools.randomZp(p);
		BigInteger beta = Tools.randomZp(p);
		Element delta = Y.pow(alpha).mul(Z.pow(beta));
		BigInteger eta = Tools.hash_p((Y.toString() + Z.toString() + M.toString() + delta.toString()).getBytes());
		BigInteger alpha_hat = d.multiply(eta).add(alpha);
		BigInteger beta_hat = r.multiply(eta).add(beta);

		// PSP verify
		if (!Omega.contains(temphash.toString())) {
			eta = Tools.hash_p((Y.toString() + Z.toString() + M.toString() + delta.toString()).getBytes());
			Element left = M.pow(eta).mul(delta);
			Element right = Y.pow(alpha_hat).mul(Z.pow(beta_hat));
			if (left.isEqual(right)) {
				Omega.add(temphash.toString());
				v = Tools.randomZp(p);
				BigInteger temppow = (v.add(a).add(mu)).modInverse(p);
				W = (X.mul(M)).pow(temppow);
				System.out.println("Server Subscribe Success");
			}
		}
		// user verify
		Element left = pairing.pairing(W, A.mul(g.pow(v.add(mu))));
		Element right = pairing.pairing(X.mul(M), g);
		if (left.isEqual(right)) {
			System.out.println("User Subscribe Success");
		}

		/*
		 * ------------------------------------
		 */

//		// PLT register
//		BigInteger b = Tools.randomZp(p);
//		Element B = g.pow(b);
//
//		// PSP verify
//		Element R_ab = B.pow(a.modInverse(p));

		/*
		 * ------------------------------------
		 */
		for (int z = 0; z < nlist.length; z++) {
			int n = nlist[z];
			// user authenticate
			long usercosttime = 0;
			long servercosttime = 0;
			for (int k = 0; k < NUM; k++) {
				long t1 = System.currentTimeMillis();
				Element U = g.pow(d.add(mu).modInverse(p));

				BigInteger alpha_1 = Tools.randomZp(p);
				BigInteger alpha_2 = Tools.randomZp(p);
				BigInteger beta_1 = alpha_1.multiply(v);
				BigInteger beta_2 = alpha_2.multiply(v);
				Element W_1 = (Y.pow(alpha_1)).mul(Z.pow(alpha_2));
				Element W_2 = W.mul(Z.pow(alpha_1));
				BigInteger rho_v = Tools.randomZp(p);
				BigInteger rho_d = Tools.randomZp(p);
				BigInteger rho_r = Tools.randomZp(p);
				BigInteger rho_alpha_1 = Tools.randomZp(p);
				BigInteger rho_alpha_2 = Tools.randomZp(p);
				BigInteger rho_beta_1 = Tools.randomZp(p);
				BigInteger rho_beta_2 = Tools.randomZp(p);
				Element delta_1 = (Y.pow(rho_alpha_1)).mul(Z.pow(rho_alpha_2));
				BigInteger rho_v_negate = rho_v.negate().mod(p);
				Element delta_2 = W_1.pow(rho_v_negate).mul(Y.pow(rho_beta_1)).mul(Z.pow(rho_beta_2));
				Element delta_3 = U.pow(rho_d);
				Element delta_4 = pairing.pairing(W_2, g).pow(rho_v_negate).mul(eYg.pow(rho_d))
						.mul(eZA.pow(rho_alpha_1)).mul(eZgmu.pow(rho_alpha_1)).mul(eZg.pow(rho_r.add(rho_beta_1)));

				eta = Tools.hash_p(
						(X.toString() + Y.toString() + Z.toString() + W_1.toString() + W_2.toString() + U.toString()
								+ delta_1.toString() + delta_2.toString() + delta_3.toString() + delta_4.toString())
										.getBytes());

				BigInteger rho_v_hat = v.multiply(eta).add(rho_v);
				BigInteger rho_d_hat = d.multiply(eta).add(rho_d);
				BigInteger rho_r_hat = r.multiply(eta).add(rho_r);
				BigInteger rho_alpha_1_hat = alpha_1.multiply(eta).add(rho_alpha_1);
				BigInteger rho_alpha_2_hat = alpha_2.multiply(eta).add(rho_alpha_2);
				BigInteger rho_beta_1_hat = beta_1.multiply(eta).add(rho_beta_1);
				BigInteger rho_beta_2_hat = beta_2.multiply(eta).add(rho_beta_2);
				long t2 = System.currentTimeMillis();
				usercosttime = usercosttime + (t2 - t1);

				//calculate communication costs and storage costs
				commsize = delta_1.toBytes().length + delta_2.toBytes().length + delta_3.toBytes().length
						+ delta_4.toBytes().length + W_1.toBytes().length + W_1.toBytes().length + U.toBytes().length
						+ rho_v_hat.toByteArray().length + rho_d_hat.toByteArray().length
						+ rho_r_hat.toByteArray().length + rho_alpha_1_hat.toByteArray().length
						+ rho_alpha_2_hat.toByteArray().length + rho_beta_1_hat.toByteArray().length
						+ rho_beta_2_hat.toByteArray().length;
				long pairingparametersize = 0;
				try {
					pairingparametersize = Files.size(new File("pairingparams.properties").toPath());
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					System.out.println("The file cannot be found: pairingparams.properties");
					System.exit(-1);
				}
				long pubparamsize = g.toBytes().length + X.toBytes().length+Y.toBytes().length+Z.toBytes().length+A.toBytes().length+mu.toByteArray().length;
				userstore = pairingparametersize + pubparamsize +W.toBytes().length + v.toByteArray().length + d.toByteArray().length + r.toByteArray().length;
				serverstore = pairingparametersize + pubparamsize + a.toByteArray().length;
				for (int kk = 0; kk<n-1;kk++)
				{
					serverstore = serverstore + temphash.toByteArray().length+ U.toBytes().length;
				}
				
				// PSP verify
				long t3 = System.currentTimeMillis();
				eta = Tools.hash_p(
						(X.toString() + Y.toString() + Z.toString() + W_1.toString() + W_2.toString() + U.toString()
								+ delta_1.toString() + delta_2.toString() + delta_3.toString() + delta_4.toString())
										.getBytes());
				Element left_1 = W_1.pow(eta).mul(delta_1);
				Element right_1 = Y.pow(rho_alpha_1_hat).mul(Z.pow(rho_alpha_2_hat));
				Element left_2 = pairing.getG1().newOneElement().pow(eta).mul(delta_2);
				BigInteger rho_v_hat_negate = rho_v_hat.negate().mod(p);
				Element right_2 = W_1.pow(rho_v_hat_negate).mul(Y.pow(rho_beta_1_hat)).mul(Z.pow(rho_beta_2_hat));
				BigInteger mu_negate = mu.negate().mod(p);
				Element left_3 = (g.mul(U.pow(mu_negate))).pow(eta).mul(delta_3);
				Element right_3 = U.pow(rho_d_hat);
				Element left_4 = (pairing.pairing(W_2, A.mul(g.pow(mu))).mul(eXg.invert())).pow(eta).mul(delta_4);
				Element right_4 = pairing.pairing(W_2, g).pow(rho_v_hat_negate).mul(eYg.pow(rho_d_hat))
						.mul(eZA.pow(rho_alpha_1_hat)).mul(eZgmu.pow(rho_alpha_1_hat))
						.mul(eZg.pow(rho_r_hat.add(rho_beta_1_hat)));

				if (left_1.isEqual(right_1) & left_2.isEqual(right_2) & left_3.isEqual(right_3)
						& left_4.isEqual(right_4)) {
					// System.out.println("User Authentication Success");
					result = true;
				}
				long t4 = System.currentTimeMillis();
				System.out.print((t4 - t3) + " ");
				servercosttime = servercosttime + (t4 - t3);
			}
			System.out.println("Scheme, reserved number:" + n);
//			System.out.println("userstore:" + userstore);
//			System.out.println("serverstore:" + serverstore);
//			System.out.println("commsize:" + commsize);
//			System.out.println("user:" + usercosttime / NUM + "ms");
//			System.out.println("Server:" + servercosttime / NUM + "ms");
//			System.out.println(result);
		}
	}
}
