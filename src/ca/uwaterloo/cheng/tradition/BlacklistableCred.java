package ca.uwaterloo.cheng.tradition;

import java.io.File;
import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Files;
import java.security.SecureRandom;

import ca.uwaterloo.cheng.utils.Tools;
import it.unisa.dia.gas.jpbc.Element;
import it.unisa.dia.gas.jpbc.Pairing;
import it.unisa.dia.gas.jpbc.PairingParameters;
import it.unisa.dia.gas.plaf.jpbc.pairing.PairingFactory;

public class BlacklistableCred {
	private static Pairing pairing = PairingFactory.getPairing("pairingparams.properties");
	private static PairingParameters pp = PairingFactory.getPairingParameters("pairingparams.properties");

	public static void main(String[] args) {
		int result = -1;
		int NUM = 1; //repeat times
		long commsize = 0;
		long userstore = 0;
		long serverstore = 0;

		// keygen
		BigInteger p = pp.getBigInteger("r");
		Element g_0 = pairing.getG1().newRandomElement().getImmutable();
		Element g_1 = pairing.getG1().newRandomElement().getImmutable();
		Element g_2 = pairing.getG1().newRandomElement().getImmutable();
		Element h_0 = pairing.getG2().newRandomElement().getImmutable();
		BigInteger gamma = Tools.randomZp(p);
		Element w = h_0.pow(gamma);

		// register
		BigInteger x = Tools.randomZp(p);
		BigInteger ypri = Tools.randomZp(p);
		Element C = g_1.pow(x).mul(g_2.pow(ypri)).getImmutable();
		// proof is omitted

		BigInteger e = Tools.randomZp(p);
		BigInteger ypp = Tools.randomZp(p);
		Element A = (g_0.mul(C).mul(g_2.pow(ypp))).pow(e.add(gamma).modInverse(p)).getImmutable();
		BigInteger y = ypri.add(ypp).mod(p);
		Element left = pairing.pairing(A, w.mul(h_0.pow(e))).getImmutable();
		Element right = pairing.pairing(g_0.mul(g_1.pow(x)).mul(g_2.pow(y)), h_0).getImmutable();
		System.out.println(left.isEqual(right));

		// generate blacklist
		int[] nlist = {10, 20, 50, 1000}; // reserved user number
		for (int z = 0; z < nlist.length; z++) {
			int n = nlist[z];
			String[] BL1 = new String[n];
			Element[] BL2 = new Element[n];

			for (int i = 0; i < n; i++) {
				BigInteger temp = BigInteger.ZERO;
				do {
					temp = Tools.randomZp(p);
				} while (temp.compareTo(x) == 0);
				BigInteger temp1 = new BigInteger(128, new SecureRandom());
				BL1[i] = temp1.toString();
				BL2[i] = Tools.hash_g((BL1[i] + "Bob").getBytes()).pow(temp).getImmutable();
			}
			// user authenticate
			long usercosttime = 0;
			long servercosttime = 0;
			for (int k = 0; k < NUM; k++) {
				long t1 = System.currentTimeMillis();
				Element[] b = new Element[n];
				for (int i = 0; i < n; i++) {
					b[i] = Tools.hash_g((BL1[i] + "Bob").getBytes()).getImmutable();
					Element t = b[i].pow(x).getImmutable();
					if (t.isEqual(BL2[i])) {
						System.out.println("In blacklist");
					}
				}

				BigInteger temp1 = new BigInteger(128, new SecureRandom());
				String s = temp1.toString();
				Element newb = Tools.hash_g((s + "Bob").getBytes()).getImmutable();
				Element newt = newb.pow(x).getImmutable();
				BigInteger rho_1 = Tools.randomZp(p);
				BigInteger rho_2 = Tools.randomZp(p);
				BigInteger rho_3 = Tools.randomZp(p);
				BigInteger rho_4 = Tools.randomZp(p);
				Element A_1 = g_1.pow(rho_1).mul(g_2.pow(rho_2)).getImmutable();
				Element A_2 = A.mul(g_2.pow(rho_1)).getImmutable();
				Element A_3 = g_1.pow(rho_3).mul(g_2.pow(rho_4)).getImmutable();
				Element[] Ahat = new Element[n];
				for (int i = 0; i < n; i++) {
					Ahat[i] = (b[i].pow(x).mul(BL2[i].invert())).pow(rho_3).getImmutable();
				}
				// BigInteger alpha_1 = rho_1.multiply(e);
				// BigInteger alpha_2 = rho_2.multiply(e);
				// BigInteger beta_3 = rho_3.multiply(x);
				// BigInteger beta_4 = rho_4.multiply(x);

				// Element e_0 = pairing.pairing(g_0, h_0).getImmutable();
				// Element e_1 = pairing.pairing(g_1, h_0).getImmutable();
				// Element e_2 = pairing.pairing(g_2, h_0).getImmutable();

				BigInteger re = Tools.randomZp(p);
				BigInteger rx = Tools.randomZp(p);
				BigInteger ry = Tools.randomZp(p);
				BigInteger rrho1 = Tools.randomZp(p);
				BigInteger rrho2 = Tools.randomZp(p);
				BigInteger rrho3 = Tools.randomZp(p);
				BigInteger rrho4 = Tools.randomZp(p);
				BigInteger ralpha1 = Tools.randomZp(p);
				BigInteger ralpha2 = Tools.randomZp(p);
				BigInteger rbeta3 = Tools.randomZp(p);
				BigInteger rbeta4 = Tools.randomZp(p);

				Element T_1 = g_1.pow(rrho1).mul(g_2.pow(rrho2)).getImmutable();
				Element T_2 = A_1.pow(re.negate().mod(p)).mul(g_1.pow(ralpha1)).mul(g_2.pow(ralpha2)).getImmutable();
				Element T_3 = g_1.pow(rrho3).mul(g_2.pow(rrho4)).getImmutable();
				Element T_4 = A_3.pow(rx.negate().mod(p)).mul(g_1.pow(rbeta3)).mul(g_2.pow(rbeta4)).getImmutable();
				Element T_5 = pairing.pairing(A_2, h_0).pow(re.negate().mod(p)).mul(pairing.pairing(g_1, h_0).pow(rx))
						.mul(pairing.pairing(g_2, h_0).pow(ry)).mul(pairing.pairing(g_2, w).pow(rrho1))
						.mul(pairing.pairing(g_2, h_0).pow(ralpha1)).getImmutable();
				Element[] That = new Element[n];
				for (int i = 0; i < n; i++) {
					That[i] = b[i].pow(rbeta3).mul(BL2[i].pow(rrho3.negate().mod(p))).getImmutable();
				}
				Element BigT = newb.pow(rbeta3).mul(newt.pow(rrho3.negate().mod(p))).getImmutable();
				String str = "";
				for (int i = 0; i < n; i++) {
					str = str + That[i].toString() + Ahat[i].toString();
				}

				BigInteger c = Tools
						.hash_p((A_1.toString() + A_2.toString() + A_3.toString() + T_1.toString() + T_2.toString()
								+ T_3.toString() + T_4.toString() + T_5.toString() + BigT.toString() + str).getBytes());

				BigInteger se = re.subtract(c.multiply(e)).mod(p);
				BigInteger sx = rx.subtract(c.multiply(x)).mod(p);
				BigInteger sy = ry.subtract(c.multiply(y)).mod(p);
				BigInteger srho1 = rrho1.subtract(c.multiply(rho_1)).mod(p);
				BigInteger srho2 = rrho2.subtract(c.multiply(rho_2)).mod(p);
				BigInteger srho3 = rrho3.subtract(c.multiply(rho_3)).mod(p);
				BigInteger srho4 = rrho4.subtract(c.multiply(rho_4)).mod(p);
				BigInteger salpha1 = ralpha1.subtract(c.multiply(rho_1).multiply(e)).mod(p);
				BigInteger salpha2 = ralpha2.subtract(c.multiply(rho_2).multiply(e)).mod(p);
				BigInteger sbeta3 = rbeta3.subtract(c.multiply(rho_3).multiply(x)).mod(p);
				BigInteger sbeta4 = rbeta4.subtract(c.multiply(rho_4).multiply(x)).mod(p);
				long t2 = System.currentTimeMillis();
				usercosttime = usercosttime + (t2 - t1);
				
				//calculate the communication overheads and storage costs
				commsize = A_1.toBytes().length+A_2.toBytes().length+A_3.toBytes().length;
				for (int cc=0; cc<n;cc++)
				{
					commsize = commsize + Ahat[cc].toBytes().length+BL1[cc].getBytes().length+BL2[cc].toBytes().length;
				}
				commsize = commsize + c.toByteArray().length+se.toByteArray().length
						+sx.toByteArray().length+sy.toByteArray().length
						+srho1.toByteArray().length+srho2.toByteArray().length
						+srho3.toByteArray().length+srho4.toByteArray().length
						+salpha1.toByteArray().length+salpha2.toByteArray().length
						+sbeta3.toByteArray().length+sbeta3.toByteArray().length;
				
				long pairingparametersize = 0;
				try {
					pairingparametersize = Files.size(new File("pairingparams.properties").toPath());
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
					System.out.println("The file cannot be found: pairingparams.properties");
					System.exit(-1);
				}
				long pubparamsize = h_0.toBytes().length + g_0.toBytes().length+g_1.toBytes().length+g_2.toBytes().length+w.toBytes().length;
				userstore = pubparamsize+ +pairingparametersize+x.toByteArray().length+y.toByteArray().length+A.toBytes().length+e.toByteArray().length;
				serverstore = pubparamsize+pairingparametersize+gamma.toByteArray().length;
				for(int kk = 0; kk < n ; kk++)
				{
					serverstore = serverstore + BL1[kk].getBytes().length+BL2[kk].toBytes().length;
				}
				// server verify
				long t3 = System.currentTimeMillis();
				Element Tpri1 = g_1.pow(srho1).mul(g_2.pow(srho2)).mul(A_1.pow(c)).getImmutable();
				Element Tpri2 = A_1.pow(se.negate().mod(p)).mul(g_1.pow(salpha1)).mul(g_2.pow(salpha2)).getImmutable();
				Element Tpri3 = g_1.pow(srho3).mul(g_2.pow(srho4)).mul(A_3.pow(c)).getImmutable();
				Element Tpri4 = A_3.pow(sx.negate().mod(p)).mul(g_1.pow(sbeta3)).mul(g_2.pow(sbeta4)).getImmutable();
				Element temp = pairing.pairing(A_2, w).mul(pairing.pairing(g_0, h_0).invert()).getImmutable();
				Element Tpri5 = pairing.pairing(A_2, h_0).pow(se.negate().mod(p)).mul(pairing.pairing(g_1, h_0).pow(sx))
						.mul(pairing.pairing(g_2, h_0).pow(sy)).mul(pairing.pairing(g_2, w).pow(srho1))
						.mul(pairing.pairing(g_2, h_0).pow(salpha1)).mul(temp.pow(c));
				Element[] Thatpri = new Element[n];
				for (int i = 0; i < n; i++) {
					Thatpri[i] = b[i].pow(sbeta3).mul(BL2[i].pow(srho3.negate().mod(p))).mul(Ahat[i].pow(c))
							.getImmutable();
				}

				Element BigTpri = newb.pow(sbeta3).mul(newt.pow(srho3.negate().mod(p))).getImmutable();

				String newstr = "";
				for (int i = 0; i < n; i++) {
					newstr = newstr + Thatpri[i].toString() + Ahat[i].toString();
				}

				BigInteger newc = Tools
						.hash_p((A_1.toString() + A_2.toString() + A_3.toString() + Tpri1.toString() + Tpri2.toString()
								+ Tpri3.toString() + Tpri4.toString() + Tpri5.toString() + BigTpri.toString() + newstr)
										.getBytes());
				result = newc.compareTo(c);
				long t4 = System.currentTimeMillis();
				servercosttime = servercosttime + (t4 - t3);
			}
			System.out.println("BLAC, the reserved number:"+ n);
			System.out.println("userstore:" + userstore);
			System.out.println("serverstore:" + serverstore);
			System.out.println("commsize:"+ commsize);
			System.out.println("user:" + usercosttime / NUM + "ms");
			System.out.println("Server:" + servercosttime / NUM + "ms");
			System.out.println(result);
		}

	}

}
