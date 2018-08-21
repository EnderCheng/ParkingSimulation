package ca.uwaterloo.cheng.tradition;

import java.io.File;
import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Files;
import java.security.SecureRandom;

import ca.uwaterloo.cheng.utils.Tools;
import it.unisa.dia.gas.jpbc.Element;
//import it.unisa.dia.gas.jpbc.Element;
import it.unisa.dia.gas.jpbc.Pairing;
import it.unisa.dia.gas.jpbc.PairingParameters;
import it.unisa.dia.gas.plaf.jpbc.pairing.PairingFactory;

public class BlacklistableAccum {
	private static Pairing pairing = PairingFactory.getPairing("pairingparams.properties");
	private static PairingParameters pp = PairingFactory.getPairingParameters("pairingparams.properties");

	public static void main(String[] args) {
		boolean result1 = false, result2 = false, result3 = false, result4 = false, result5 = false;
		int NUM = 100; // repeat times
		long commsize = 0;
		long userstore = 0;
		long serverstore = 0;

		BigInteger p = pp.getBigInteger("r");

		Element g_0 = pairing.getG1().newRandomElement().getImmutable();
		Element g_1 = pairing.getG1().newRandomElement().getImmutable();
		Element g_2 = pairing.getG1().newRandomElement().getImmutable();
		Element h_0 = pairing.getG2().newRandomElement().getImmutable();
		BigInteger gamma = Tools.randomZp(p);
		Element w = h_0.pow(gamma);
		BigInteger pprime = new BigInteger(1024, 80, new SecureRandom());
		BigInteger qprime = new BigInteger(1024, 80, new SecureRandom());
		BigInteger N = pprime.multiply(qprime);

		BigInteger gg = null;
		BigInteger hh = null;
		BigInteger order = (qprime.subtract(BigInteger.ONE)).multiply(pprime.subtract(BigInteger.ONE));
		while (true) {
			BigInteger temph = Tools.randomZStarN(N);
			BigInteger a = Tools.randomZStarN(N);
			BigInteger aModP = a.mod(pprime);
			BigInteger aModQ = a.mod(qprime);
			BigInteger hModP = temph.mod(pprime);
			BigInteger hModQ = temph.mod(qprime);
			if (aModP == BigInteger.ZERO || aModP == BigInteger.ONE || aModP == (pprime.subtract(BigInteger.ONE)))
				continue;
			if (aModQ == BigInteger.ZERO || aModQ == BigInteger.ONE || aModQ == (qprime.subtract(BigInteger.ONE)))
				continue;
			if (hModP == BigInteger.ZERO || hModP == BigInteger.ONE || hModP == (pprime.subtract(BigInteger.ONE)))
				continue;
			if (hModQ == BigInteger.ZERO || hModQ == BigInteger.ONE || hModQ == (qprime.subtract(BigInteger.ONE)))
				continue;
			gg = a.pow(2).mod(N);
			hh = temph.pow(2).mod(N);
			break;
		}
		BigInteger a = null;
		BigInteger b = null;
		BigInteger mupirme = null;
		BigInteger x = null;

		// Element left = pairing.pairing(A, w.mul(h_0.pow(e))).getImmutable();
		// Element right = pairing.pairing(g_0.mul(g_1.pow(x)).mul(g_2.pow(rpri)),
		// h_0).getImmutable();
		// System.out.println(left.isEqual(right));
		int l = 160;
		int[] nlist = { 10, 20, 50, 1000 }; // reserved user number
		// user authentication
		for (int zzz = 0; zzz < nlist.length; zzz++) {
			int n = nlist[zzz];
			long usercosttime = 0;
			long servercosttime = 0;

			// generate the accumulator and signature 
			BigInteger[] muall = new BigInteger[n];
			BigInteger rpri = null;
			BigInteger e = null;
			Element A = null;
			while (true) {
				BigInteger mu = BigInteger.ONE;
				for (int i = 0; i < n; i++) {
					muall[i] = new BigInteger(l, 120, new SecureRandom());
					mu = mu.multiply(muall[i]);
				}
				x = new BigInteger(l, 120, new SecureRandom()).mod(order).mod(p);
				mupirme = mu.mod(order);
				BigInteger vals[] = Tools.gcd(mupirme, x, order);
				a = vals[1];
				b = vals[2];
				if (a.multiply(mupirme).add(b.multiply(x)).mod(order).compareTo(BigInteger.ONE) == 0) {
					rpri = Tools.randomZp(p);
					Element c1 = g_1.pow(x).mul(g_2.pow(rpri));
					e = Tools.randomZp(p);
					A = (g_0.mul(c1)).pow(e.add(gamma).modInverse(p)).getImmutable();
					break;
				}
			}

			for (int kk = 0; kk < NUM; kk++) {
				long t1 = System.currentTimeMillis();
				BigInteger mu = BigInteger.ONE;
				for (int i = 0; i < n; i++) {
					mu = mu.multiply(muall[i]);
				}
				mupirme = mu.mod(order);
				BigInteger vals[] = Tools.gcd(mupirme, x, order);
				a = vals[1];
				b = vals[2];
				if (a.multiply(mupirme).add(b.multiply(x)).mod(order).compareTo(BigInteger.ONE) != 0) {
					System.out.println("Error: not co-prime");
					System.exit(-1);
				}
				// System.out.println(a.multiply(mupirme).add(b.multiply(x)).mod(order));
				int alength = a.bitLength();
				if (alength > N.bitLength() / 2) {
					while (true) {
						BigInteger temp = new BigInteger(l - 1, new SecureRandom());
						BigInteger k = (a.subtract(temp)).divide(x);
						BigInteger aprime = a.subtract(k.multiply(x));
						int tempalength = aprime.bitLength();
						// System.out.println(tempalength);
						if (aprime.compareTo(BigInteger.ONE) == 1 && tempalength < N.bitLength() / 2) {
							a = aprime;
							b = b.add(k.multiply(mupirme));
							break;
						}
					}
				}
				BigInteger d = gg.modPow(b.negate().mod(order), N);
				BigInteger c = gg.modPow(mupirme, N);
				// BigInteger left1 = c.modPow(a, N).mod(N);
				// BigInteger right1 = d.modPow(x, N).multiply(gg).mod(N);
				// if (left1.equals(right1))
				// System.out.println("success");
				BigInteger ww = new BigInteger(l, new SecureRandom());
				BigInteger rx = new BigInteger(l, new SecureRandom());
				BigInteger ra = new BigInteger(l, new SecureRandom());
				BigInteger rw = new BigInteger(l, new SecureRandom());
				BigInteger rz = new BigInteger(l, new SecureRandom());
				BigInteger re = new BigInteger(l, new SecureRandom());
				BigInteger cx = gg.modPow(x, N).multiply(hh.modPow(rx, N)).mod(N);
				BigInteger ca = gg.modPow(a, N).multiply(hh.modPow(ra, N)).mod(N);
				BigInteger cd = d.multiply(gg.modPow(ww, N)).mod(N);
				BigInteger cw = gg.modPow(ww, N).multiply(hh.modPow(rw, N)).mod(N);
				BigInteger z = x.multiply(ww).mod(order);
				BigInteger cz = gg.modPow(z, N).multiply(hh.modPow(rz, N)).mod(N);
				BigInteger ce = cd.modPow(x, N).multiply(hh.modPow(re, N)).mod(N);
				// if(ce.equals(d.modPow(x, N)
				// .multiply(gg.modPow(z, N))
				// .multiply(hh.modPow(re, N)).mod(N)))
				// {
				// System.out.println("success");
				// }
				// if (ce.equals(gg.modInverse(N).multiply(c.modPow(a, N)).multiply(gg.modPow(z,
				// N)).multiply(hh.modPow(re, N))
				// .mod(N))) {
				// System.out.println("success");
				// }

				BigInteger xhat = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger rxhat = new BigInteger(l, new SecureRandom()).mod(order);
				Element M = g_1.pow(x).mul(g_2.pow(rx)).getImmutable();
				Element delta_11 = g_1.pow(xhat).mul(g_2.pow(rxhat)).getImmutable();
				BigInteger delta_12 = gg.modPow(xhat, N).multiply(hh.modPow(rxhat, N)).mod(N);
				BigInteger eta_1 = Tools.hash_p((delta_11.toString() + delta_12.toString() + g_1.toString()
						+ g_2.toString() + gg.toString() + hh.toString()).getBytes());
				BigInteger sxhat = (x.multiply(eta_1).add(xhat)).mod(order);
				BigInteger srxhat = (rx.multiply(eta_1).add(rxhat)).mod(order);

				BigInteger xhat2 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger rxhat2 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger rehat2 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger delta_21 = cd.modPow(xhat2, N).multiply(hh.modPow(rehat2, N)).mod(N);
				BigInteger delta_22 = gg.modPow(xhat2, N).multiply(hh.modPow(rxhat2, N)).mod(N);
				BigInteger eta_2 = Tools.hash_p((cx.toString() + ce.toString() + delta_21.toString()
						+ delta_22.toString() + cd.toString() + gg.toString() + hh.toString()).getBytes());
				BigInteger srehat2 = (re.multiply(eta_2).add(rehat2)).mod(order);
				BigInteger sxhat2 = (x.multiply(eta_2).add(xhat2)).mod(order);
				BigInteger srxhat2 = (rx.multiply(eta_2).add(rxhat2)).mod(order);

				BigInteger xhat3 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger ahat3 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger zhat3 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger rahat3 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger rzhat3 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger rehat3 = new BigInteger(l, new SecureRandom()).mod(order);

				BigInteger delta_31 = c.modPow(ahat3, N).multiply(gg.modPow(zhat3, N)).multiply(hh.modPow(rehat3, N))
						.mod(N);
				BigInteger delta_32 = gg.modPow(ahat3, N).multiply(hh.modPow(rahat3, N)).mod(N);
				BigInteger delta_33 = gg.modPow(zhat3, N).multiply(hh.modPow(rzhat3, N)).mod(N);
				BigInteger delta_34 = cd.modPow(xhat3, N).multiply(hh.modPow(rehat3, N)).mod(N);
				;
				BigInteger eta_3 = Tools.hash_p((delta_31.toString() + delta_32.toString() + delta_33.toString()
						+ delta_34.toString() + ce.toString() + gg.toString() + c.toString() + hh.toString()
						+ ca.toString() + cz.toString() + cd.toString()).getBytes());

				BigInteger sxhat3 = (x.multiply(eta_3).add(xhat3)).mod(order);
				BigInteger sahat3 = (a.multiply(eta_3).add(ahat3)).mod(order);
				BigInteger szhat3 = (z.multiply(eta_3).add(zhat3)).mod(order);

				BigInteger srehat3 = (re.multiply(eta_3).add(rehat3)).mod(order);
				BigInteger srahat3 = (ra.multiply(eta_3).add(rahat3)).mod(order);
				BigInteger srzhat3 = (rz.multiply(eta_3).add(rzhat3)).mod(order);

				BigInteger rho = rz.subtract(x.multiply(rw)).mod(order);
				BigInteger zhat4 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger what4 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger xhat4 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger rzhat4 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger rwhat4 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger rxhat4 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger rhohat4 = new BigInteger(l, new SecureRandom()).mod(order);

				BigInteger delta_41 = gg.modPow(zhat4, N).multiply(hh.modPow(rzhat4, N)).mod(N);
				BigInteger delta_42 = gg.modPow(what4, N).multiply(hh.modPow(rwhat4, N)).mod(N);
				BigInteger delta_43 = gg.modPow(xhat4, N).multiply(hh.modPow(rxhat4, N)).mod(N);
				BigInteger delta_44 = cw.modPow(xhat4, N).multiply(hh.modPow(rhohat4, N)).mod(N);
				BigInteger eta_4 = Tools
						.hash_p((delta_41.toString() + delta_42.toString() + delta_43.toString() + delta_44.toString()
								+ cz.toString() + gg.toString() + cw.toString() + hh.toString() + cx.toString())
										.getBytes());

				BigInteger szhat4 = (z.multiply(eta_4).add(zhat4)).mod(order);
				BigInteger swhat4 = (ww.multiply(eta_4).add(what4)).mod(order);
				BigInteger sxhat4 = (x.multiply(eta_4).add(xhat4)).mod(order);

				BigInteger srzhat4 = (rz.multiply(eta_4).add(rzhat4)).mod(order);
				BigInteger srwhat4 = (rw.multiply(eta_4).add(rwhat4)).mod(order);
				BigInteger srxhat4 = (rx.multiply(eta_4).add(rxhat4)).mod(order);
				BigInteger srhohat4 = (rho.multiply(eta_4).add(rhohat4)).mod(order);

				BigInteger ahat5 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger xhat5 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger rahat5 = new BigInteger(l, new SecureRandom()).mod(order);
				BigInteger rxhat5 = new BigInteger(l, new SecureRandom()).mod(order);

				BigInteger delta_51 = gg.modPow(ahat5, N).multiply(hh.modPow(rahat5, N)).mod(N);
				BigInteger delta_52 = gg.modPow(xhat5, N).multiply(hh.modPow(rxhat5, N)).mod(N);

				// CL signature
				BigInteger r1 = Tools.randomZp(p);
				BigInteger r2 = Tools.randomZp(p);
				BigInteger sigma1 = r1.multiply(e).mod(p);
				BigInteger sigma2 = r2.multiply(e).mod(p);
				Element A1 = g_1.pow(r1).mul(g_2.pow(r2)).getImmutable();
				Element A2 = A.mul(g_2.pow(r1)).getImmutable();

				BigInteger r1hat = Tools.randomZp(p);
				BigInteger r2hat = Tools.randomZp(p);
				BigInteger ehat = Tools.randomZN(p);
				BigInteger xhat55 = Tools.randomZp(p);
				BigInteger rprihat = Tools.randomZp(p);
				BigInteger sigma1hat = Tools.randomZp(p);
				BigInteger sigma2hat = Tools.randomZp(p);

				Element delta_53 = g_1.pow(r1hat).mul(g_2.pow(r2hat)).getImmutable();
				Element delta_54 = A1.pow(ehat.negate().mod(p)).mul(g_1.pow(sigma1hat)).mul(g_2.pow(sigma2hat))
						.getImmutable();
				Element delta_55 = pairing.pairing(A2, h_0).pow(ehat.negate().mod(p))
						.mul(pairing.pairing(g_1, h_0).pow(xhat55)).mul(pairing.pairing(g_2, h_0).pow(sigma1hat))
						.mul(pairing.pairing(g_2, w).pow(r1hat)).mul(pairing.pairing(g_2, h_0).pow(rprihat))
						.getImmutable();

				BigInteger eta_5 = Tools
						.hash_p((delta_51.toString() + delta_52.toString() + delta_53.toString() + delta_54.toString()
								+ delta_55.toString() + gg.toString() + cx.toString() + hh.toString() + ca.toString())
										.getBytes());
				BigInteger sahat5 = (a.multiply(eta_5).add(ahat5)).mod(order);
				BigInteger sxhat5 = (x.multiply(eta_5).add(xhat5)).mod(order);
				BigInteger srahat5 = (ra.multiply(eta_5).add(rahat5)).mod(order);
				BigInteger srxhat5 = (rx.multiply(eta_5).add(rxhat5)).mod(order);

				BigInteger sr1hat = r1hat.subtract(eta_5.multiply(r1)).mod(p);
				BigInteger sr2hat = r2hat.subtract(eta_5.multiply(r2)).mod(p);
				BigInteger sehat = ehat.subtract(eta_5.multiply(e)).mod(p);
				BigInteger sxhat55 = xhat55.subtract(eta_5.multiply(x)).mod(p);
				BigInteger srprihat = rprihat.subtract(eta_5.multiply(rpri)).mod(p);
				BigInteger ssigma1hat = sigma1hat.subtract(eta_5.multiply(sigma1)).mod(p);
				BigInteger ssigma2hat = sigma2hat.subtract(eta_5.multiply(sigma2)).mod(p);
				long t2 = System.currentTimeMillis();
				usercosttime = usercosttime + (t2 - t1);
				commsize = 0;
				for(int i=0;i<n;i++)
				{
					commsize = commsize + muall[i].toByteArray().length;
				}
				commsize = commsize + delta_11.toBytes().length+delta_12.toByteArray().length+delta_21.toByteArray().length+
						delta_22.toByteArray().length+delta_31.toByteArray().length+delta_32.toByteArray().length+
						delta_33.toByteArray().length+delta_34.toByteArray().length+delta_41.toByteArray().length+
						delta_42.toByteArray().length+delta_43.toByteArray().length+delta_44.toByteArray().length+
						delta_51.toByteArray().length+delta_52.toByteArray().length+delta_53.toBytes().length+
						delta_54.toBytes().length+delta_55.toBytes().length+
						eta_1.toByteArray().length+eta_2.toByteArray().length+eta_3.toByteArray().length+
						eta_4.toByteArray().length+eta_5.toByteArray().length;
				commsize = commsize + cx.toByteArray().length+ca.toByteArray().length+cd.toByteArray().length+
						cw.toByteArray().length+ce.toByteArray().length+cz.toByteArray().length;
				commsize = commsize + M.toBytes().length+sxhat.toByteArray().length+srxhat.toByteArray().length+
						srehat2.toByteArray().length+srxhat2.toByteArray().length+sxhat2.toByteArray().length+
						sxhat3.toByteArray().length+sahat3.toByteArray().length+szhat3.toByteArray().length+
						srehat3.toByteArray().length+srahat3.toByteArray().length+srzhat3.toByteArray().length+
						szhat4.toByteArray().length+swhat4.toByteArray().length+sxhat4.toByteArray().length+
						srzhat4.toByteArray().length+srwhat4.toByteArray().length+srxhat4.toByteArray().length+
						srhohat4.toByteArray().length+sahat5.toByteArray().length+sxhat5.toByteArray().length+
						srahat5.toByteArray().length+srxhat5.toByteArray().length+sxhat55.toByteArray().length+
						sr1hat.toByteArray().length+sr2hat.toByteArray().length+sehat.toByteArray().length+
						srprihat.toByteArray().length+ssigma1hat.toByteArray().length+ssigma2hat.toByteArray().length;
				
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
				userstore = pubparamsize+ +pairingparametersize+x.toByteArray().length+rpri.toByteArray().length+A.toBytes().length+e.toByteArray().length+a.toByteArray().length+d.toByteArray().length+gg.toByteArray().length+hh.toByteArray().length;
				serverstore = pubparamsize+pairingparametersize+gamma.toByteArray().length+gg.toByteArray().length+hh.toByteArray().length;
				for(int kkk = 0; kkk < n ; kkk++)
				{
					serverstore = serverstore + muall[kkk].toByteArray().length;
				}
 
				// server
				long t3 = System.currentTimeMillis();
				eta_1 = Tools.hash_p((delta_11.toString() + delta_12.toString() + g_1.toString() + g_2.toString()
						+ gg.toString() + hh.toString()).getBytes());
				eta_2 = Tools.hash_p((cx.toString() + ce.toString() + delta_21.toString() + delta_22.toString()
						+ cd.toString() + gg.toString() + hh.toString()).getBytes());
				eta_3 = Tools.hash_p((delta_31.toString() + delta_32.toString() + delta_33.toString()
						+ delta_34.toString() + ce.toString() + gg.toString() + c.toString() + hh.toString()
						+ ca.toString() + cz.toString() + cd.toString()).getBytes());
				eta_4 = Tools
						.hash_p((delta_41.toString() + delta_42.toString() + delta_43.toString() + delta_44.toString()
								+ cz.toString() + gg.toString() + cw.toString() + hh.toString() + cx.toString())
										.getBytes());
				eta_5 = Tools
						.hash_p((delta_51.toString() + delta_52.toString() + delta_53.toString() + delta_54.toString()
								+ delta_55.toString() + gg.toString() + cx.toString() + hh.toString() + ca.toString())
										.getBytes());

				if (M.pow(eta_1).mul(delta_11).isEqual(g_1.pow(sxhat).mul(g_2.pow(srxhat)))
						&& ((cx.modPow(eta_1, N).multiply(delta_12)).mod(N))
								.equals((gg.modPow(sxhat, N).multiply(hh.modPow(srxhat, N))).mod(N))) {
					// System.out.println("verify success");
					result1 = true;
				}

				if (((ce.modPow(eta_2, N).multiply(delta_21)).mod(N))
						.equals((cd.modPow(sxhat2, N).multiply(hh.modPow(srehat2, N))).mod(N))
						&& ((cx.modPow(eta_2, N).multiply(delta_22)).mod(N))
								.equals((gg.modPow(sxhat2, N).multiply(hh.modPow(srxhat2, N))).mod(N))) {
					// System.out.println("verify success2");
					result2 = true;
				}
				if ((ce.multiply(gg)).modPow(eta_3, N).multiply(delta_31).mod(N).equals(
						c.modPow(sahat3, N).multiply(gg.modPow(szhat3, N)).multiply(hh.modPow(srehat3, N)).mod(N))
						&& ca.modPow(eta_3, N).multiply(delta_32).mod(N)
								.equals(gg.modPow(sahat3, N).multiply(hh.modPow(srahat3, N)).mod(N))
						&& cz.modPow(eta_3, N).multiply(delta_33).mod(N)
								.equals(gg.modPow(szhat3, N).multiply(hh.modPow(srzhat3, N)).mod(N))
						&& ce.modPow(eta_3, N).multiply(delta_34).mod(N)
								.equals(cd.modPow(sxhat3, N).multiply(hh.modPow(srehat3, N)).mod(N))) {
					// System.out.println("verify success3");
					result3 = true;
				}
				if (cz.modPow(eta_4, N).multiply(delta_41).mod(N)
						.equals(gg.modPow(szhat4, N).multiply(hh.modPow(srzhat4, N)).mod(N))
						&& cw.modPow(eta_4, N).multiply(delta_42).mod(N)
								.equals(gg.modPow(swhat4, N).multiply(hh.modPow(srwhat4, N)).mod(N))
						&& cx.modPow(eta_4, N).multiply(delta_43).mod(N)
								.equals(gg.modPow(sxhat4, N).multiply(hh.modPow(srxhat4, N)).mod(N))
						&& cz.modPow(eta_4, N).multiply(delta_44).mod(N)
								.equals(cw.modPow(sxhat4, N).multiply(hh.modPow(srhohat4, N)).mod(N))) {
					// System.out.println("verify success4");
					result4 = true;
				}

				Element Tpri1 = g_1.pow(sr1hat).mul(g_2.pow(sr2hat)).mul(A1.pow(eta_5)).getImmutable();
				Element Tpri2 = A1.pow(sehat.negate().mod(p)).mul(g_1.pow(ssigma1hat)).mul(g_2.pow(ssigma2hat))
						.getImmutable();
				Element temp = pairing.pairing(A2, w).mul(pairing.pairing(g_0, h_0).invert()).getImmutable();
				Element Tpri3 = pairing.pairing(A2, h_0).pow(sehat.negate().mod(p))
						.mul(pairing.pairing(g_1, h_0).pow(sxhat55)).mul(pairing.pairing(g_2, h_0).pow(ssigma1hat))
						.mul(pairing.pairing(g_2, w).pow(sr1hat)).mul(pairing.pairing(g_2, h_0).pow(srprihat))
						.mul(temp.pow(eta_5));

				BigInteger neweta_5 = Tools
						.hash_p((delta_51.toString() + delta_52.toString() + Tpri1.toString() + Tpri2.toString()
								+ Tpri3.toString() + gg.toString() + cx.toString() + hh.toString() + ca.toString())
										.getBytes());

				if (ca.modPow(eta_5, N).multiply(delta_51).mod(N)
						.equals(gg.modPow(sahat5, N).multiply(hh.modPow(srahat5, N)).mod(N))
						&& cx.modPow(eta_5, N).multiply(delta_52).mod(N)
								.equals(gg.modPow(sxhat5, N).multiply(hh.modPow(srxhat5, N)).mod(N))
						&& neweta_5.equals(eta_5)) {
					// System.out.println("verify success5");
					result5 = true;
				}
				long t4 = System.currentTimeMillis();
				servercosttime = servercosttime + (t4 - t3);
			}
			System.out.println("BLACAccum, the reserved number:" + n);
			System.out.println("userstore:" + userstore);
			System.out.println("serverstore:" + serverstore);
			System.out.println("commsize:" + commsize);
			System.out.println("user:" + usercosttime / NUM + "ms");
			System.out.println("Server:" + servercosttime / NUM + "ms");
			System.out.println(result1 && result2 && result3 && result4 & result5);
		}
	}
}
