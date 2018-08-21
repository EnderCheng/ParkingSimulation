package ca.uwaterloo.cheng.tradition;

import java.io.File;
import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Files;

import ca.uwaterloo.cheng.utils.Tools;
import it.unisa.dia.gas.jpbc.Element;
import it.unisa.dia.gas.jpbc.Pairing;
import it.unisa.dia.gas.jpbc.PairingParameters;
import it.unisa.dia.gas.plaf.jpbc.pairing.PairingFactory;

public class LinkableSig {
	private static Pairing pairing = PairingFactory.getPairing("pairingparams.properties");
	private static PairingParameters pp = PairingFactory.getPairingParameters("pairingparams.properties");

	public static void main(String[] args) {
		BigInteger p = pp.getBigInteger("r");
		Element g = pairing.getG1().newRandomElement().getImmutable();
		Element h = pairing.getG1().newRandomElement().getImmutable();
		long commsize = 0;
		long userstore = 0;
		long serverstore_10 = 0;
		long serverstore_20 = 0;
		long serverstore_50 = 0;
		long serverstore_1000 = 0;
		int[] nlist = {100, 200, 500, 10000}; // total user number
		int NUM = 1; //repeat times
		int result = -1;
		// keygen
		for (int z = 0; z < nlist.length; z++) {
			int n = nlist[z];
			BigInteger[] x = new BigInteger[n];
			BigInteger[] y = new BigInteger[n];
			Element[] Z = new Element[n];
			for (int i = 0; i < n; i++) {
				x[i] = Tools.randomZp(p);
				y[i] = Tools.randomZp(p);
				Z[i] = g.pow(x[i]).mul(h.pow(y[i])).getImmutable();
			}

			// user sign
			long usercosttime = 0;
			long servercosttime = 0;
			for (int k = 0; k < NUM; k++) {
				long t1 = System.currentTimeMillis();
				String event = "Order";
				Element e = Tools.hash_g(event.getBytes()).getImmutable();
				Element t = e.pow(x[0]).getImmutable();
				BigInteger[] c = new BigInteger[n];
				for (int i = 0; i < n - 1; i++) {
					c[i] = Tools.randomZp(p);
				}
				BigInteger rx = Tools.randomZp(p);
				BigInteger ry = Tools.randomZp(p);
				Element K = g.pow(rx).mul(h.pow(ry)).getImmutable();
				for (int i = 0; i < n - 1; i++) {
					K = K.mul(Z[i + 1].pow(c[i])).getImmutable();
				}
				BigInteger call = BigInteger.ZERO;
				for (int i = 0; i < n - 1; i++) {
					call = call.add(c[i]);
				}
				Element Kpri = e.pow(rx).mul(t.pow(call)).getImmutable();
				String hashstring = "";
				for (int i = 0; i < n; i++) {
					hashstring = hashstring + Z[i].toString();
				}
				BigInteger hash = Tools.hash_p(
						(hashstring + event + t.toString() + "parking" + K.toString() + Kpri.toString()).getBytes());
				BigInteger cnew = hash.subtract(call).mod(p);
				// System.out.println(call.add(cnew).mod(p).compareTo(hash));
				BigInteger xhat = rx.subtract(cnew.multiply(x[0])).mod(p);
				BigInteger yhat = ry.subtract(cnew.multiply(y[0])).mod(p);
				long t2 = System.currentTimeMillis();
				usercosttime = usercosttime + (t2 - t1);
				
				//calculate communication overheads and storage costs
				commsize = 0;
				for (int cc = 0;  cc < n-1; cc++)
				{
					commsize = commsize + c[cc].toByteArray().length;
				}
				commsize = commsize + xhat.toByteArray().length + yhat.toByteArray().length+cnew.toByteArray().length+t.toBytes().length; 
				long pairingparametersize = 0;
				try {
					pairingparametersize = Files.size(new File("pairingparams.properties").toPath());
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
					System.out.println("The file cannot be found: pairingparams.properties");
					System.exit(-1);
				}
				long pubparamsize = h.toBytes().length + g.toBytes().length;
				userstore = pairingparametersize + pubparamsize+ x[0].toByteArray().length+y[0].toByteArray().length;
				serverstore_10 = pairingparametersize + pubparamsize;
				serverstore_20 = pairingparametersize + pubparamsize;
				serverstore_50 = pairingparametersize + pubparamsize;
				serverstore_1000 = pairingparametersize + pubparamsize;
				for(int kk = 0; kk<n;kk++)
				{
					userstore = userstore + Z[kk].toBytes().length;
					serverstore_10 = serverstore_10 + Z[kk].toBytes().length;
					serverstore_20 = serverstore_10 + Z[kk].toBytes().length;
					serverstore_50 = serverstore_10 + Z[kk].toBytes().length;
					serverstore_1000 = serverstore_10 + Z[kk].toBytes().length;
				}
				for(int kk = 0; kk<10;kk++)
				{
					serverstore_10 = serverstore_10 + t.toBytes().length;
				}
				for(int kk = 0; kk<20;kk++)
				{
					serverstore_20 = serverstore_20 + t.toBytes().length;
				}
				for(int kk = 0; kk<50;kk++)
				{
					serverstore_50 = serverstore_50 + t.toBytes().length;
				}
				
				for(int kk = 0; kk<1000;kk++)
				{
					serverstore_1000 = serverstore_1000 + t.toBytes().length;
				}
				
				// server verify
				long t3 = System.currentTimeMillis();
				Element temp1 = g.pow(xhat).mul(h.pow(yhat)).mul(Z[0].pow(cnew)).getImmutable();

				hashstring = "";
				for (int i = 0; i < n; i++) {
					hashstring = hashstring + Z[i].toString();
				}

				for (int i = 0; i < n - 1; i++) {
					temp1 = temp1.mul(Z[i + 1].pow(c[i])).getImmutable();
				}

				call = cnew;
				for (int i = 0; i < n - 1; i++) {
					call = call.add(c[i]);
				}

				Element temp2 = e.pow(xhat).mul(t.pow(call)).getImmutable();
				hash = Tools
						.hash_p((hashstring + event + t.toString() + "parking" + temp1.toString() + temp2.toString())
								.getBytes());
				result = call.mod(p).compareTo(hash);
				long t4 = System.currentTimeMillis();
				servercosttime = servercosttime + (t4 - t3);
			}
			System.out.println("LSR, total number:" + n);
			System.out.println("userstore:" + userstore);
			System.out.println("serverstore, 10 reserved number:" + serverstore_10);
			System.out.println("serverstore, 20 reserved number:" + serverstore_20);
			System.out.println("serverstore, 50 reserved number:" + serverstore_50);
			System.out.println("serverstore, 1000 reserved number:" + serverstore_1000);
			System.out.println("commsize:" + commsize);
			System.out.println("user:" + usercosttime / NUM + "ms");
			System.out.println("Server:" + servercosttime / NUM + "ms");
			System.out.println(result);
		}
	}
}
