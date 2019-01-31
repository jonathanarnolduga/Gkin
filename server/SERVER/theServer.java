import java.util.*;
import java.net.*;
import java.io.*;

/**
 *
 */ 
public class theServer {

	/**
	 *
	 */
	public static void main (String [] args) {
		if( args.length != 2 ) {
			System.err.println("USAGE: provide an integer which is the port number we will bind to.");
			System.err.println("       The third parameter is the path to the SOLVER.");
			System.err.println("Example:");
			System.err.println("        java theServer 8888 \"tmpSolver\";  USE THIS FOR TESTING");
			System.err.println("        java theServer 8888 kin;            THIS IS THE REAL KINSOLVER");
			System.exit(-1);
		}
		new theServer(Integer.parseInt(args[0]), args[1]);		
	}

	/**
	 * Binds to the specified port and loops forever accepting connections.
	 * Each new connection is handled by a separate thread.
	 * @param portNum the port number we will bind to.
	 */
	public theServer(int portNum, String solverBin) {
		Random rand = new Random(System.nanoTime());
		String ff = new String(rand.nextInt(100000) + "_" + System.nanoTime()+ ".txt");		// TEMP FILE NAME
		try {
       			ServerSocket serverSocket = new ServerSocket(portNum);
			while( true ) {
       				Socket connectSock = serverSocket.accept();

				//theServerChild sc = new theServerChild(connectSock, ff, solverBin);
				theServerChild sc = new theServerChild(connectSock, "kin.i01", solverBin);	// HARDCODED NAME !!!!!
				sc.start();
			}
		}
     		catch (SocketException se) {
       			System.out.println("Can't create socket. " + se);
			System.exit(-1);
		}
     		catch (IOException ioe) {
       			System.out.println("IO Exception in Socket. " + ioe);
			System.exit(-1);
     		}
	}

}




class theServerChild extends Thread {
  	Socket connectSock;
	String tmp_inputFile = null;
	String solverBin = null;

	/**
	 * We will handle the new connection.
	 * @param connectSock is the connected socket to the client.
	 * @param thmp_inputFile is the file name of the incoming data, the client sends us data and we store it in this file.
	 */
	theServerChild(Socket connectSock, String tmp_inputFile, String solverBin) {
		this.connectSock = connectSock;
		this.tmp_inputFile = tmp_inputFile;
		this.solverBin = solverBin;
	}

	public void run() {
		BufferedReader in = null;
		PrintStream out = null;
		PrintWriter fw = null;
		BufferedReader br = null;
     		try {
       			in = new BufferedReader(new InputStreamReader(connectSock.getInputStream()));
       			out = new PrintStream(connectSock.getOutputStream());

			fw = new PrintWriter(new FileWriter(tmp_inputFile));

       			while (true) {
         			String aLine = in.readLine();
				if( aLine == null ) {
					// System.out.println("Other side terminated OR done sending us stuff.");
					// System.out.println("Don't think we are supposed to ever see this line.");
					break;
				}

         			if (aLine.equalsIgnoreCase("THE END")) {
					break;
         			}
				else {
					fw.println(aLine);		// store this to a file
				}
      			}
			if( fw != null ) {
				fw.close();
			}

			/*
			 * At this point we have received the GUI generated Input file to the Solver.
			 * What we need to do next:
			 * 	1. Run the Solver
			 * 	2. Open the file(s) that the solver generated.
			 * 	3. Send back the generated file
			 * 	4. Close the file(s)
			 *
			 */

			// HERE the server does calculations 
			/*
			 * 1. start a checker thread
			 * 2. Open the exec
			 * 3. wait till it is done
			 * 4. done
			 * 5. in the mean time if the checker receives a terminate signal from client
			 *    we should abort the solver
			 */	

			Checker ch = new Checker(this, connectSock, in);
			ch.start();






			try {
				//String [] cmdArgs = { "java", "tmpSolver" };
				String [] cmdArgs = { solverBin };

				Process process = null;
				process = Runtime.getRuntime().exec(cmdArgs, null, new File("."));
				try {
					int exitVal = process.waitFor();

					if( exitVal != 0 ) {
						System.out.println("SERVER error in exec (return code): " + exitVal);
						return;
					}
				}
				catch (InterruptedException e) {
					// user on client pressed the interrupt button
					// and the checker in this file throu an exception because of the failed read (test failed read).
					// We catch the exception and we kill the solver process.
					// 
					//
					//System.out.println("SERVER Interrupted: " + e);
					process.destroy();
					return;
				}
			}
			catch(IOException e1) {
				System.out.println("SERVER IO error: " + e1);
				return;
			}

			//System.out.println("DONE CALCULATIONS	DONE CALCULATIONS");







			/*
			 * and here sends results back
			 */
			br = new BufferedReader( new FileReader("kin.o02"));	// LEON for now use this hardcoded name



			// send the results
			String oneLine = null;
			while( (oneLine = br.readLine() ) != null) {
				out.println(oneLine);
			}

     		}
     		catch (IOException ioe) {
       			System.out.println("SERVER IO Exception. " + ioe);
			// return;
     		}
     		catch (Exception sc) {
       			System.out.println("SERVER LEON Exception. " + sc);
			// return;
     		}

		finally {
			try {
				File fin_tmp = new File(tmp_inputFile);

				out.println("THE END");		// indicate the end of the transaction

				fin_tmp.delete();
				if( in != null )  in.close();
				if( out != null ) out.close();
				connectSock.close();

				// LEON HERE WE SHOULD DELETE THE FILE NAME pointed to by br.
				if(br != null) br.close();
			}
			catch(IOException e) {
				System.out.println("in finally (server): " + e);
			}
			//System.out.println("DONE TRANSACTION");
		}
	}
}

class Checker extends Thread {
	Thread parentth = null;
	Socket s = null;
	BufferedReader in = null;
	Checker(Thread parentth, Socket s, BufferedReader in) {
		this.parentth = parentth;
		this.s = s;
		this.in = in;
	}

	public void run() {
		while( true ) {
			try {
				Thread.sleep(500);

				char []buf = new char[2];
				int h = in.read(buf, 0, 1);
				//System.out.println("return of read: " + h);
				parentth.interrupt();
				return;
			}
     			catch (IOException ioe) {
				//Possibly because the solver finished before we even do the read() above, which is fine
				//
				//
       				//System.out.println("LEON SERVER CHECKER IO Exception. " + ioe);
				parentth.interrupt();		// in this case I don't think this hurts. (don't think we need it though).
				return;
     			}
			catch(InterruptedException ie) {
				System.out.println("SERVER CHECKER INTERRUPTED");	// go back in the loop
			}
		}
	}

}
