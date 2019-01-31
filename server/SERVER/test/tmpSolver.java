import java.io.*;

/**
 * Simulating the real solver.
 */ 
public class tmpSolver {

	/**
	 *
	 */
	public static void main (String [] args) { new tmpSolver(); }

	public tmpSolver() {
		long LL = 0L;
		for(int gg=0; gg<100000; gg++) 
			for(int kk=0; kk < 20000; kk++) 
				if ( LL == 987654 ) 
					LL +=2033L;
		System.out.println("LL is: " + LL);
	}
}
