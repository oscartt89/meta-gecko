import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.File;
import java.util.Random;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.*;
import org.biojava.bio.seq.*;

/**
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @license all bugs reserved (CopyLeft)
 */
public class Generador {
	// FIELDS
	private static final char[] _nucleotides = {'A','C','G','T'};
	private static final int _numOfNucleotides = _nucleotides.length;
	private static final String _fakeHeader = "> Fake header for metagenome";
	
	// METHODS
	/**
	 * Method used to generate a metagenome file with a number of reads with length specified.
	 * Sequence are randomly generated and read headers are all the same (a fake header). 
	 * @param file where dictionary will be written.
	 * @param numReads number of reads that will be written.
	 * @param rLength length of the reads.
	 */
	public static void generateMetagenomes(String file,int numReads, int rLength){
		if(numReads <= 0) throw new IllegalArgumentException("Num reads can't be zero or a negative number.");
		if(rLength <= 0) throw new IllegalArgumentException("Length of reads can't be zero or a negative number.");
		
		try{
			PrintWriter writer = new PrintWriter(new File(file));
				
			// Write reads
			for(int i=0; i<numReads; ++i){
				// Write fake header
				writer.println(_fakeHeader);
				// Write random nucleotides
				writer.println(randomSeq(rLength));
			}			
			writer.close();
		}catch(IOException ioe){
			ioe.printStackTrace();
		}		
	}
	
	
	/**
	 * This method is used to generate a random sequence of nucleotides.
	 * @param length of the sequence.
	 * @return a string with the sequence.
	 */
	public static String randomSeq(int length){
		String seq = new String();
		Random rdm = new Random();
		
		for(int j=0; j<length; ++j){
			seq += _nucleotides[rdm.nextInt(_numOfNucleotides)]; // random nucleotide
		}
		return seq;
	}
	
	
	/**
	 * This function is used to generate a set of reads that are subsequences of a genome file given.
	 * @param out file prineter where reads will be written.
	 * @param genome file, must be in FASTA format.
	 * @param matches number of subsequences to take.
	 * @param minL minimum length of read. If is impossible to generate a large enough read an exception will be thrown.
	 */
	public static void coincidentReads(PrintWriter out,String genome,int matches,int minL){
		try{
			// Open genome file stream
			BufferedReader br = new BufferedReader(new FileReader(new File(genome)));
			SequenceIterator seqIt = (SequenceIterator)SeqIOTools.fileToBiojava(SeqIOConstants.FASTA_DNA, br); 
			int aux = 20;
			
			// Start to read
			Sequence seq = seqIt.nextSequence();
						
			// Check for exceptions
			if(seq.length() < minL){
				throw new IllegalArgumentException("Full genome sequence isn't larger that minimum length specified.");
			}
			else if(seq.length() == minL){
				for(int i=0; i<matches; ++i)
					out.write("> FULL "+genome+"\n"+seq+"\n");
			}else{
				// Sequence stored on line variable.
				// Now take subsequences as reads
				Random rdm = new Random();
				int rdm_index, diff;
				for(int i=0; i<matches; ++i){
					rdm_index = rdm.nextInt(seq.length()); // Random index
					diff = (rdm_index + minL) >= seq.length()? Math.abs((rdm_index + minL)-seq.length()) : 0;
					out.write("> Subsequence of "+genome);
					if(diff != 0)
						out.write(" ("+(rdm_index-diff)+",)\n"+seq.subStr(rdm_index - diff,seq.length()-1).toUpperCase()+"\n");
					else
						out.write(" ("+rdm_index+","+(rdm_index+minL)+")\n"+seq.subStr(rdm_index, rdm_index+minL).toUpperCase()+"\n");
				}
			}
		}catch(IOException ioe){
			ioe.printStackTrace();
		}catch(BioException be){
			be.printStackTrace();
		}
	}
	
	
	// MAIN
	public static void main(String[] args){
// USE THIS IF YOU WANT GENERATE A SET OF RANDOM METAGENOMES IN A RANGE OF LENGTH AND NUM_READS
//		if(args.length != 7){
//			System.err.println("Bad call error. Use: Executable out minReads maxReads minL maxL Rdif Ldif\nNOTE: dont't write extension of file in outFileBase.");
//		}else{
//			int minR = Integer.valueOf(args[1]), maxR = Integer.valueOf(args[2]);
//			int minL = Integer.valueOf(args[3]), maxL = Integer.valueOf(args[4]);
//			int difR = Integer.valueOf(args[5]), difL = Integer.valueOf(args[6]);
//			
//			if(difR <= 0 | difL <= 0) throw new IllegalArgumentException("Increment can't be zero or a negative number.");
//			if(maxR < minR) throw new IllegalArgumentException("MaxR can't be less than minR.");
//			if(maxL < minL) throw new IllegalArgumentException("MaxL can't be less than minL.");
//			
//			
//			// Generate metagenomes
//			for(int i=minR; i<=maxR; i+=difR)
//				for(int j=minL; j<=maxL;j+=difL)
//					generateMetagenomes(args[0]+"_R"+i+"L"+j+".fasta", i, j);
//			// If an exception is thrown stop the execution
//		}
	
// USE THIS IF YOU WANT GENERATE A METAGENOME FROM A GENOME GIVEN
		if(args.length != 4){
			System.err.println("Bad call error. Use: Executable out genome numReads minLength.");
		}else{
			try{
				PrintWriter writer = new PrintWriter(new File(args[0]));
				coincidentReads(writer,args[1], Integer.valueOf(args[2]), Integer.valueOf(args[3]));
				writer.close();
			}catch(IOException ioe){
				ioe.printStackTrace();
			}
		}
	}// END MAIN
}// END CLASS
