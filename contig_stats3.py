### Boas Pucker ###
### b.pucker@tu-bs.de ###

## based on script contig_stats.py Pucker et al., 2016. doi:10.1371/journal.pone.0164321

import re, sys, os
from operator import itemgetter

# --- end of imports --- #

__version__ = "v3.0"

__citation__ = "Oliveira et al. (2025)"	#update

__usage__ = """ 
				Calculation of assembly statistics (""" + __version__ + """):
				python3 contig_stats3.py
				--in <INPUT_FASTA_FILENAME>
				
				optional:
				--min <MIN_CONTIG_LENGTH> [1000]
				--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
				--exp <EXPRESSION_FILE(normalized)>
				
				bug reports and feature requests: b.pucker@tu-bs.de
				Please cite: """ + __citation__ + """
			"""

def calculate_formal_contig_stats_expX( filename, expX ):
	"""! @brief calculates some formal stats of the given multiple fasta file (assembly)
	
		@param filename (string) full path to a assembly output file (multiple fasta file)
		
		@return (dictionary) contains all formal stats of the analyzed assembly
		
		@author Boas Pucker
	"""
	
	sys.stdout.write( "calculation of formal assembly stats ... please wait!\n" )
	sys.stdout.flush()
	
	number_of_bases_without_N = 0	#counts all bases without N
	number_of_gc = 0		#counts occurences of G or C in sequence
	contig_lengths = []		#lengths of all contigs in the assembly; used for calculation of min, max and mean
	exp_contig_lengths = []
	
	with open( filename, 'r' ) as f:
		header = f.readline().strip()[1:]
		if " " in header:
			header = header.split(' ')[0]
		line = f.readline()
		sequence = []
		counter = 1
		while line:
			if line[0] == '>':	#new header => evaluate current sequence and set back to empty string
				sequence = "".join( sequence )
				for base in sequence.upper():
					if base == 'G' or base == 'C':
						number_of_gc += 1
						number_of_bases_without_N += 1
					elif base == 'A' or base == 'T':
						number_of_bases_without_N += 1
				contig_lengths.append( len( sequence ) )
				try:
					expX[ header ]
					exp_contig_lengths.append( len( sequence ) )
				except KeyError:
					pass
				sequence = []
				header = line.strip()[1:]
				if " " in header:
					header = header.split(' ')[0]
			else:
				sequence.append( line.strip() )
			line = f.readline()
			counter += 1
			if counter % 1000 == 0:
				sys.stdout.write( str( counter/1000 ) + ' x1000 lines processed\n' )
				sys.stdout.flush()
		#place block from new header here again (for last sequence in file)
		sequence = "".join( sequence )
		for base in sequence.upper():
			if base == 'G' or base == 'C':
				number_of_gc += 1
				number_of_bases_without_N += 1
			elif base == 'A' or base == 'T':
				number_of_bases_without_N += 1
		contig_lengths.append( len( sequence ) )
		try:
			expX[ header ]
			exp_contig_lengths.append( len( sequence ) )
		except KeyError:
			pass
	
	# --- calculate remaining stats --- #
	number_of_contigs = len( contig_lengths )	#counts number of contigs / scaffolds in this assembly
	total_number_of_bases = sum( contig_lengths )	#counts all bases in the assembyl
	mean_contig_length = total_number_of_bases / number_of_contigs	#average contig lengths
	minimal_contig_length = min( contig_lengths )
	maximal_contig_length = max( contig_lengths )

	# --- sort list of contig length decreasing --- #
	sorted_contig_lengths = sorted( contig_lengths )[::-1]	#invert to get it decreasing
	N25 = False
	N50 = False
	N75 = False
	N90 = False
	
	cum_length = total_number_of_bases
	
	for contig_length in sorted_contig_lengths:
		cum_length -= contig_length
		if cum_length <= 0.1 * total_number_of_bases:
			if not N90:
				N90 = contig_length
		elif cum_length <= 0.25 * total_number_of_bases:
			if not N75:
				N75 = contig_length
		elif cum_length <= 0.5 * total_number_of_bases:
			if not N50:
				N50 = contig_length
		elif cum_length <= 0.75 * total_number_of_bases:
			if not N25:
				N25 = contig_length
	
	# --- sort list of expression-filtered contig length decreasing --- #
	sorted_exp_contig_lengths = sorted( exp_contig_lengths )[::-1]	#invert to get it decreasing
	EN25 = False
	EN50 = False
	EN75 = False
	EN90 = False
	
	total_exp_bases = sum(  exp_contig_lengths )
	cum_exp_length = sum(  exp_contig_lengths )
	
	
	for contig_length in sorted_exp_contig_lengths:
		cum_exp_length -= contig_length
		if cum_exp_length <= 0.1 * total_exp_bases:
			if not EN90:
				EN90 = contig_length
		elif cum_exp_length <= 0.25 * total_exp_bases:
			if not EN75:
				EN75 = contig_length
		elif cum_exp_length <= 0.5 * total_exp_bases:
			if not EN50:
				EN50 = contig_length
		elif cum_exp_length <= 0.75 * total_exp_bases:
			if not EN25:
				EN25 = contig_length
	
	
	stats = { 	'number_of_contigs': number_of_contigs,
			'mean_contig_length': mean_contig_length,
			'minimal_contig_length': minimal_contig_length,
			'maximal_contig_length': maximal_contig_length,
			'total_number_of_bases': total_number_of_bases,
			'number_of_bases_without_N': number_of_bases_without_N,
			'gc_content': float( number_of_gc ) /number_of_bases_without_N,
			'N25': N25,
			'N50': N50,
			'N75': N75,
			'N90': N90,
			'EN25': EN25,
			'EN50': EN50,
			'EN75': EN75,
			'EN90': EN90
		 }
	
	sys.stdout.write(  "calculation of formal assembly stats done.\n" )
	sys.stdout.flush()
	return stats


def write_NExp_evaluation_to_file( stats_outputfile, formal_stats, assembly_name, percent_cutoff ):
	"""! @brief writes all calculated evaluation results to file
		
		@param prefix (string) path to the ouput loction of all files
		
		@param outputfile (string)) only name of file for result output
		
		@param formal_stats (dictionary) contains some statistics about the assembly
		
		@param assembly_name (string) the name of the currently processed assembly
	"""
	
	sys.stdout.write( "writing results to file ... please wait!\n" )
	sys.stdout.flush()
	
	with open( stats_outputfile, 'w' ) as out:
		out.write( 'assembly name: ' + assembly_name + '\n\n' )
		
		out.write( 'number of contigs:\t' + str( formal_stats['number_of_contigs'] ) + '\n' )
		out.write( 'average contig length:\t' + str( formal_stats['mean_contig_length'] ) + '\n' )
		out.write( 'minimal contig length:\t' + str( formal_stats['minimal_contig_length'] ) + '\n' )
		out.write( 'maximal contig length:\t' + str( formal_stats['maximal_contig_length'] ) + '\n\n' )
		
		out.write( 'total number of bases:\t' + str( formal_stats['total_number_of_bases'] ) + '\n' )
		out.write( 'total number of bases without Ns:\t' + str( formal_stats['number_of_bases_without_N'] ) + '\n' )
		out.write( 'GC content:\t' + str( formal_stats['gc_content'] ) + '\n\n' )
		
		out.write( 'N25:\t' + str( formal_stats['N25'] ) + '\n' )
		out.write( 'N50:\t' + str( formal_stats['N50'] ) + '\n' )
		out.write( 'N75:\t' + str( formal_stats['N75'] ) + '\n' )
		out.write( 'N90:\t' + str( formal_stats['N90'] ) + '\n\n' )
		
		out.write( "E" + str( percent_cutoff ) + 'N25:\t' + str( formal_stats['EN25'] ) + '\n' )
		out.write( "E" + str( percent_cutoff ) + 'N50:\t' + str( formal_stats['EN50'] ) + '\n' )
		out.write( "E" + str( percent_cutoff ) + 'N75:\t' + str( formal_stats['EN75'] ) + '\n' )
		out.write( "E" + str( percent_cutoff ) + 'N90:\t' + str( formal_stats['EN90'] ) + '\n\n' )
		
	sys.stdout.write( "all results written to file.\n" )
	sys.stdout.flush()


def clean_assembly_file( input_file, output_file, cutoff ):
	"""! @brief removes small contigs and cleans contig name to 'contig_<INTEGER>' """
	
	sys.stdout.write( "cleaning contig names and removing small contigs ... please wait!\n" )
	sys.stdout.flush()
	with open( output_file, "w" ) as out:
		with open( input_file, "r" ) as f:
			line = f.readline()
			try:
				try:
					try:
						try:
							try:
								try:
									header = re.findall( "contig_\d+", line )[0]
								except:
									header = re.findall( "contig\d+", line )[0]
							except:
								header = re.findall( "scaffold\d+", line )[0]
						except:
							header = re.findall( "C\d+", line )[0]
					except:
						header = re.findall( "NODE_\d+", line )[0]
				except:
					header = re.findall( "seq\d+", line )[0]
			except:
				header = line.strip()[1:]
			seq = []
			while line:
				if line[0] == '>':
					seq = "".join( seq )
					if len( seq ) >= cutoff:
						out.write( '>' + header + '\n' + seq + '\n' )
					seq = []
					try:
						try:
							try:
								try:
									try:
										try:
											header = re.findall( "contig_\d+", line )[0]
										except:
											header = re.findall( "contig\d+", line )[0]
									except:
										header = re.findall( "scaffold\d+", line )[0]
								except:
									header = re.findall( "C\d+", line )[0]
							except:
								header = re.findall( "NODE_\d+", line )[0]
						except:
							header = re.findall( "seq\d+", line )[0]
					except:
						header = line.strip()[1:]
				else:
					seq.append( line.strip() )
				line = f.readline()
			seq = "".join( seq )
			if len( seq ) >= cutoff:
				out.write( '>' + header + '\n' + seq + '\n' )


def load_expression( exp_file, percent_cutoff ):
	"""! @brief load expression """
	
	# --- load expression --- #
	exp = {}
	with open( exp_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			value = sum( map( float, parts[1:] ) ) / len( parts[1:] )
			exp.update( { parts[0]: value } )
			line = f.readline()
	
	# --- get top X % of contigs based on expression --- #
	genes_for_sorting = []
	for key in exp.keys():
		genes_for_sorting.append( { 'contig': key, 'val': exp[ key ] } )
	total_exp = sum( exp.values() )
	cutoff = ( total_exp / 100.0 ) * percent_cutoff
	contigs_sorted_by_exp = sorted( genes_for_sorting, key=itemgetter('val') )[::-1]
	valid_contigs = {}
	counter = 0
	
	for each in contigs_sorted_by_exp:
		if counter < cutoff:
			counter += each['val']
			valid_contigs.update( { each['contig']: None } )
		else:
			return valid_contigs
	return valid_contigs


def main( arguments ):
	"""! @brief runs all parts of this script """
	
	if '--input' in arguments:
		raw_assembly_file = arguments[ arguments.index( '--input' ) + 1 ]
	else:
		raw_assembly_file = arguments[ arguments.index( '--in' ) + 1 ]
	
	if '--min_contig_len' in arguments:
		cutoff = int( arguments[ arguments.index( '--min_contig_len' ) + 1 ] )
	elif '--min' in arguments:
		cutoff = int( arguments[ arguments.index( '--min' ) + 1 ] )
	else:
		cutoff = 1000
	
	if '--out' in arguments:
		output_dir = arguments[ arguments.index( '--out' ) + 1 ]
		if output_dir[-1] != '/':
			output_dir += "/"
		if not os.path.exists( output_dir ):
			os.makedirs( output_dir )
		clean_assembly_filename = output_dir + raw_assembly_file.split('/')[-1] + '_trimmed.fasta'
		stats_outputfile = output_dir + raw_assembly_file.split('/')[-1] + '_stats.txt'
	else:
		clean_assembly_filename = raw_assembly_file + '_trimmed.fasta'
		stats_outputfile = clean_assembly_filename.replace( "_trimmed.fasta", "_stats.txt" )
	
	if '--exp' in arguments:
		exp_file = arguments[ arguments.index( '--exp' ) + 1 ]
		percent_cutoff = 90
		expX = load_expression( exp_file, percent_cutoff )
	else:
		percent_cutoff = 90
		expX = {}
	
	# --- cleaning assembly --- #
	clean_assembly_file( raw_assembly_file, clean_assembly_filename, cutoff )
	
	# --- calculating assembly stats --- #
	formal_assembly_stats = calculate_formal_contig_stats_expX( clean_assembly_filename, expX )
	assembly_name = '.'.join( clean_assembly_filename.split('/')[-1].split('.')[:-1] )	
	
	# ---- write all results of the evaluation to file --- #
	write_NExp_evaluation_to_file( stats_outputfile, formal_assembly_stats, assembly_name, percent_cutoff )


if '--input' in sys.argv:
	main( sys.argv )
elif '--in' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
