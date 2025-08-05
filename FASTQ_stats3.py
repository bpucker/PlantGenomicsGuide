### Boas Pucker ###
### b.pucker@tu-bs.de ###

__version__ = "v0.41"


__usage__ = """
		Calculation of FASTQ statististics (""" + __version__ + """):
		
		python3 FASTQ_stats3.py
		--in <FULL_PATH_TO_FASTQ_FILE> |	--in_dir <FULL_PATH_TO_DIRECTORY>
		
		optional:
		--rfig <READ_LEN_HIST_FIGURE_FILE>
		--cutoff <READ_LEN_CUTOFF_FOR_PLOT_IN_KB>[100]
		--qfig <QUALITY_VS_READ_LEN_FIGURE_FILE>
		--lencut <UPPER_READ_LENGTH_CUTOFF_FOR_PLOT_IN_KB>[200]
		--qualcut <UPPER_QUAL_CUTOFF_FOR_PLOT>[40]
		--quality <ACTIVATES_QUALITY_ANALYSIS>
		
		bug reports and feature requests: b.pucker@tu-bs.de
		"""

import sys, gzip, glob
try:
	import matplotlib.pyplot as plt
	import matplotlib.ticker as ticker
	import pandas as pd
	import seaborn as sns
except ImportError:
	sys.stdout.write( "WARNING: matplotlib import failed - figure plotting not possible.\n" )
	sys.stdout.flush()


# --- end of imports --- #

def calculate_n50( total_length ):
	"""! @brief calculate N50 based on list of given read lengths """
	
	total = sum( total_length )
	half = total / 2.0
	sorted_lengths = sorted( total_length )[::-1]
	counter = 0
	i = 0
	for i in range( len( sorted_lengths ) ):
		counter += sorted_lengths[ i ]
		if counter >= half:
			return sorted_lengths[ i ]
	sys.exit( "ERROR: no read lengths detected." )


def calc_avg_qual( qual ):
	"""! @brief quality """
	
	quality_values = []
	for q in qual:
		quality_values.append( ord( q ) - 33.0 )
	return sum( quality_values ) / float( len( quality_values ) )


def analyze_FASTQ( filename, qual_status ):
	"""! @brief analysis of FASTQ file """
	
	gzip_state = False
	try:
		file_extension = filename.split('.')[-1]
		if file_extension in [ "gz", "gzip", "GZ", "GZIP" ]:
			gzip_state = True
	except:
		pass
	
	if not gzip_state:
		with open( filename, "r" ) as f:
			total_length = []
			average_quality = []
			total_GC = []
			line = f.readline()	#header
			while line:
				seq = f.readline().strip().upper()
				total_length.append( len( seq ) )
				total_GC.append( seq.count('C') )
				total_GC.append( seq.count('G') )
				f.readline()	#useless line
				qual = f.readline()	#quality line
				if qual_status:
					average_quality.append( calc_avg_qual( qual ) )
				line = f.readline()
			
			total_len = sum( total_length )
			total_gc = sum( total_GC )
			n50 = calculate_n50( total_length )
			
			sys.stdout.write(  filename + "\n" )
			sys.stdout.write(  "total number of nucleotides:\t" + str( total_len / 1000000000.0 ) + " Gbp\n" )
			sys.stdout.write(  "N50: " + str( n50 ) + "\n" )
			sys.stdout.write( "number of reads: " + str( len( total_length ) ) + "\n" )
			
			sys.stdout.write(  "average read length:\t" + str( total_len / float( len( total_length ) ) ) + "\n" )
			sys.stdout.write(  "GC content:\t" + str( total_gc / float( total_len ) ) + "\n" )
			
			if qual_status:
				if len( average_quality ) > 0:
					sum_avg_qual = sum( average_quality ) / len( average_quality )
				else:
					sum_avg_qual = 0
				sys.stdout.write( "average read quality score: " + str( sum_avg_qual ) + "\n" )
			
			sys.stdout.flush()
		
	else:
		with gzip.open( filename, "rb" ) as f:
			total_length = []
			average_quality = []
			total_GC = []
			line = f.readline().decode("utf-8")	#header
			while line:
				seq = f.readline().decode("utf-8").strip().upper()
				total_length.append( len( seq ) )
				total_GC.append( seq.count('C') )
				total_GC.append( seq.count('G') )
				f.readline()	#useless line
				qual = f.readline().decode("utf-8")	#quality line
				if qual_status:
					average_quality.append( calc_avg_qual( qual ) )
				line = f.readline().decode("utf-8")
			
			total_len = sum( total_length )
			total_gc = sum( total_GC )
			n50 = calculate_n50( total_length )
			
			sys.stdout.write( filename + "\n" )
			sys.stdout.write( "total number of nucleotides:\t" + str( total_len / 1000000000.0 ) + " Gbp\n" )
			sys.stdout.write( "N50: " + str( n50 ) + "\n" )
			sys.stdout.write( "number of reads: " + str( len( total_length ) ) + "\n" )
			
			sys.stdout.write( "average read length:\t" + str( total_len / float( len( total_length ) ) ) + "\n" )
			sys.stdout.write( "GC content:\t" + str( total_gc / float( total_len ) ) + "\n" )
			
			if qual_status:
				if len( average_quality ) > 0:
					sum_avg_qual = sum( average_quality ) / len( average_quality )
				else:
					sum_avg_qual = 0
				sys.stdout.write( "average read quality score: " + str( sum_avg_qual ) + "\n" )
			
			sys.stdout.flush()
	return total_length, average_quality


def generate_read_len_hist( total_length, fig_file, len_cutoff ):
	"""! @brief generate histogram of read length distribution """
	
	# --- preprocess data --- #
	bins = []
	for i in range( len_cutoff ):
		bins.append( [] )
	for length in total_length:
		if length > (len_cutoff*1000):
			bins[-1].append( length )
		else:
			bins[ int( length/1000.0 ) ].append( length )
	df_input = []
	labels = []
	for idx, each in enumerate( bins ):
		df_input.append( [ str( idx ) + "-" + str(idx+1) + "kb", sum( each ) / 1000000.0 ] )
		labels.append( str( idx ) + "-" + str(idx+1) + "kb" )
	df = pd.DataFrame( data=df_input, index=labels, columns=["label", "value"] )
	
	# --- generate figure --- #
	fig, ax = plt.subplots()
	sns.barplot( data=df, x="label", y="value", color='steelblue' )
	ax.set_xlabel( "read length [kb]" )
	ax.set_ylabel( "total sequence amount [Mbp]" )
	
	ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
	ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
	ax.tick_params(axis='x', rotation=90)
	
	plt.tight_layout()
	
	fig.savefig( fig_file, dpi=300 )


def generate_quality_vs_read_len_figure( total_length, average_quality, figfile, max_qual_cut, max_len_cut ):
	"""! @brief generate quality vs. read length figure """
	
	xvalues, yvalues = [], []
	for idx, val in enumerate( average_quality ):
		xvalues.append( min( [ total_length[ idx ] / 1000.0, max_len_cut ] ) )	#unit is kb
		yvalues.append( min( [ val, max_qual_cut ] ) )
	
	sns.set_theme(style="dark")
	fig, ax = plt.subplots()
	sns.scatterplot(x=xvalues, y=yvalues, s=5, color=".15")
	sns.histplot(x=xvalues, y=yvalues, bins=100, pthresh=.1, cmap="Greens")	#cmap="mako"
	sns.kdeplot(x=xvalues, y=yvalues, levels=5, color="r", linewidths=1)
	
	ax.set_xlabel( "read length [kb]" )
	ax.set_ylabel( "quality (Phred score)" )
	ax.set_title( "number of reads: " + str( len( xvalues ) ) )
	
	ax.set_xlim( 0, max_len_cut+1 )
	ax.set_ylim( 0, max_qual_cut+1 )
	
	plt.tight_layout()
	
	fig.savefig( figfile, dpi=300 )


def main( arguments ):
	"""! @brief runs everything """
	
	if '--qfig' in arguments:
		qual_status = True
	else:
		if '--quality' in arguments:
			qual_status = True
		else:
			qual_status = False
	
	if '--in_file' in arguments or '--in' in arguments:	#single file mode
		if '--in_file' in arguments:
			input_file = arguments[ arguments.index( '--in_file' )+1 ]
		elif '--in' in arguments:
			input_file = arguments[ arguments.index( '--in' )+1 ]
		total_length, average_quality = analyze_FASTQ( input_file, qual_status )
		if '--rfig' in arguments:
			read_len_fig_file = arguments[ arguments.index( '--rfig' )+1 ]
			if '--cutoff' in arguments:
				len_cutoff = int( arguments[ arguments.index( '--cutoff' )+1 ] )
			else:
				len_cutoff = 100
			generate_read_len_hist( total_length, read_len_fig_file, len_cutoff )
		if '--qfig' in arguments:
			quality_vs_read_len_fig_file = arguments[ arguments.index( '--qfig' )+1 ]
			
			if '--qualcut' in arguments:
				max_qual_cut = int( arguments[ arguments.index( '--qualcut' )+1 ] )
			else:
				max_qual_cut = 40
			if '--lencut' in arguments:
				max_len_cut = int( arguments[ arguments.index( '--lencut' )+1 ] )
			else:
				max_len_cut = 200
			
			generate_quality_vs_read_len_figure( total_length, average_quality, quality_vs_read_len_fig_file, max_qual_cut, max_len_cut )
	
	else:	#folder analysis mode
		directory = arguments[ arguments.index( '--in_dir' )+1 ]
		if directory[-1] != '/':
			directory += "/"
		input_files = []
		extensions = [ ".fq", ".fastq", ".fq.gzip", ".fastq.gzip", ".fq.gz", ".fastq.gz", ".FQ", ".FASTQ", ".FQ.GZIP", ".FASTQ.GZIP", ".FQ.GZ", ".FASTQ.GZ" ]
		for extension in extensions:
			input_files += glob.glob( directory + '*' + extension )
		for filename in input_files:
			try:
				total_length, average_quality = analyze_FASTQ( filename, qual_status )
			except:
				sys.stdout.write( "ERROR while processing " + filename + "\n" )
				sys.stdout.flush()


if '--in_file' in sys.argv or '--in_dir' in sys.argv or '--in' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
