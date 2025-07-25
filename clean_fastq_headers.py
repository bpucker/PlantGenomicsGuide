### Boas Pucker ###
### pucker@uni-bonn.de ###
### v0.15 ###

__usage__ = """
					python2 clean_fastq_headers.py
					--in <FASTQ_FILE (INPUT)>
					--out <FASTQ_FILE (OUTPUT)>
					
					WARNING: this script expects gzip compressed files
					"""

import gzip, os, sys

# --- end of imports --- #

def clean_fastq_header( fastq_file_in, fastq_file_out ):
	"""! @brief extraction of header and sequence from fastq file for fasta construction """
	
	if ".gz" in fastq_file_in:
		with open( fastq_file_out, "w" ) as out:
			with gzip.open( fastq_file_in, "rb" ) as f:
				line = f.readline()
				while line:
					out.write( line.split( '\t')[0] + "\n" )
					out.write( f.readline() )
					out.write( f.readline() )
					out.write( f.readline() )
					line = f.readline()
	else:
		with open( fastq_file_out, "w" ) as out:
			with open( fastq_file_in, "r" ) as f:
				line = f.readline()
				while line:
					out.write( line.split( '\t')[0] + "\n" )
					out.write( f.readline() )
					out.write( f.readline() )
					out.write( f.readline() )
					line = f.readline()


def main( arguments ):
	"""! @brief run everything """
	
	fastq_file_in = arguments[ arguments.index('--in')+1 ]
	fastq_file_out = arguments[ arguments.index('--out')+1 ]
	
	clean_fastq_header( fastq_file_in, fastq_file_out )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

