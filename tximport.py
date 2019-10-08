## 
import gffutils
import sys
import os
import tempfile
import optparse
import subprocess

CHUNK_SIZE = 2**20

# CALL:
#tximport.py -gff_file $gff_source.gff_key.fields.path --out_file $gene_level_values --sample sample1 --sample sample2 --sample sampleN


# create transcript - gene table in tmp_dir
def get_tx2gen_table(gff_file,tmp_file):
    out_file = open(tmp_file, 'w')
    db = gffutils.create_db(gff_file, 'test.db', force=True)
    db = gffutils.FeatureDB('test.db', keep_order=True)
    for gene in db.features_of_type('gene'):
        for child in db.children(gene, order_by='start'):
           if child.featuretype=='mRNA' or child.featuretype=='rRNA':
              out_file.write(child.id + '\t' + gene.id + '\n')


def main():
    #parse samples list, gff file and output
    parser = optparse.OptionParser()
    parser.add_option('-s', '--sample', action='append', dest='samples_list',   help='Add repeated values to a list' )
    parser.add_option( '-o', '--out_file', dest='out_file', action='store', type="string", default=None )
    parser.add_option( '-i', '--input_format', dest='input_format', action='store', type="string", default=None )  
              # if none then there should also be parameters for:
    #            <param name="tx_id_col" type="text" label="Name of the txID columns"/>
    #            <param name="abundance_col" type="text" label="Name of the abundance column"/>
    #            <param name="counts_col" type="text" label="Name of the counts column"/>
    #            <param name="length_col" type="text" label="Name of the length column"/>

    parser.add_option( '-c', '--col_name', dest='col_name', action='store', type="string", default=None )
    parser.add_option( '-f', '--gff_file', dest='gff_file', action='store', type="string", default=None )
    parser.add_option( '-p', '--base_dir', dest='base_dir', action='store', type="string", default=None )
    parser.add_option( '-t', '--tx2gene', dest='tx2gene_table', action='store', type="string", default=None )
    parser.add_option( '-m', '--out_mode', dest='out_mode', action='store', type="string", default=None )  # out_mode is either 'individual' 'merge'
    (options, args) = parser.parse_args()
 
    #gff_file=options.gff_file
    #samples=options.samples_list
    
    #create tmp file to save tx-gene table
    if options.tx2gene_table:
	gene_tx_table=options.tx2gene_table
    else:
    	tmp_dir = tempfile.mkdtemp( prefix='tmp-gene-tx-table' )
	# convert to SAM
    	gene_tx_table = os.path.join( tmp_dir, 'gene_table.csv' )
	get_tx2gen_table(options.gff_file, gene_tx_table)
    
    scriptR=  os.path.join( options.base_dir, 'tximport.r') 
    #scriptsR = options.base_dir + '/' + 'tximport.r'
    # prepare call to R script
    args= ['Rscript']
    args.append(scriptR)
    args.append(gene_tx_table)
    args.append(options.out_file)
    #args.append(options.ou)
    #if options.out_mode:
    #for o in options.out_file:
    #   args.append(o)
    #args.append(options.out_file)
    for s in options.samples_list:
       args.append(s)
    #print args
    tmp_stderr = tempfile.NamedTemporaryFile( prefix = "tmp-stderr" )
    return_code = subprocess.call( args=args, shell=False, stderr=tmp_stderr.fileno() )
    #return_code = subprocess.call( args=args, shell=False, stderr=None)
    if return_code:
        tmp_stderr.flush()
        tmp_stderr.seek(0)
        print >> sys.stderr, "Error in process call"
        while True:
            chunk = tmp_stderr.read( CHUNK_SIZE )
            if not chunk:
                break
            sys.stderr.write( chunk )
        sys.exit( return_code )
    tmp_stderr.close()

   
if __name__ == "__main__":
    main()


