from flask import Flask, render_template, request, jsonify
from werkzeug.utils import secure_filename
import os
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max file size

# Create uploads directory if it doesn't exist
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

# Example: First 1000 bases of SARS-CoV-2 reference genome (for demo)
REFERENCE_GENOME = (
    "ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT"
    "GTTCTCTAAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACT"
    "CACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCT"
    "TCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTCG"
    "GATCTCTTGTAGATCTGTTCTCTAAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCT"
    "GCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACG"
    "AGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCA"
    "GCACATCTAGGTTTCGATCTCTTGTAGATCTGTTCTCTAAAACGAACTTTAAAATCTGTGT"
    "GGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGT"
    "CGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGT"
    "TGCAGCCGATCATCAGCACATCTAGGTTTCGATCTCTTGTAGATCTGTTCTCTAAAACGAA"
)

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in {'fasta', 'fa'}

def find_mutations(reference_seq, sample_seq):
    mutations = []
    for i, (ref_base, sample_base) in enumerate(zip(reference_seq, sample_seq)):
        if ref_base != sample_base and ref_base != '-' and sample_base != '-':
            mutations.append({
                'position': i + 1,
                'reference': ref_base,
                'sample': sample_base
            })
    return mutations

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_files():
    try:
        if 'files[]' not in request.files:
            logger.error("No files part in request")
            return jsonify({'error': 'No files provided'}), 400
        
        files = request.files.getlist('files[]')
        logger.debug(f"Number of files received: {len(files)}")
        
        if len(files) < 2 or len(files) > 3:
            logger.error(f"Invalid number of files: {len(files)}")
            return jsonify({'error': 'Please upload 2-3 FASTA files'}), 400

        results = []
        for file in files:
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                logger.debug(f"Processing file: {filename}")
                filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                file.save(filepath)
                
                try:
                    record = next(SeqIO.parse(filepath, "fasta"))
                    sample_name = record.id
                    sample_seq = str(record.seq)
                    logger.debug(f"Sample {sample_name} sequence length: {len(sample_seq)}")
                    
                    # Truncate or pad sample to match reference length for demo
                    sample_seq = sample_seq[:len(REFERENCE_GENOME)].ljust(len(REFERENCE_GENOME), '-')
                    aligner = PairwiseAligner()
                    alignments = aligner.align(REFERENCE_GENOME, sample_seq)
                    best_alignment = alignments[0]
                    mutations = find_mutations(best_alignment[0], best_alignment[1])
                    logger.debug(f"Found {len(mutations)} mutations in {sample_name}")
                    
                    results.append({
                        'sample_name': sample_name,
                        'mutations': mutations
                    })
                    os.remove(filepath)
                except Exception as e:
                    logger.error(f"Error processing {filename}: {str(e)}")
                    return jsonify({'error': f'Error processing {filename}: {str(e)}'}), 400
            else:
                logger.error(f"Invalid file: {file.filename if file else 'None'}")
                return jsonify({'error': f'Invalid file: {file.filename if file else "None"}'}), 400
        
        return jsonify({'results': results})
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
        return jsonify({'error': f'Unexpected error: {str(e)}'}), 500

if __name__ == '__main__':
    app.run(debug=True) 