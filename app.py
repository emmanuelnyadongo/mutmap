from flask import Flask, render_template, request, jsonify
from werkzeug.utils import secure_filename
import os
from Bio import SeqIO
from Bio.Align import PairwiseAligner

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
    if 'files[]' not in request.files:
        return jsonify({'error': 'No files provided'}), 400
    
    files = request.files.getlist('files[]')
    if len(files) < 2 or len(files) > 3:
        return jsonify({'error': 'Please upload 2-3 FASTA files'}), 400

    results = []
    for file in files:
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)
            try:
                record = next(SeqIO.parse(filepath, "fasta"))
                sample_name = record.id
                sample_seq = str(record.seq)
                # Truncate or pad sample to match reference length for demo
                sample_seq = sample_seq[:len(REFERENCE_GENOME)].ljust(len(REFERENCE_GENOME), '-')
                aligner = PairwiseAligner()
                alignments = aligner.align(REFERENCE_GENOME, sample_seq)
                best_alignment = alignments[0]
                mutations = find_mutations(best_alignment[0], best_alignment[1])
                results.append({
                    'sample_name': sample_name,
                    'mutations': mutations
                })
                os.remove(filepath)
            except Exception as e:
                return jsonify({'error': f'Error processing {filename}: {str(e)}'}), 400
    
    return jsonify({'results': results})

if __name__ == '__main__':
    app.run(debug=True) 