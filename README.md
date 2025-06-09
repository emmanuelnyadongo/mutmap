# Viral Genome Mutation Visualizer

A web application for visualizing mutations in viral genome sequences, with a focus on SARS-CoV-2 and other RNA viruses.

## Bioinformatics Significance

This tool addresses several key aspects of viral genome analysis:

1. **Mutation Tracking**: Essential for monitoring viral evolution and identifying emerging variants
2. **Sequence Alignment**: Compares multiple viral sequences against a reference genome
3. **SNP Detection**: Identifies single nucleotide polymorphisms (SNPs) that may affect viral behavior
4. **Visualization**: Provides intuitive representation of mutation patterns

### Applications in Infectious Disease Research

- **Variant Surveillance**: Track emerging viral variants and their mutations
- **Outbreak Analysis**: Compare sequences from different outbreaks to identify transmission patterns
- **Vaccine Development**: Monitor mutations that might affect vaccine efficacy
- **Treatment Planning**: Identify mutations that could impact antiviral drug effectiveness

## Features

- Upload 2-3 viral genome sequences in FASTA format
- Automatic alignment against reference genome
- SNP detection and visualization
- Interactive genome visualization
- Detailed mutation table with position and base change information

## Technical Requirements

- Python 3.8+
- Biopython
- Modern web browser

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/viral-mutation-visualizer.git
cd viral-mutation-visualizer
```

2. Install Python dependencies:
```bash
pip install -r requirements.txt
```

3. Run the application:
```bash
python app.py
```

4. Open your browser and navigate to `http://localhost:5000`

## Usage

1. Prepare your FASTA files containing viral genome sequences
2. Upload 2-3 sequences through the web interface
3. View the mutation visualization and detailed results table
4. Download results if needed

## Project Structure

```
viral-mutation-visualizer/
├── app.py              # Main Flask application
├── requirements.txt    # Python dependencies
├── static/            # Static files (CSS, JS)
│   ├── css/
│   └── js/
├── templates/         # HTML templates
└── utils/            # Utility functions for sequence analysis
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details. 