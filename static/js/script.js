document.addEventListener('DOMContentLoaded', () => {
    const uploadForm = document.getElementById('uploadForm');
    const resultsSection = document.querySelector('.results-section');
    const genomeVisualization = document.getElementById('genomeVisualization');
    const mutationsTable = document.getElementById('mutationsTable').querySelector('tbody');

    uploadForm.addEventListener('submit', async (e) => {
        e.preventDefault();
        
        const files = document.getElementById('files').files;
        if (files.length < 2 || files.length > 3) {
            alert('Please select 2-3 FASTA files');
            return;
        }

        const formData = new FormData();
        for (let file of files) {
            formData.append('files[]', file);
        }

        try {
            const response = await fetch('/upload', {
                method: 'POST',
                body: formData
            });

            if (!response.ok) {
                throw new Error('Upload failed');
            }

            const data = await response.json();
            displayResults(data.results);
        } catch (error) {
            alert('Error: ' + error.message);
        }
    });

    function displayResults(results) {
        // Clear previous results
        genomeVisualization.innerHTML = '';
        mutationsTable.innerHTML = '';
        
        // Show results section
        resultsSection.style.display = 'block';

        // Calculate genome length (assuming SARS-CoV-2 length of ~30,000 bases)
        const genomeLength = 30000;
        const visualizationWidth = genomeVisualization.offsetWidth;

        // Create visualization for each sample
        results.forEach((result, sampleIndex) => {
            const sampleColor = getSampleColor(sampleIndex);
            
            result.mutations.forEach(mutation => {
                // Create mutation marker
                const marker = document.createElement('div');
                marker.className = 'mutation-marker';
                marker.style.left = `${(mutation.position / genomeLength) * 100}%`;
                marker.style.backgroundColor = sampleColor;
                marker.setAttribute('data-mutation', 
                    `${result.sample_name}: ${mutation.reference}→${mutation.sample} at position ${mutation.position}`);
                
                genomeVisualization.appendChild(marker);

                // Add to table
                const row = document.createElement('tr');
                row.innerHTML = `
                    <td>${result.sample_name}</td>
                    <td>${mutation.position}</td>
                    <td>${mutation.reference}</td>
                    <td>${mutation.sample}</td>
                `;
                mutationsTable.appendChild(row);
            });
        });
    }

    function getSampleColor(index) {
        const colors = ['#e74c3c', '#3498db', '#2ecc71'];
        return colors[index % colors.length];
    }
}); 