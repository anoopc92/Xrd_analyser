🔬 XRD Analysis Pipeline

A comprehensive, end-to-end X-ray Diffraction (XRD) Analysis Pipeline built with Python and deployed as an interactive web app using Streamlit. This tool processes XRD data from various file formats, performs robust peak detection, Voigt-profile fitting, Scherrer and Williamson-Hall analysis, and phase matching with the Materials Project database. Visualize results with interactive Plotly charts and download detailed JSON reports.
🚀 Try the Live App on Streamlit Cloud (Insert your Streamlit app URL here after deployment)

✨ Features

Heterogeneous File Support: Load XRD data from CSV, TXT, XY, UXD, and RAW formats.
Robust Preprocessing: Smooth signals, remove backgrounds, and handle noisy data with asymmetric least squares.
Peak Detection & Fitting: Identify peaks and fit them with Voigt profiles using lmfit.
Crystallite Analysis: Compute crystallite size (Scherrer) and size/strain (Williamson-Hall).
Phase Matching: Match experimental peaks against Materials Project references (optional, requires API key).
Interactive Visualizations: Generate multi-panel Plotly plots for patterns, peak fits, Williamson-Hall analysis, and phase matches.
Provenance Tracking: Log all analysis steps for reproducibility.
Streamlit Web App: Upload files, tweak hyperparameters, and view results in a user-friendly interface.
Batch Processing: Analyze multiple files at once (CLI mode).
Hyperparameter Optimization: Fine-tune preprocessing and peak detection parameters.


📸 Screenshots
Coming soon: Add a screenshot of the Streamlit app showing the interactive plot or phase match table!

🛠️ Installation
Prerequisites

Python 3.8+
Git
(Optional) Materials Project API key for phase matching (get one here)

Steps

Clone the Repository:
git clone https://github.com/YOUR_USERNAME/xrd-analysis-app.git
cd xrd-analysis-app


Create a Virtual Environment (recommended):
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate


Install Dependencies:
pip install -r requirements.txt

The requirements.txt includes:
streamlit
numpy
scipy
lmfit
pymatgen
mp-api
plotly


(Optional) Set Up Materials Project API Key:

For phase matching, add your API key to .streamlit/secrets.toml:MP_API_KEY = "your_api_key_here"


Alternatively, input the API key via the Streamlit app.




🚀 Usage
Running the Streamlit App Locally

Start the app:streamlit run app.py


Open your browser to http://localhost:8501.
Upload an XRD file (CSV, TXT, XY, UXD, or RAW).
(Optional) Enter a Materials Project API key, chemical formulas (e.g., Si,SiO2), or material IDs (e.g., mp-149).
Adjust hyperparameters (e.g., smoothing window, peak prominence) in the sidebar.
Click "Run Analysis" to view results, including:
Peak detection summary
Crystallite size and strain
Phase matching table (if API key provided)
Interactive Plotly plot
Downloadable JSON report



Running in CLI Mode
For batch processing or programmatic use:
from xrd_pipeline import XRDAnalysisPipeline, batch_analyze

pipeline = XRDAnalysisPipeline(mp_api_key="YOUR_API_KEY")
state = pipeline.run(
    filepath="sample_xrd.csv",
    material_ids=["mp-149"],
    output_dir="xrd_outputs"
)

Example Dataset
A synthetic dataset (sample_xrd.csv) is included for testing:

Format: Two columns (2theta,intensity)
2θ range: 10° to 90°
Peaks at 20°, 30°, 45°, 60°, and 75° with noise and background


🌐 Deploying to Streamlit Cloud

Push to GitHub:

Ensure app.py, requirements.txt, and (optionally) sample_xrd.csv are in your repo.
Push to a public GitHub repository:git add .
git commit -m "Initial commit"
git push origin main




Deploy on Streamlit Cloud:

Go to share.streamlit.io and sign in with GitHub.
Select your repo (xrd-analysis-app), branch (main), and main file (app.py).
Deploy! You'll get a URL like https://xrd-analysis-app.streamlit.app.
Update this README with the live URL.


Secure API Key:

Add MP_API_KEY to Streamlit Cloud's secrets management instead of hardcoding it.




📊 Example Output

Peaks Detected: Number of peaks found with positions, intensities, and FWHM.
Crystallite Analysis:
Scherrer size (e.g., 50.23 nm)
Williamson-Hall size (e.g., 48.12 nm) and strain (e.g., 0.0023)


Phase Matches (if API key provided):
Example: Si (mp-149) with a match score of 0.92 (9/10 peaks matched)


Plots: Multi-panel Plotly figure showing the XRD pattern, peak fits, Williamson-Hall plot, and phase match table.


🧪 Testing

Use the provided sample_xrd.csv to test the app or CLI.
Upload the file in the Streamlit app or run:pipeline.run("sample_xrd.csv", material_ids=["mp-149"], output_dir="xrd_outputs")


Check the xrd_outputs folder for JSON reports and HTML plots.


🤝 Contributing
Contributions are welcome! To contribute:

Fork the repository.
Create a feature branch: git checkout -b feature-name.
Commit changes: git commit -m "Add feature".
Push to the branch: git push origin feature-name.
Open a pull request.

Please include tests and update documentation as needed.

📜 License
This project is licensed under the MIT License. See the LICENSE file for details.

📬 Contact
For questions or suggestions, open an issue or contact [YOUR_EMAIL_OR_GITHUB_PROFILE].
🌟 Star this repo if you find it useful!
