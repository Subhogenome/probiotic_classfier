🧬 ePredict: Genomic Classifier for Probiotic and Pathogenic Detection

ePredict is a Dockerized machine learning pipeline that classifies genomic sequences (FASTA format) into Probiotic or Pathogenic categories with up to 92% accuracy. Designed for microbiome and bioinformatics researchers, this tool processes raw genome sequences and visualizes prediction probabilities in real-time.

🚀 Key Features

✅ Accepts both draft and complete genomes

📦 Built with LightGBM, BioPython, Streamlit, and Plotly

🔬 Extracts k-mer (1–8) frequency vectors from DNA

📊 Visualizes predicted probabilities using interactive bar plots

🐳 Dockerized for consistent deployment: subho08/epredict:latest

🧪 (Optional) AMR resistance analysis with AMRfinder

🧑‍🔬 How It Works

Upload a .fasta or .fna file

Choose whether your genome is draft or complete

Preprocessing removes invalid characters and formats your sequence

K-mer encoding generates frequency features

Trained LightGBM model predicts Probiotic or Pathogenic

Interactive Plotly chart shows result, with CSV export option

🧪 Tech Stack

Component

Technology

Model

LightGBM (92% Accuracy)

UI

Streamlit + Plotly

Genomics

BioPython

Packaging

Docker

Visuals

Plotly Express

Sequence Parsing

FASTA file input

📦 Run via Docker

# Pull from Docker Hub
docker pull subho08/epredict:latest

# Run the container
docker run -p 8501:8501 subho08/epredict:latest

# Open in browser
http://localhost:8501

📁 Project Structure

epredict/
├── finalized_model.sav       # Trained LightGBM model
├── epredict_app.py           # Main Streamlit application
├── requirements.txt
├── Dockerfile
└── output/                   # Generated results (fasta, CSV, plots)

📈 Sample Output

Genome Name

Probiotic Score

Pathogenic Score

Lactobacillus_plantarum

0.94

0.06

Clostridium_difficile

0.03

0.97

🧪 Use Cases

💊 Probiotic product development

🧬 Gut microbiome classification

🔬 Pathogen screening for drug discovery

🧠 ML research on biological sequences

📄 Citation / Credit

This tool was created by Subho Chatterjee as part of a bioinformatics + machine learning research pipeline.

GitHub: subho08

DockerHub: subho08/epredict

If you use this project in your work, feel free to cite or credit in your documentation 🙌

📬 Contributing / Feedback

Open to feature requests, performance improvements, or suggestions. Submit a GitHub issue or pull request if you’d like to contribute!

🧠 Future Work

✅ Batch processing of multiple genomes

📊 Integrate SHAP explainability

🔄 Add AMR gene output overlay

☁️ Deploy via Streamlit Cloud or Hugging Face Spaces

Let the genome speak. Predict with precision. 🧬

