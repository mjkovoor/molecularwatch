# MolecularWatch

MolecularWatch is a cloud-native application for predicting drug-like molecular properties from SMILES strings. It features a FastAPI backend for scientific computation, a React frontend for user interaction, and CI/CD pipelines for deployment on Google Cloud Run.

---

## Features

- Solubility prediction using the ESOL method
- Lipinski rule-based drug-likeness scoring
- Molecular structure rendering via RDKit
- FastAPI backend with REST endpoints
- React-based frontend interface
- Metrics endpoint compatible with Prometheus
- Docker-based deployment with GitHub Actions CI/CD

---

## Tech Stack

- **Backend:** Python, FastAPI, RDKit
- **Frontend:** React, TypeScript, Vite
- **Deployment:** Docker, GitHub Actions, Google Cloud Run
- **Monitoring:** Prometheus (via FastAPI Instrumentator)

---

## Local Development

### Backend

```bash
# Set up the conda environment
conda create -n molwatch python=3.10
conda activate molwatch
conda install -c conda-forge rdkit fastapi uvicorn prometheus-fastapi-instrumentator starlette
pip install -r requirements.txt

# Run the backend
uvicorn app.main:app --reload