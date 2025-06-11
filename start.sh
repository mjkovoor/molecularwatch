#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate molwatch
exec uvicorn app.main:app --host 0.0.0.0 --port ${PORT:-8080}