'''
The following is the Dockerfile for local deployment
'''
# FROM continuumio/miniconda3

# # Create conda environment and install RDKit + FastAPI stack
# RUN conda create -y -n molwatch python=3.10 \
#   && conda install -y -n molwatch -c conda-forge rdkit fastapi uvicorn prometheus-fastapi-instrumentator starlette \
#   && conda clean --all --yes

# # Use conda environment shell for subsequent commands
# SHELL ["conda", "run", "-n", "molwatch", "/bin/bash", "-c"]

# # Set working directory
# WORKDIR /app

# # Copy only necessary files (avoid Docker cache invalidation)
# COPY requirements.txt /app/requirements.txt
# WORKDIR /app
# RUN pip install -r requirements.txt

# # Expose port
# EXPOSE 8000

# # Run the FastAPI app using uvicorn from the conda environment (port: 8000 if local, 8080 for cloud-based)
# ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "molwatch", "uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8080"]

'''
The following is the Dockerfile for cloud-based deployment
'''

FROM continuumio/miniconda3

# Create environment and install packages
RUN conda create -y -n molwatch python=3.10 \
  && conda install -y -n molwatch -c conda-forge rdkit fastapi uvicorn prometheus-fastapi-instrumentator starlette \
  && conda clean --all --yes

# Set environment so that future RUN commands use conda env
ENV PATH /opt/conda/envs/molwatch/bin:$PATH

# Set working directory
WORKDIR /app

# Copy app code
COPY . /app

# Install any pip-based dependencies
RUN pip install -r requirements.txt

# Expose the correct port (Cloud Run uses 8080)
EXPOSE 8080

# Start script (see below)
COPY start.sh /start.sh
RUN chmod +x /start.sh

# Use shell entrypoint for proper env activation
ENTRYPOINT ["/bin/bash", "/start.sh"]