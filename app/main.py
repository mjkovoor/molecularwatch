from fastapi import FastAPI, Query
from prometheus_fastapi_instrumentator import Instrumentator
from starlette.responses import Response
from app.predict import predict_properties
import time

# Rest for visualization
from fastapi.responses import HTMLResponse
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import base64

# For Cross-Origin Resource Sharing issues
from fastapi.middleware.cors import CORSMiddleware


app = FastAPI()

Instrumentator().instrument(app).expose(app)

origins = [
    "http://localhost:5173",
    # add any other origins you want to allow
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
async def read_root():
    return {"message": "Welcome to MolecularWatch"}

@app.get("/predict")
async def predict(smiles: str = Query(..., description="Molecule in SMILES format")):
    return predict_properties(smiles)

@app.get("/health")
def health_check():
    return {"status": "ok"}

'''
This page offers a visualization of our molecule
'''
@app.get("/molecule_page", response_class=HTMLResponse)
def molecule_html(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES string.")
    
    # Generate PNG image in memory
    img = Draw.MolToImage(mol, size=(300, 300))
    buf = BytesIO()
    img.save(buf, format='PNG')
    img_bytes = buf.getvalue()
    img_base64 = base64.b64encode(img_bytes).decode("utf-8")

    # Return HTML page with image
    html = f"""
    <html>
        <head><title>Molecule Viewer</title></head>
        <body>
            <h2>Structure for SMILES: {smiles}</h2>
            <img src="data:image/png;base64,{img_base64}" alt="Molecule Image"/>
        </body>
    </html>
    """
    return HTMLResponse(content=html)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("app.main:app", host="0.0.0.0", port=8080) # port=8000 for local testing