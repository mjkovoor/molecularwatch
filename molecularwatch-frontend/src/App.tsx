import { useState } from 'react';

type PredictionResult = {
  ["molecular weight"]: number;
  logP: number;
  H_bond_donors: number;
  H_bond_acceptors: number;
  TPSA: number;
  rotatable_bonds: number;
  aromatic_rings: number;
  esol_solubility: string;
  lipinski: string[];
  likely_toxic: boolean;
  toxicophores_detected: string[];
  validated_smiles: string;
  molecule_image_base64: string;
};

function App() {
  const [smiles, setSmiles] = useState('');
  const [result, setResult] = useState<PredictionResult | null>(null);
  const [error, setError] = useState('');

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    try {
      const response = await fetch(`http://localhost:8000/predict?smiles=${encodeURIComponent(smiles)}`);
      if (!response.ok) {
        throw new Error('Prediction failed');
      }
      const data = await response.json();
      setResult(data);
      setError('');
    } catch (err) {
      setError('Error fetching prediction.');
      setResult(null);
    }
  };

  return (
    <div style={{ padding: '2rem', fontFamily: 'Arial, sans-serif' }}>
      <h1>Molecular Property Predictor</h1>
      <form onSubmit={handleSubmit}>
        <input
          type="text"
          value={smiles}
          onChange={(e) => setSmiles(e.target.value)}
          placeholder="Enter SMILES string"
          style={{ width: '300px', padding: '0.5rem' }}
        />
        <button type="submit" style={{ marginLeft: '1rem', padding: '0.5rem 1rem' }}>
          Predict
        </button>
      </form>

      {error && <p style={{ color: 'red' }}>{error}</p>}

      {result && (
        <div style={{ marginTop: '2rem' }}>
          <h3>Results:</h3>
          <p><strong>Molecular Weight:</strong> {result["molecular weight"]}</p>
          <p><strong>LogP:</strong> {result.logP}</p>
          <p><strong>H-bond Donors:</strong> {result.H_bond_donors}</p>
          <p><strong>H-bond Acceptors:</strong> {result.H_bond_acceptors}</p>
          <p><strong>TPSA:</strong> {result.TPSA}</p>
          <p><strong>Rotatable Bonds:</strong> {result.rotatable_bonds}</p>
          <p><strong>Aromatic Rings:</strong> {result.aromatic_rings}</p>
          <p><strong>ESOL Solubility:</strong> {result.esol_solubility}</p>

          <p><strong>Lipinski Evaluation:</strong></p>
          <ul>
            <li>Violations: {result.lipinski.violations}</li>
          </ul>

          <p><strong>Likely Toxic?</strong> {result.likely_toxic ? 'Yes' : 'No'}</p>

          {result.toxicophores_detected.length > 0 && (
            <>
              <p><strong>Toxicophores Detected:</strong></p>
              <ul>
                {result.toxicophores_detected.map((t, idx) => (
                  <li key={idx}>{t}</li>
                ))}
              </ul>
            </>
          )}

          <p><strong>Validated SMILES:</strong> {result.validated_smiles}</p>

          {result.molecule_image_base64 && (
            <img
              src={`data:image/png;base64,${result.molecule_image_base64}`}
              alt="Molecule Structure"
              style={{ maxWidth: '300px', marginTop: '1rem', border: '1px solid #ccc' }}
            />
          )}
        </div>
      )}
    </div>
  );
}

export default App;