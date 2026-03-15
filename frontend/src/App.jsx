import { useState } from 'react'
import ResultsTable from './components/ResultsTable'
import ResultsMap from './components/ResultsMap'
import './App.css'

function App() {
  const [haplogroup, setHaplogroup] = useState('')
  const [file, setFile] = useState(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState(null)
  const [data, setData] = useState(null)

  const handleSubmit = async (e) => {
    e.preventDefault()
    if (!haplogroup.trim() || !file) return

    setLoading(true)
    setError(null)
    setData(null)

    const formData = new FormData()
    formData.append('haplogroup', haplogroup.trim())
    formData.append('file', file)
    formData.append('limit', '1000')

    try {
      const res = await fetch('/api/match', {
        method: 'POST',
        body: formData,
      })
      const json = await res.json()
      if (!res.ok) {
        setError(json.error || 'Request failed')
      } else {
        setData(json)
      }
    } catch (err) {
      setError('Could not connect to backend. Is the server running?')
    } finally {
      setLoading(false)
    }
  }

  return (
    <div className="app">
      <h1>Y-HaploMatch</h1>
      <p className="subtitle">
        Compare your Y chromosome DNA against ancient individuals from the AADR database
      </p>

      <form className="upload-form" onSubmit={handleSubmit}>
        <div className="form-row">
          <div className="form-group">
            <label htmlFor="haplogroup">Y Haplogroup (e.g. R-M269, R1a2)</label>
            <input
              id="haplogroup"
              type="text"
              value={haplogroup}
              onChange={(e) => setHaplogroup(e.target.value)}
              placeholder="R-M269"
            />
          </div>
          <div className="form-group">
            <label htmlFor="file">Genotype file (PLINK .raw or tab-separated)</label>
            <input
              id="file"
              type="file"
              accept=".txt,.raw,.csv,.tsv"
              onChange={(e) => setFile(e.target.files[0])}
            />
          </div>
          <button
            type="submit"
            className="btn btn-primary"
            disabled={loading || !haplogroup.trim() || !file}
          >
            {loading ? 'Matching...' : 'Find Matches'}
          </button>
        </div>
      </form>

      {loading && (
        <div className="loading">
          <span className="spinner" />
          Comparing genotypes...
        </div>
      )}

      {error && <div className="error">{error}</div>}

      {data && (
        <>
          {data.resolution_note && (
            <div className="warning">{data.resolution_note}</div>
          )}

          <div className="results-summary">
            <div className="summary-card">
              Resolved haplogroup: <strong>{data.resolved_haplogroup}</strong>
            </div>
            {data.resolved_snp && (
              <div className="summary-card">
                Resolved by SNP: <strong>{data.resolved_snp}</strong>
              </div>
            )}
            <div className="summary-card">
              Y SNPs found from your file: <strong>{data.user_snps_parsed}</strong>
            </div>
            <div className="summary-card">
              Matches found: <strong>{
                data.results
                  ? data.results.filter(r => {
                      const shared = Array.isArray(r.shared_mutations) ? r.shared_mutations.length : r.snps_matched;
                      return shared > 0 && r.snps_compared > 0;
                    }).length
                  : 0
              }</strong>
            </div>
          </div>

          <ResultsMap results={data.results} />
          <ResultsTable results={data.results} />
        </>
      )}
    </div>
  )
}

export default App
