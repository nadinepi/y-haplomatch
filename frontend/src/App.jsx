import { useState } from 'react'
import ResultsTable from './components/ResultsTable'
import ResultsMap from './components/ResultsMap'
import './App.css'

function countVisibleMatches(results) {
  if (!results) return 0

  return results.filter((row) => {
    const shared = Array.isArray(row.shared_mutations) ? row.shared_mutations.length : row.snps_matched
    return shared > 0 && row.snps_compared > 0
  }).length
}

function App() {
  const [haplogroup, setHaplogroup] = useState('')
  const [file, setFile] = useState(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState(null)
  const [data, setData] = useState(null)
  const [selectedResultId, setSelectedResultId] = useState(null)
  const [mapFocusTick, setMapFocusTick] = useState(0)

  const selectResult = (resultId, focusMap = false) => {
    setSelectedResultId(resultId)
    if (focusMap) {
      setMapFocusTick((tick) => tick + 1)
    }
  }

  const runMatch = async (nextHaplogroup, nextFile = file) => {
    if (!nextHaplogroup.trim() || !nextFile) return
    setLoading(true)
    setError(null)
    setData(null)
    setSelectedResultId(null)
    setMapFocusTick(0)

    const formData = new FormData()
    formData.append('haplogroup', nextHaplogroup.trim())
    formData.append('file', nextFile)
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

  const handleSubmit = async (e) => {
    e.preventDefault()
    if (!haplogroup.trim() || !file) return
    await runMatch(haplogroup, file)
  }

  const handleSuggestedRerun = async () => {
    const next = data?.suggested_haplogroup?.haplogroup
    if (!next || !file) return
    setHaplogroup(next)
    await runMatch(next, file)
  }

  const loadDemoUser = async (runNow = false) => {
    try {
      const res = await fetch('/demo_i1_user.txt')
      const text = await res.text()
      const demoFile = new File([text], 'demo_i1_user.txt', { type: 'text/plain' })
      setHaplogroup('I1')
      setFile(demoFile)
      if (runNow) {
        await runMatch('I1', demoFile)
      }
    } catch (err) {
      setError('Could not load the built-in demo user.')
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
            <label htmlFor="file">Genotype file (txt, csv, tsv, or raw)</label>
            <input
              id="file"
              type="file"
              accept=".txt,.raw,.csv,.tsv"
              onChange={(e) => {
                setFile(e.target.files[0] || null)
              }}
            />
          </div>
          <button
            type="submit"
            className="btn btn-primary"
            disabled={loading || !haplogroup.trim() || !file}
          >
            {loading ? 'Matching...' : 'Find Matches'}
          </button>
          <button
            type="button"
            className="btn"
            onClick={() => loadDemoUser(true)}
            disabled={loading}
          >
            Run Demo User
          </button>
        </div>
      </form>

      {loading && (
        <div className="loading">
          <span className="spinner" />
          Matching your Y positions to ISOGG SNPs...
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
              Matches found: <strong>{countVisibleMatches(data.results)}</strong>
            </div>
          </div>

          {data.suggested_haplogroup && (
            <div className="warning" style={{ marginTop: 12 }}>
              Your uploaded SNPs may better support haplogroup <strong>{data.suggested_haplogroup.haplogroup}</strong>.
              {' '}
              <button
                type="button"
                className="btn btn-primary"
                onClick={handleSuggestedRerun}
                disabled={loading}
                style={{ marginLeft: 10 }}
              >
                Rerun with this haplogroup
              </button>
            </div>
          )}

          <ResultsMap
            results={data.results}
            selectedResultId={selectedResultId}
            mapFocusTick={mapFocusTick}
            onSelectResult={(resultId) => selectResult(resultId)}
          />
          <ResultsTable
            results={data.results}
            selectedResultId={selectedResultId}
            onSelectResult={(resultId) => selectResult(resultId, true)}
          />
        </>
      )}
    </div>
  )
}

export default App
