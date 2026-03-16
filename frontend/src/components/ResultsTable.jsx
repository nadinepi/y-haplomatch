import { useState } from 'react'

const COLUMNS = [
  { key: 'id', label: 'Individual ID' },
  { key: 'y_haplogroup_clean', label: 'Haplogroup' },
  { key: 'country', label: 'Country' },
  { key: 'full_date_range', label: 'Date Range' },
  { key: 'snps_compared', label: 'SNPs Compared' },
  { key: 'snps_matched', label: 'Shared Mutations' },
]

// get a unique list of shared mutation labels for this result
const getSharedList = (result) => {
  if (!Array.isArray(result.shared_mutations)) return []
  const seen = new Set()
  const cleaned = []
  for (const mut of result.shared_mutations) {
    const label = typeof mut === 'string' ? mut.trim() : String(mut)
    if (!label || seen.has(label)) continue
    seen.add(label)
    cleaned.push(label)
  }
  return cleaned
}


function ResultsTable({ results }) {
  const [sortKey, setSortKey] = useState('snps_matched')
  const [sortAsc, setSortAsc] = useState(false)
  const [expandedRowIds, setExpandedRowIds] = useState(new Set())
  const [page, setPage] = useState(1)
  const PAGE_SIZE = 50

  // only keep results that actually have some shared mutations and were compared
  const filtered = results.filter(r => {
    // shared_mutations can be an array or a number (snps_matched)
    const shared = Array.isArray(r.shared_mutations) ? r.shared_mutations.length : r.snps_matched;
    return shared > 0 && r.snps_compared > 0;
  });

  if (!filtered || filtered.length === 0) {
    return <div className="info">No matching individuals found.</div>;
  }

  const sorted = [...filtered].sort((a, b) => {
    let va = a[sortKey];
    let vb = b[sortKey];
    if (va == null) va = sortAsc ? Infinity : -Infinity;
    if (vb == null) vb = sortAsc ? Infinity : -Infinity;
    if (typeof va === 'string') return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
    return sortAsc ? va - vb : vb - va;
  });

  // pagination stuff so we don't show a million rows at once
  const totalPages = Math.ceil(sorted.length / PAGE_SIZE);
  const paged = sorted.slice((page - 1) * PAGE_SIZE, page * PAGE_SIZE);

  // handle sorting when you click a column
  const handleSort = (key) => {
    if (key === sortKey) {
      setSortAsc(!sortAsc)
    } else {
      setSortKey(key)
      setSortAsc(false)
    }
    setPage(1)
  }

  // go to a different page
  const handlePageChange = (newPage) => {
    setPage(newPage)
  }

  return (
    <div className="results-table-wrapper">
      <table className="results-table">
        <thead>
          <tr>
            <th style={{ width: 60 }}>#</th>
            {COLUMNS.map((col) => (
              <th
                key={col.key}
                onClick={col.sortable !== false ? () => handleSort(col.key) : undefined}
                style={col.sortable === false ? { cursor: 'default' } : {}}
              >
                {col.label}
                {sortKey === col.key && (sortAsc ? ' ▲' : ' ▼')}
              </th>
            ))}
            <th style={{ width: 120, textAlign: 'center' }}>SNP List</th>
          </tr>
        </thead>
        <tbody>
          {paged.map((r, i) => (
            <>
              <tr key={r.id} className={r.low_overlap ? 'low-overlap-row' : ''} style={{ position: 'relative' }}>
                <td>{(page - 1) * PAGE_SIZE + i + 1}</td>
                <td>{r.id}</td>
                <td>{r.y_haplogroup_clean}</td>
                <td>{r.country || '—'}</td>
                <td>
                  {r.full_date_range || '—'}
                </td>
                <td>{r.snps_compared ?? '—'}</td>
                <td>{r.snps_matched ?? '—'}</td>
                <td style={{ textAlign: 'center', verticalAlign: 'middle' }}>
                  <button
                    onClick={() => setExpandedRowIds(prev => {
                      const newSet = new Set(prev)
                      if (newSet.has(r.id)) {
                        newSet.delete(r.id)
                      } else {
                        newSet.add(r.id)
                      }
                      return newSet
                    })}
                    aria-label={expandedRowIds.has(r.id) ? 'Hide SNP list' : 'Show SNP list'}
                    className="details-toggle"
                  >
                    {expandedRowIds.has(r.id) ? 'Hide' : 'Show'}
                  </button>
                </td>
              </tr>
              {expandedRowIds.has(r.id) && (
                <tr key={r.id + '-details'}>
                  <td colSpan={COLUMNS.length + 2} className="details-row-cell">
                    <div className="details-panel">
                      <div className="details-title">Shared mutation labels</div>
                      {getSharedList(r).length > 0 ? (
                        <div className="mutation-chip-list">
                          {getSharedList(r).map((mut) => {
                            let display = mut;
                            if (typeof mut === 'string') {
                              // Replace any leading RS/rs/Rs/rS with lowercase 'rs'
                              display = mut.replace(/^\s*rs/i, 'rs');
                            }
                            return (
                              <span key={mut} className="mutation-chip">{display}</span>
                            );
                          })}
                        </div>
                      ) : (
                        <span className="details-empty">No shared mutation labels available.</span>
                      )}
                    </div>
                  </td>
                </tr>
              )}
            </>
          ))}
        </tbody>
      </table>
      {/* pagination controls at the bottom */}
      {totalPages > 1 && (
        <div className="pagination-controls-modern" style={{
          marginTop: 24,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          gap: 16,
          fontFamily: 'Inter, Roboto, Arial, sans-serif',
        }}>
          <button
            onClick={() => handlePageChange(1)}
            disabled={page === 1}
            className="modern-btn"
            aria-label="First page"
          >⏮</button>
          <button
            onClick={() => handlePageChange(page - 1)}
            disabled={page === 1}
            className="modern-btn"
            aria-label="Previous page"
          >◀</button>
          <span style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
            <span style={{ fontWeight: 500, fontSize: 15 }}>Page</span>
            <input
              type="number"
              min={1}
              max={totalPages}
              value={page}
              onChange={e => {
                let val = Number(e.target.value)
                if (val >= 1 && val <= totalPages) handlePageChange(val)
              }}
              style={{
                width: 48,
                textAlign: 'center',
                borderRadius: 6,
                border: '1.5px solid #bdbdbd',
                padding: '4px 6px',
                fontSize: 15,
                boxShadow: '0 1px 2px rgba(0,0,0,0.04)',
                outline: 'none',
                transition: 'border 0.2s',
              }}
            />
            <span style={{ fontWeight: 500, fontSize: 15 }}>of {totalPages}</span>
          </span>
          <button
            onClick={() => handlePageChange(page + 1)}
            disabled={page === totalPages}
            className="modern-btn"
            aria-label="Next page"
          >▶</button>
          <button
            onClick={() => handlePageChange(totalPages)}
            disabled={page === totalPages}
            className="modern-btn"
            aria-label="Last page"
          >⏭</button>
        </div>
      )}
    </div>
  )
}

export default ResultsTable
