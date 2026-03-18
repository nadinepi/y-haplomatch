import React, { useEffect, useMemo, useRef, useState } from 'react'
import { CircleMarker, MapContainer, Marker, Popup, TileLayer } from 'react-leaflet'
import L from 'leaflet'
import 'leaflet/dist/leaflet.css'

function colorFromSharedMutations(val, min, max) {
  if (max === min) return '#4ade80'
  const t = (val - min) / (max - min)
  const red = [239, 68, 68]
  const orange = [251, 146, 60]
  const yellow = [250, 204, 21]
  const green = [34, 197, 94]

  let color
  if (t <= 0.33) {
    const localT = t / 0.33
    color = red.map((c, i) => Math.round(c * (1 - localT) + orange[i] * localT))
  } else if (t <= 0.66) {
    const localT = (t - 0.33) / 0.33
    color = orange.map((c, i) => Math.round(c * (1 - localT) + yellow[i] * localT))
  } else {
    const localT = (t - 0.66) / 0.34
    color = yellow.map((c, i) => Math.round(c * (1 - localT) + green[i] * localT))
  }
  return `rgb(${color[0]},${color[1]},${color[2]})`
}

function formatBP(bp) {
  if (bp == null) return '—'
  const year = Math.round(bp)
  if (year > 1950) return `${year - 1950} BCE`
  return `${1950 - year} CE`
}

function ResultsMap({ results, selectedResultId, mapFocusTick, onSelectResult }) {
  const [expanded, setExpanded] = useState({})
  const mapRef = useRef(null)
  const markerRefs = useRef({})
  const lastHandledFocusTick = useRef(0)
  const MAX_SHOWN = 5

  const filtered = useMemo(() => {
    return results.filter((r) => {
      const shared = Array.isArray(r.shared_mutations) ? r.shared_mutations.length : r.snps_matched
      return shared > 0 && r.snps_compared > 0 && r.lat != null && r.lon != null
    })
  }, [results])

  const grouped = useMemo(() => {
    const next = {}
    for (const row of filtered) {
      const key = `${row.lat},${row.lon}`
      if (!next[key]) next[key] = []
      next[key].push(row)
    }
    return next
  }, [filtered])

  useEffect(() => {
    if (!selectedResultId) return
    if (mapFocusTick === lastHandledFocusTick.current) return

    const selected = filtered.find((r) => r.id === selectedResultId)
    if (!selected) return

    const key = `${selected.lat},${selected.lon}`
    setExpanded((prev) => ({ ...prev, [key]: true }))

    if (mapRef.current) {
      mapRef.current.setView([selected.lat, selected.lon], Math.max(mapRef.current.getZoom(), 5), {
        animate: false,
      })
    }

    const marker = markerRefs.current[key]
    if (marker && marker.openPopup) {
      marker.openPopup()
    }

    lastHandledFocusTick.current = mapFocusTick
  }, [filtered, selectedResultId, mapFocusTick])

  if (filtered.length === 0) return null

  const avgLat = filtered.reduce((s, r) => s + r.lat, 0) / filtered.length
  const avgLon = filtered.reduce((s, r) => s + r.lon, 0) / filtered.length

  const sharedCounts = filtered.map((r) => Array.isArray(r.shared_mutations) ? r.shared_mutations.length : r.snps_matched)
  const minShared = Math.min(...sharedCounts)
  const maxShared = Math.max(...sharedCounts)

  return (
    <div className="map-container">
      <h2 style={{ margin: 0 }}>Geographic Distribution</h2>
      <div className="map-wrapper" style={{ position: 'relative' }}>
        <MapContainer
          center={[avgLat, avgLon]}
          zoom={3}
          style={{ height: '100%', width: '100%' }}
          ref={mapRef}
        >
          <TileLayer
            attribution='&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a>'
            url="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png"
          />

          {Object.entries(grouped).map(([key, group]) => {
            const first = group[0]
            const isExpanded = expanded[key]
            const shown = isExpanded ? group : group.slice(0, MAX_SHOWN)
            const maxGroupShared = Math.max(...group.map((ind) => Array.isArray(ind.shared_mutations) ? ind.shared_mutations.length : ind.snps_matched))
            const borderColor = colorFromSharedMutations(maxGroupShared, minShared, maxShared)

            if (group.length > 1) {
              const iconHtml = `
                <div style="
                  display: flex;
                  align-items: center;
                  justify-content: center;
                  width: 28px;
                  height: 28px;
                  background: #fff;
                  color: #222;
                  font-size: 15px;
                  font-weight: 700;
                  border: 3px solid ${borderColor};
                  border-radius: 50%;
                  box-shadow: 0 1px 4px rgba(0,0,0,0.12);
                ">
                  ${group.length}
                </div>
              `

              const icon = L.divIcon({
                html: iconHtml,
                className: '',
                iconSize: [28, 28],
                iconAnchor: [14, 14],
                popupAnchor: [0, -14],
              })

              return (
                <Marker
                  key={key}
                  position={[first.lat, first.lon]}
                  icon={icon}
                  ref={(layer) => {
                    if (layer) markerRefs.current[key] = layer
                  }}
                >
                  <Popup>
                    <div className="map-popup">
                      <div className="map-popup-title">{group.length} individuals at this location</div>
                      <div className="map-popup-subtitle">
                        Click a sample to highlight it in the table.
                      </div>

                      <div className="map-popup-list">
                        {shown.map((ind) => (
                          <button
                            key={ind.id}
                            type="button"
                            className={`map-popup-card ${selectedResultId === ind.id ? 'map-popup-card-selected' : ''}`}
                            onClick={() => onSelectResult && onSelectResult(ind.id)}
                          >
                            <div className="map-popup-card-header">
                              <strong>{ind.id}</strong>
                              <span>{ind.y_haplogroup_clean || 'Unknown haplogroup'}</span>
                            </div>
                            <div className="map-popup-card-body">
                              <span>{ind.country || 'Unknown country'}</span>
                              <span>{ind.full_date_range || `~${formatBP(ind.date_mean)}`}</span>
                              <span>{ind.snps_matched} shared mutation{ind.snps_matched === 1 ? '' : 's'}</span>
                            </div>
                          </button>
                        ))}
                      </div>

                      {group.length > MAX_SHOWN && (
                        <button
                          type="button"
                          className="map-popup-link"
                          onClick={() => setExpanded((prev) => ({ ...prev, [key]: !prev[key] }))}
                        >
                          {isExpanded ? 'Show fewer' : `Show all ${group.length}`}
                        </button>
                      )}
                    </div>
                  </Popup>
                </Marker>
              )
            }

            const row = group[0]
            const shared = Array.isArray(row.shared_mutations) ? row.shared_mutations.length : row.snps_matched
            return (
              <CircleMarker
                key={key}
                center={[row.lat, row.lon]}
                radius={selectedResultId === row.id ? 10 : 7}
                fillColor={colorFromSharedMutations(shared, minShared, maxShared)}
                color={selectedResultId === row.id ? '#1d4ed8' : '#fff'}
                weight={selectedResultId === row.id ? 3 : 1}
                fillOpacity={0.85}
                ref={(layer) => {
                  if (layer) markerRefs.current[key] = layer
                }}
                eventHandlers={{
                  click: () => {
                    if (onSelectResult) {
                      onSelectResult(row.id)
                    }
                  },
                }}
              >
                <Popup>
                  <div className="map-popup">
                    <div className="map-popup-title">{row.id}</div>
                    <div className="map-popup-subtitle">
                      {row.lineage_relation === 'exact' ? 'Exact branch match' : row.lineage_relation === 'broader' ? 'Broader ancestor match' : row.lineage_relation === 'downstream' ? 'Downstream branch match' : 'Position-based match'}
                    </div>
                    <div className="map-popup-facts">
                      <div><strong>Haplogroup:</strong> {row.y_haplogroup_clean || 'Unknown'}</div>
                      <div><strong>Country:</strong> {row.country || 'Unknown'}</div>
                      <div><strong>Date:</strong> {row.full_date_range || `~${formatBP(row.date_mean)}`}</div>
                      <div><strong>Shared mutations:</strong> {row.snps_matched}</div>
                    </div>
                  </div>
                </Popup>
              </CircleMarker>
            )
          })}
        </MapContainer>

        <div style={{ position: 'absolute', top: 16, right: 16, background: 'rgba(255,255,255,0.92)', borderRadius: 10, boxShadow: '0 2px 8px rgba(0,0,0,0.10)', padding: '10px 16px 8px 16px', zIndex: 1000, minWidth: 180 }}>
          <div style={{ fontWeight: 500, fontSize: 14, marginBottom: 4 }}>Shared Mutations</div>
          <div style={{ width: 140, height: 14, background: 'linear-gradient(to right, #ef4444, #fb923c 33%, #facc15 66%, #22c55e)', borderRadius: 7, border: '1px solid #ccc', marginBottom: 2 }} />
          <div style={{ width: 140, display: 'flex', justifyContent: 'space-between', fontSize: 12 }}>
            <span>{minShared} (fewest)</span>
            <span>{maxShared} (most)</span>
          </div>
          <div style={{ fontSize: 12, color: '#444', marginTop: 6 }}>
            <span style={{ fontWeight: 500, fontSize: 13 }}>ⓘ</span> Number in circle = # of individuals
          </div>
        </div>
      </div>
    </div>
  )
}

export default ResultsMap
