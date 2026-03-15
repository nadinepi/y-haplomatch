import React, { useState } from 'react';
import { MapContainer, TileLayer, CircleMarker, Popup, Marker } from 'react-leaflet';
import L from 'leaflet';
import 'leaflet/dist/leaflet.css'


// this just makes a nice color gradient from red to green for the dots
function colorFromSharedMutations(val, min, max) {
  if (max === min) return '#4ade80'; // if all values are the same, just use green
  const t = (val - min) / (max - min);
  // these are the main colors we blend between
  const red = [239, 68, 68];      // kinda bright red
  const orange = [251, 146, 60];  // orange vibes
  const yellow = [250, 204, 21];  // yellow, like a highlighter
  const green = [34, 197, 94];    // chill green

  let color;
  if (t <= 0.33) {
    // blend from red to orange
    const localT = t / 0.33;
    color = red.map((c, i) => Math.round(c * (1 - localT) + orange[i] * localT));
  } else if (t <= 0.66) {
    // blend from orange to yellow
    const localT = (t - 0.33) / 0.33;
    color = orange.map((c, i) => Math.round(c * (1 - localT) + yellow[i] * localT));
  } else {
    // blend from yellow to green
    const localT = (t - 0.66) / 0.34;
    color = yellow.map((c, i) => Math.round(c * (1 - localT) + green[i] * localT));
  }
  return `rgb(${color[0]},${color[1]},${color[2]})`;
}

function formatBP(bp) {
  if (bp == null) return '—';
  const year = Math.round(bp);
  if (year > 1950) return `${year - 1950} BCE`;
  return `${1950 - year} CE`;
}

function ResultsMap({ results }) {

  // only keep results that actually have some shared mutations and coordinates
  const filtered = results.filter(r => {
    const shared = Array.isArray(r.shared_mutations) ? r.shared_mutations.length : r.snps_matched;
    return shared > 0 && r.snps_compared > 0 && r.lat != null && r.lon != null;
  });
  if (filtered.length === 0) return null;
  // get the average lat/lon so we can center the map
  const avgLat = filtered.reduce((s, r) => s + r.lat, 0) / filtered.length;
  const avgLon = filtered.reduce((s, r) => s + r.lon, 0) / filtered.length;

  // figure out the min and max shared mutations for the color scale
  const sharedCounts = filtered.map(r => Array.isArray(r.shared_mutations) ? r.shared_mutations.length : r.snps_matched);
  const minShared = Math.min(...sharedCounts);
  const maxShared = Math.max(...sharedCounts);

  // group samples by their lat/lon so we can stack them
  const groupKey = (r) => `${r.lat},${r.lon}`;
  const grouped = {};
  for (const r of filtered) {
    const key = groupKey(r);
    if (!grouped[key]) grouped[key] = [];
    grouped[key].push(r);
  }

  // this keeps track of which popups are open
  const [expanded, setExpanded] = useState({});
  const MAX_SHOWN = 5;

  return (
    <div className="map-container">
      <h2 style={{ margin: 0 }}>Geographic Distribution</h2>
      <div className="map-wrapper" style={{ position: 'relative' }}>
        <MapContainer
          center={[avgLat, avgLon]}
          zoom={3}
          style={{ height: '100%', width: '100%' }}
        >
          <TileLayer
            attribution='&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a>'
            url="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png"
          />
          {Object.entries(grouped).map(([key, group]) => {
            // just grab the first sample for color and popup stuff
            const r = group[0];
            const isExpanded = expanded[key];
            const shown = isExpanded ? group : group.slice(0, MAX_SHOWN);
            // if there's more than one sample at a spot, show a number in a circle
            if (group.length > 1) {
              // use the highest shared mutations in the group for the color
              const maxGroupShared = Math.max(...group.map(ind => Array.isArray(ind.shared_mutations) ? ind.shared_mutations.length : ind.snps_matched));
              const borderColor = colorFromSharedMutations(maxGroupShared, minShared, maxShared);
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
              `;
              const icon = L.divIcon({
                html: iconHtml,
                className: '',
                iconSize: [28, 28],
                iconAnchor: [14, 14],
                popupAnchor: [0, -14],
              });
              return (
                <Marker key={key} position={[r.lat, r.lon]} icon={icon}>
                  <Popup>
                    <div style={{ fontWeight: 600, marginBottom: 4 }}>
                      {`${group.length} individuals at this site:`}
                    </div>
                    <ul style={{ margin: 0, paddingLeft: 18, maxHeight: 180, overflowY: isExpanded ? 'auto' : 'hidden' }}>
                      {shown.map((ind) => (
                        <li key={ind.id}>
                          <strong>{ind.id}</strong> — Haplogroup: {ind.y_haplogroup_clean}
                          {ind.country && <> — {ind.country}</>}
                          {ind.date_mean != null && <> — ~{formatBP(ind.date_mean)}</>}
                        </li>
                      ))}
                    </ul>
                    {!isExpanded && group.length > MAX_SHOWN && (
                      <button
                        style={{ marginTop: 6, fontSize: 13, color: '#2563eb', background: 'none', border: 'none', cursor: 'pointer', textDecoration: 'underline' }}
                        onClick={e => {
                          e.stopPropagation();
                          setExpanded((prev) => ({ ...prev, [key]: true }));
                        }}
                      >
                        Show all ({group.length})
                      </button>
                    )}
                    {isExpanded && group.length > MAX_SHOWN && (
                      <button
                        style={{ marginTop: 6, fontSize: 13, color: '#2563eb', background: 'none', border: 'none', cursor: 'pointer', textDecoration: 'underline' }}
                        onClick={e => {
                          e.stopPropagation();
                          setExpanded((prev) => ({ ...prev, [key]: false }));
                        }}
                      >
                        Show less
                      </button>
                    )}
                  </Popup>
                </Marker>
              );
            } else {
              // just color the dot by how many shared mutations
              const shared = Array.isArray(r.shared_mutations) ? r.shared_mutations.length : r.snps_matched;
              return (
                <CircleMarker
                  key={key}
                  center={[r.lat, r.lon]}
                  radius={7}
                  fillColor={colorFromSharedMutations(shared, minShared, maxShared)}
                  color="#fff"
                  weight={1}
                  fillOpacity={0.85}
                >
                  <Popup>
                    <div style={{ fontWeight: 600, marginBottom: 4 }}>{r.id}</div>
                    Haplogroup: {r.y_haplogroup_clean}<br />
                    {r.country && <>Country: {r.country}<br /></>}
                    {r.date_mean != null && <>Date: ~{formatBP(r.date_mean)}</>}
                  </Popup>
                </CircleMarker>
              );
            }
          })}
        </MapContainer>
        {/* floating legend overlay */}
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
  );
}

export default ResultsMap
