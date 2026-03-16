import os
from flask import Flask, request, jsonify
from flask_cors import CORS

from matching import match_with_broadening, resolve_haplogroup_simple
from parse_user import parse_user_file

app = Flask(__name__)
CORS(app)

DB_PATH = os.path.join(os.path.dirname(__file__), '..', 'data', 'ydna.db')


@app.route('/api/match', methods=['POST'])
def match():
    print("Received match request")
    # get haplogroup input from user
    user_input = request.form.get('haplogroup', '').strip()
    if not user_input:
        return jsonify({'error': 'Missing haplogroup'}), 400

    # resolve haplogroup using SNP reference
    resolved, snp, fallback = resolve_haplogroup_simple(DB_PATH, user_input)
    if not resolved:
        return jsonify({'error': f'Could not resolve haplogroup: {user_input}'}), 404

    # get uploaded genotype file
    uploaded = request.files.get('file')
    if not uploaded:
        return jsonify({'error': 'Missing genotype file'}), 400

    file_content = uploaded.read().decode('utf-8', errors='replace')
    user_genotypes, user_labels = parse_user_file(file_content, DB_PATH)
    if not user_genotypes:
        return jsonify({
            'error': 'Sorry, this tool did not recognize any SNPs from your file in our Y SNP database.'
        }), 400

    # limit for number of results
    try:
        limit = int(request.form.get('limit', 1000))
    except (TypeError, ValueError):
        limit = 1000
    limit = min(max(limit, 1), 1000)

    # broaden if needed
    results, final_haplo, levels_broadened = match_with_broadening(
        DB_PATH,
        resolved,
        user_genotypes,
        user_labels=user_labels,
        limit=limit,
        min_snps=1,  # allow matches with as few as 1 SNP
    )

    # resolution notes
    notes = []
    if fallback:
        notes.append(f'Input "{user_input}" interpreted as haplogroup "{resolved}".')
    if levels_broadened > 0:
        notes.append(f'No matches at "{resolved}", broadened {levels_broadened} level(s) to "{final_haplo}".')
    if not results:
        notes.append(f'No matches found even after broadening to "{final_haplo}". Try a broader haplogroup input.')

    response = {
        'resolved_haplogroup': resolved,
        'resolved_snp': snp,
        'user_snps_parsed': len(user_genotypes),
        'count': len(results),
        'results': results,
        'broadened_to': final_haplo if levels_broadened > 0 else None,
        'levels_broadened': levels_broadened,
        'resolution_note': ' '.join(notes) if notes else None
    }

    return jsonify(response)


@app.route('/api/individual/<individual_id>')
def get_individual(individual_id):
    # get full details for one ancient individual, including genotypes
    import sqlite3
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    cur.execute("SELECT * FROM individuals WHERE id = ?", (individual_id,))
    row = cur.fetchone()
    if not row:
        conn.close()
        return jsonify({'error': 'Individual not found'}), 404

    individual = dict(row)

    # add their genotypes
    cur.execute(
        "SELECT snp_name, position, value, is_derived FROM genotypes WHERE individual_id = ?",
        (individual_id,),
    )
    individual['genotypes'] = {
        r['snp_name']: {
            'position': r['position'],
            'value': r['value'],
            'is_derived': r['is_derived'],
        }
        for r in cur
    }
    conn.close()
    return jsonify(individual)


if __name__ == '__main__':
    port = int(os.environ.get('PORT', '5000'))
    app.run(debug=True, host='0.0.0.0', port=port)
