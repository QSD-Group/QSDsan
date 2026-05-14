# -*- coding: utf-8 -*-
"""
build_wrrf_rst.py
Generate SVG blocks for wrrf_interactive.rst from WERF_WRRF_diagrams.html.

Usage:
  python build_wrrf_rst.py                 # regenerate B1 SVG only (default)
  python build_wrrf_rst.py B1 B2           # regenerate specific config SVGs
  python build_wrrf_rst.py all             # regenerate all 18 SVGs
  python build_wrrf_rst.py --update-data   # inject wrrf_flow_data.json into const FD
  python build_wrrf_rst.py all --update-data  # both

Each config's <div class="wrrf-diagram" id="XX">...</div> block in the RST
is replaced with a freshly processed version from the source HTML.

Transformation pipeline (applied to every config):
  1. Remove Disinfection unit; extend Effluent line to FC/MBR output
  2. Remove smalllabel unit-name annotations
  3. Replace <style> block with canonical RST fonts (32/21/18 px)
  4. subtitle y: 60 → 78
  5. LIQUID TREATMENT y: 110 → 190
  6. Raw WW → RWW
  7. FC → SC (Final → Secondary Clarifier)
  8. Remove note annotations (process-specific bullet points)
  9. Remove bioreactor description section labels (keep LIQUID/SOLIDS TREATMENT only)
 10. Replace legend with canonical translate(36,630) version

Config-specific post-processing:
  G1/G2/G3: remove PE split label; shorten 80%/20% flow labels; add External C label
  H1:       add Coagulant + External C labels; move MLR y=140→170
  N1:       add External C label; move MLR y=140→170
  N2:       add External C label; move MLR y=140→170 (unique context to avoid I1/I2/I3)
"""

import re, os, sys, json, textwrap

HERE    = os.path.dirname(os.path.abspath(__file__))
REPO    = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
IN_HTML = os.path.join(REPO, 'WRRF_configs', 'WERF_WRRF_diagrams.html')
FD_JSON = os.path.join(HERE, 'wrrf_flow_data.json')
OUT_RST = os.path.join(HERE, 'wrrf_interactive.rst')

ALL_CONFIGS = ['B1', 'B2', 'B3', 'C1', 'C2', 'C3', 'E2', 'E2P',
               'F1', 'G1', 'G2', 'G3', 'H1', 'I1', 'I2', 'I3', 'N1', 'N2']

# ── canonical RST <style> block ──────────────────────────────────────────────
RST_CSS = """\
<style>
    .bg          { fill: #fbfbf8; }
    .title       { font: 700 32px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #1a1a1a; }
    .subtitle    { font: 400 21px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #555; }
    .smalllabel  { font: 500 21px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #444; text-anchor: middle; }
    .streamlabel { font: 500 21px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #444; }
    .zonelabel   { font: 700 21px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #fff; text-anchor: middle; }
    .zonelabel-d { font: 700 21px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #1a1a1a; text-anchor: middle; }
    .section     { font: 700 18px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #888; text-anchor: middle; letter-spacing:2px; }
    .note        { font: 400 21px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #555; }
    .dim         { font: 400 18px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #888; text-anchor: middle; }
    .unmodeled   { fill: #ececec; stroke: #aaa; stroke-width:1; stroke-dasharray:3,2; }

    .anaerobic   { fill: #6b4a8a; stroke: #4a3463; stroke-width:1; }
    .anoxic      { fill: #2f7d6e; stroke: #1f5b50; stroke-width:1; }
    .aerobic     { fill: #4a90c0; stroke: #2e6b95; stroke-width:1; }
    .clarifier   { fill: #e0d4b8; stroke: #8c7a55; stroke-width:1.2; }
    .digester-an { fill: #6b4a8a; stroke: #4a3463; stroke-width:1.2; }
    .digester-ae { fill: #4a90c0; stroke: #2e6b95; stroke-width:1.2; }
    .thickener   { fill: #c8bb9c; stroke: #8c7a55; stroke-width:1; }
    .mbr         { fill: #2e6b95; stroke: #1a4e75; stroke-width:1; }
    .chem        { fill: #d97a4a; stroke: #a85a30; stroke-width:1; }
    .dw          { fill: #b8a878; stroke: #7a6a40; stroke-width:1; }

    .mainflow    { stroke: #1a1a1a; stroke-width:1.6; fill:none; }
    .recycle     { stroke: #c2272d; stroke-width:1.4; fill:none; stroke-dasharray: 5,3; }
    .sludge      { stroke: #7a5a30; stroke-width:1.5; fill:none; }
    .gas         { stroke: #6b4a8a; stroke-width:1.4; fill:none; stroke-dasharray: 3,2; }
    .dose        { stroke: #d97a4a; stroke-width:1.4; fill:none; }
  </style>"""

# ── canonical legend at translate(36,630) ────────────────────────────────────
RST_LEGEND = (
    '<g transform="translate(36,630)">\n'
    '   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/>'
        '<text x="19" y="9" class="streamlabel">Anaerobic</text>\n'
    '   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/>'
        '<text x="168" y="9" class="streamlabel">Anoxic</text>\n'
    '   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/>'
        '<text x="279" y="9" class="streamlabel">Aerobic</text>\n'
    '   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/>'
        '<text x="390" y="9" class="streamlabel">Clarifier</text>\n'
    '   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/>'
        '<text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>\n'
    '   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/>'
        '<text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>\n'
    '   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/>'
        '<text x="19" y="39" class="streamlabel">Thickener</text>\n'
    '   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/>'
        '<text x="168" y="39" class="streamlabel">Dewatering</text>\n'
    '   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/>'
        '<text x="325" y="39" class="streamlabel">Membrane (MBR)</text>\n'
    '   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/>'
        '<text x="519" y="39" class="streamlabel">Chemical dose</text>\n'
    '   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/>'
        '<text x="30" y="74" class="streamlabel">Main flow</text>\n'
    '   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/>'
        '<text x="175" y="74" class="streamlabel">Recycle</text>\n'
    '   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/>'
        '<text x="303" y="74" class="streamlabel">Sludge</text>\n'
    '   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/>'
        '<text x="425" y="74" class="streamlabel">Biogas</text>\n'
    '   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/>'
        '<text x="547" y="74" class="streamlabel">Chemical addition</text>\n'
    '   </g>'
)

# ── dis-removal (same logic as build_interactive_html.py) ────────────────────
_DIS_PAT = re.compile(
    r'<line\s+x1="([\d.]+)"\s+y1="[\d.]+"\s+x2="[\d.]+"\s+y2="[\d.]+"\s+class="mainflow"[^/]*/>'
    r'<text\b[^>]*class="streamlabel"[^>]*>[^<]+</text>'
    r'<rect\b[^>]*class="unmodeled"[^>]*/>'
    r'<text\b[^>]*class="zonelabel-d">Dis\.</text>'
    r'<text\b[^>]*class="smalllabel">Disinfection</text>'
    r'<text\b[^>]*class="dim">\(not modeled\)</text>'
    r'(<line\s+x1="[\d.]+"\s+y1="[\d.]+"\s+x2="[\d.]+"\s+y2="[\d.]+"\s+class="mainflow"[^/]*/>)'
)

def _dis_sub(m):
    x1_fc    = m.group(1)
    effl_old = m.group(2)
    return re.sub(r'x1="[\d.]+"', f'x1="{x1_fc}"', effl_old, count=1)

# ── common transforms ────────────────────────────────────────────────────────
def apply_common(svg):
    # 1. remove Dis. unit
    svg = _DIS_PAT.sub(_dis_sub, svg)
    # 2. remove unit-name smalllabels (but not chemical-dose labels added later)
    svg = re.sub(r'<text\b[^>]*class="smalllabel"[^>]*>[^<]+</text>', '', svg)
    # 3. canonical CSS
    svg = re.sub(r'<style>.*?</style>', RST_CSS, svg, flags=re.DOTALL)
    # 4. subtitle y 60→78
    svg = re.sub(r'(<text x="36" y=")60(" class="subtitle">)', r'\g<1>78\2', svg)
    # 5. LIQUID TREATMENT y 110→190
    svg = svg.replace('" y="110" class="section">LIQUID TREATMENT</text>',
                      '" y="190" class="section">LIQUID TREATMENT</text>')
    # 6. Raw WW → RWW
    svg = svg.replace('>Raw WW</text>', '>RWW</text>')
    svg = svg.replace('>Raw WW (wastewater)</text>', '>RWW</text>')
    # 7. FC → SC
    svg = svg.replace('class="zonelabel-d">FC</text>', 'class="zonelabel-d">SC</text>')
    # 8. remove note annotations (process-specific bullet points that clutter the diagram)
    svg = re.sub(r'<text\b[^>]*class="note"[^>]*>[^<]+</text>', '', svg)
    # 9. remove bioreactor description section labels; keep only structural headers
    def _keep_section(m):
        return m.group(0) if m.group(1).strip() in ('LIQUID TREATMENT', 'SOLIDS TREATMENT') else ''
    svg = re.sub(r'<text\b[^>]*class="section"[^>]*>([^<]+)</text>', _keep_section, svg)
    # 10. replace legend (any translate(36,N) group)
    svg = re.sub(r'<g transform="translate\(36,\d+\)">.*?</g>', RST_LEGEND,
                 svg, flags=re.DOTALL)
    return svg

# ── config-specific post-processing ─────────────────────────────────────────
def _chem_box_label(x, y_center, box_text):
    """White zonelabel text centered inside an orange chem box."""
    return f'<text x="{x}" y="{y_center}" class="zonelabel" dominant-baseline="central">{box_text}</text>'

def _resize_chem_rect(svg, x_old, y_old, x_new, y_new, dose_x, dose_y1_old, dose_y1_new, dose_y2):
    """Enlarge a chem pill rect and update the dose line y1 to match new box bottom."""
    svg = svg.replace(
        f'<rect x="{x_old}" y="{y_old}" width="72" height="18" class="chem" rx="9"/>',
        f'<rect x="{x_new}" y="{y_new}" width="120" height="30" class="chem" rx="15"/>')
    svg = svg.replace(
        f'<line x1="{dose_x}" y1="{dose_y1_old}" x2="{dose_x}" y2="{dose_y2}" class="dose" marker-end="url(#arr-dose)"/>',
        f'<line x1="{dose_x}" y1="{dose_y1_new}" x2="{dose_x}" y2="{dose_y2}" class="dose" marker-end="url(#arr-dose)"/>')
    return svg

def post_g_series(svg):
    """G1, G2, G3: remove PE split; shorten 80/20% labels; label chem box."""
    svg = svg.replace(
        '<text x="196" y="224.0" class="streamlabel" text-anchor="end">PE split</text>', '')
    svg = svg.replace('>80 % PE → A2</text>', '>80%</text>')
    svg = svg.replace('>20 % PE → A3 (step feed)</text>', '>20%</text>')
    # chem rect: old x=202 y=134 w=72 h=18 → new x=178 y=128 w=120 h=30; dose y1: 152→158
    svg = _resize_chem_rect(svg, '202.0', '134', '178.0', '128', '238.0', '152', '158', '204')
    # add inside-box label (center y=128+15=143)
    dose_new = '<line x1="238.0" y1="158" x2="238.0" y2="204" class="dose" marker-end="url(#arr-dose)"/>'
    svg = svg.replace(dose_new, dose_new + _chem_box_label('238', '143', 'External C'))
    return svg

def post_h1(svg):
    """H1: label Coagulant + External C boxes; move MLR y=140→170."""
    # both chem rects at y=119 w=72 h=18 → y=113 w=120 h=30; dose y1: 137→143
    svg = _resize_chem_rect(svg, '53.0',  '119', '29.0',  '113', '89.0',  '137', '143', '232.0')
    svg = _resize_chem_rect(svg, '202.0', '119', '178.0', '113', '238.0', '137', '143', '204')
    coag_dose = '<line x1="89.0" y1="143" x2="89.0" y2="232.0" class="dose" marker-end="url(#arr-dose)"/>'
    extc_dose = '<line x1="238.0" y1="143" x2="238.0" y2="204" class="dose" marker-end="url(#arr-dose)"/>'
    svg = svg.replace(coag_dose, coag_dose + _chem_box_label('89',  '128', 'Coagulant'))
    svg = svg.replace(extc_dose, extc_dose + _chem_box_label('238', '128', 'External C'))
    svg = svg.replace(
        '<text x="348.0" y="140" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text>',
        '<text x="348.0" y="170" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text>')
    return svg

def post_n1(svg):
    """N1: label External C box; move MLR y=140→170."""
    svg = _resize_chem_rect(svg, '202.0', '119', '178.0', '113', '238.0', '137', '143', '204')
    extc_dose = '<line x1="238.0" y1="143" x2="238.0" y2="204" class="dose" marker-end="url(#arr-dose)"/>'
    svg = svg.replace(extc_dose, extc_dose + _chem_box_label('238', '128', 'External C'))
    svg = svg.replace(
        '<text x="350.0" y="140" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text>',
        '<text x="350.0" y="170" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text>')
    return svg

def post_n2(svg):
    """N2: label External C box; move MLR y=140→170 (unique context avoids I1/I2/I3)."""
    svg = _resize_chem_rect(svg, '170.0', '119', '146.0', '113', '206.0', '137', '143', '204')
    extc_dose = '<line x1="206.0" y1="143" x2="206.0" y2="204" class="dose" marker-end="url(#arr-dose)"/>'
    extc_lbl  = _chem_box_label('206', '128', 'External C')
    svg = svg.replace(extc_dose, extc_dose + extc_lbl)
    # use dose + label as unique context to avoid touching I1/I2/I3
    n2_ctx_old = (
        extc_dose + extc_lbl +
        '<path d="M 388.0 205 L 388.0 145 L 248.0 145 L 248.0 205" class="recycle" marker-end="url(#arr-recycle)"/>'
        '<text x="318.0" y="140" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text>'
    )
    n2_ctx_new = (
        extc_dose + extc_lbl +
        '<path d="M 388.0 205 L 388.0 145 L 248.0 145 L 248.0 205" class="recycle" marker-end="url(#arr-recycle)"/>'
        '<text x="318.0" y="170" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text>'
    )
    svg = svg.replace(n2_ctx_old, n2_ctx_new)
    return svg

_POST = {
    'G1': post_g_series, 'G2': post_g_series, 'G3': post_g_series,
    'H1': post_h1,
    'N1': post_n1,
    'N2': post_n2,
}

# ── extract one config's SVG from source HTML ────────────────────────────────
def extract_svg(html, cfg_id):
    """Return the raw SVG string (opening tag through </svg>) for cfg_id."""
    marker = f'id="{cfg_id}">'
    start  = html.find(marker)
    if start < 0:
        raise ValueError(f'Config {cfg_id!r} not found in source HTML')
    svg_start = html.find('<svg', start)
    svg_end   = html.find('</svg>', svg_start) + 6
    return html[svg_start:svg_end]

# ── process one config ────────────────────────────────────────────────────────
def process_config(cfg_id, html):
    svg = extract_svg(html, cfg_id)
    svg = apply_common(svg)
    if cfg_id in _POST:
        svg = _POST[cfg_id](svg)
    return svg

# ── build the RST div block ───────────────────────────────────────────────────
def make_rst_block(cfg_id, svg):
    """Wrap SVG in the wrrf-diagram div and indent for RST raw:: html."""
    block = f'<div class="wrrf-diagram" id="{cfg_id}">\n{svg}\n</div>'
    # add 3-space RST indentation to every line
    return textwrap.indent(block, '   ')

# ── replace one config block in the RST string ───────────────────────────────
def replace_in_rst(rst, cfg_id, new_block):
    """
    Find and replace the <div class="wrrf-diagram" id="cfg_id">...</div> block.
    B1 is special: its closing </div> is followed by the tooltip/JS block in the
    same raw:: html section, so we search for </div> then the next blank line.
    For B2-N2 each block ends with </div> immediately before </svg> in the next
    block or the references section.
    """
    open_tag = f'<div class="wrrf-diagram" id="{cfg_id}">'
    start = rst.find('   ' + open_tag)
    if start < 0:
        raise ValueError(f'Could not find div block for {cfg_id!r} in RST')

    # find the </div> that closes this specific diagram div
    # (there is exactly one </div> per diagram block, right after </svg>)
    close_tag = '   </div>'
    end = rst.find(close_tag, start)
    if end < 0:
        raise ValueError(f'Could not find closing </div> for {cfg_id!r}')
    end += len(close_tag)

    return rst[:start] + new_block + rst[end:]

# ── flow data injection ───────────────────────────────────────────────────────
def inject_fd(rst):
    """Replace the const FD = {...}; line in the RST with fresh data from FD_JSON."""
    if not os.path.exists(FD_JSON):
        print(f'[FD] {FD_JSON} not found, skipping')
        return rst
    with open(FD_JSON, encoding='utf-8') as f:
        data = json.load(f)
    fd_str = json.dumps(data, separators=(',', ':'))
    new_line = f'   const FD = {fd_str};'
    rst, n = re.subn(r'   const FD = \{.*?\};', new_line, rst)
    if n == 0:
        print('[FD] WARNING: const FD line not found in RST — not updated')
    else:
        print(f'[FD] injected flow data for {len(data)} configs')
    return rst

# ── main ─────────────────────────────────────────────────────────────────────
def main():
    args = sys.argv[1:]
    update_data = '--update-data' in args
    targets = [a for a in args if a != '--update-data'] or ['B1']
    if targets == ['all']:
        targets = ALL_CONFIGS

    with open(OUT_RST, encoding='utf-8') as f:
        rst = f.read()

    if update_data:
        rst = inject_fd(rst)

    if targets:
        with open(IN_HTML, encoding='utf-8') as f:
            source_html = f.read()
        for cfg_id in targets:
            if cfg_id not in ALL_CONFIGS:
                print(f'Unknown config {cfg_id!r}, skipping')
                continue
            print(f'Processing {cfg_id}...', end=' ')
            svg   = process_config(cfg_id, source_html)
            block = make_rst_block(cfg_id, svg)
            rst   = replace_in_rst(rst, cfg_id, block)
            print('done')

    with open(OUT_RST, 'w', encoding='utf-8') as f:
        f.write(rst)
    print(f'Written -> {OUT_RST}')

if __name__ == '__main__':
    main()
