# -*- coding: utf-8 -*-
"""
Build WERF_WRRF_interactive.html:
  - Reads WERF_WRRF_diagrams.html  (existing nice hand-crafted SVGs)
  - Reads wrrf_flow_data.json       (simulation output from simulate_wrrf.py)
  - Processing steps:
      1. Remove Dis. unit from all 18 SVGs; extend effluent line to FC/MBR output
      2. Increase font sizes directly in SVG <style> blocks
      3. Inject CSS tooltip styles
      4. Inject JS: hit-area overlays, stream + unit hover highlights, tooltips
         (legend elements are excluded from hover effects)
  - Writes WRRF_configs/WERF_WRRF_interactive.html
"""

import json, os, re

HERE     = os.path.dirname(os.path.abspath(__file__))
ROOT     = os.path.dirname(HERE)
IN_HTML  = os.path.join(ROOT, 'WRRF_configs', 'WERF_WRRF_diagrams.html')
IN_JSON  = os.path.join(ROOT, 'WRRF_configs', 'wrrf_flow_data.json')
OUT_HTML = os.path.join(ROOT, 'WRRF_configs', 'WERF_WRRF_interactive.html')

# --------------------------------------------------------------------------- #
# 1. Load inputs
# --------------------------------------------------------------------------- #
with open(IN_JSON) as f:
    flow_data = json.load(f)

with open(IN_HTML, encoding='utf-8') as f:
    html = f.read()

# --------------------------------------------------------------------------- #
# 2. Remove Disinfection unit from all 18 SVGs
#    Pattern (all on one line per SVG):
#      [SE/Permeate line x1=A] [SE/Permeate label] [Dis.rect] [Dis.labels×3]
#      [Effluent line x1=B]
#    Action: remove the short stream line + label + Dis.box + labels;
#            change Effluent line x1 from B to A (FC/MBR output).
# --------------------------------------------------------------------------- #
_DIS_PATTERN = re.compile(
    r'<line\s+x1="([\d.]+)"\s+y1="[\d.]+"\s+x2="[\d.]+"\s+y2="[\d.]+"\s+class="mainflow"[^/]*/>'
    r'<text\b[^>]*class="streamlabel"[^>]*>[^<]+</text>'
    r'<rect\b[^>]*class="unmodeled"[^>]*/>'
    r'<text\b[^>]*class="zonelabel-d">Dis\.</text>'
    r'<text\b[^>]*class="smalllabel">Disinfection</text>'
    r'<text\b[^>]*class="dim">\(not modeled\)</text>'
    r'(<line\s+x1="[\d.]+"\s+y1="[\d.]+"\s+x2="[\d.]+"\s+y2="[\d.]+"\s+class="mainflow"[^/]*/>)'
)

def _remove_dis(m):
    x1_fc    = m.group(1)  # FC/MBR output x-coordinate (from the short SE line)
    effl_old = m.group(2)  # full Effluent line element
    effl_new = re.sub(r'x1="[\d.]+"', f'x1="{x1_fc}"', effl_old, count=1)
    return effl_new

html = _DIS_PATTERN.sub(_remove_dis, html)

# --------------------------------------------------------------------------- #
# 3. Remove "Not modeled" legend entry (Dis. unit removed; entry is now unused)
# --------------------------------------------------------------------------- #
html = re.sub(
    r'<rect\b[^>]*class="unmodeled"[^/]*/><text\b[^>]*class="streamlabel">Not modeled</text>',
    '', html
)

# --------------------------------------------------------------------------- #
# 3b. Rename FC -> SC throughout (except the explanatory MBR note "no FC needed")
# --------------------------------------------------------------------------- #
html = html.replace('>Final Clarifier<', '>Secondary Clarifier<')
html = html.replace('class="zonelabel-d">FC</text>', 'class="zonelabel-d">SC</text>')
html = html.replace('<code>FC</code> Final clarifier', '<code>SC</code> Secondary clarifier')
html = html.replace(
    ' Disinfection is drawn as a dashed grey box\n  downstream of the final clarifier because it is not modeled in the EXPOsan\n  reference implementation.',
    ''
)
html = html.replace('<code>GT</code> Gravity thickener (PS)', '<code>GT</code> Gravity thickener')
html = html.replace('<code>MT</code> Mechanical thickener (WAS)', '<code>MT</code> Mechanical thickener')
html = html.replace(
    '<code>WAS</code>/<code>PS</code> Waste/primary sludge',
    '<code>PS</code> Primary sludge &nbsp; <code>WAS</code> Waste activated sludge'
)

# --------------------------------------------------------------------------- #
# 3c. Spell out abbreviations in stream labels
# --------------------------------------------------------------------------- #
html = html.replace('>WAS</text>', '>WAS (waste activated sludge)</text>')
html = html.replace('>Raw WW</text>', '>Raw WW (wastewater)</text>')

# --------------------------------------------------------------------------- #
# 4. Fix font sizes directly in SVG <style> blocks
# --------------------------------------------------------------------------- #
FONT_REPLACEMENTS = [
    ("font: 700 26px 'Inter'",   "font: 700 32px 'Inter'"),
    ("font: 400 13px 'Inter'",   "font: 400 16px 'Inter'"),
    ("font: 500 9.5px 'Inter'",  "font: 500 12px 'Inter'"),
    ("font: 700 10.5px 'Inter'", "font: 700 13px 'Inter'"),
    ("font: 700 11px 'Inter'",   "font: 700 13px 'Inter'"),
    ("font: 700 9.5px 'Inter'",  "font: 700 11px 'Inter'"),
    ("font: 400 10px 'Inter'",   "font: 400 12px 'Inter'"),
    ("font: 400 9px 'Inter'",    "font: 400 11px 'Inter'"),
]
for old, new in FONT_REPLACEMENTS:
    html = html.replace(old, new)

# --------------------------------------------------------------------------- #
# 4. CSS block (tooltip styles)
# --------------------------------------------------------------------------- #
TOOLTIP_CSS = """
<style id="wrrf-interactive">
  #wrrf-tt {
    position: fixed; pointer-events: none; z-index: 9999;
    background: #1a1a1a; color: #f5f3ee;
    padding: 10px 14px; border-radius: 6px;
    font: 13px/1.5 'Inter','Helvetica Neue',Arial,sans-serif;
    box-shadow: 0 4px 16px rgba(0,0,0,.35);
    opacity: 0; transition: opacity .1s ease;
    max-width: 260px; min-width: 160px;
  }
  #wrrf-tt.show { opacity: 1; }
  #wrrf-tt .tt-name {
    font-weight: 700; font-size: 12px; letter-spacing: .04em;
    text-transform: uppercase; margin-bottom: 6px; padding-bottom: 5px;
    border-bottom: 1px solid #333;
  }
  #wrrf-tt.tt-stream .tt-name { color: #8ac; }
  #wrrf-tt.tt-unit   .tt-name { color: #c8a870; }
  #wrrf-tt .tt-row   { display: flex; justify-content: space-between; gap: 14px; line-height: 1.7; }
  #wrrf-tt .tt-lbl   { color: #888; font-size: 12px; }
  #wrrf-tt .tt-val   { font-weight: 600; font-size: 12px; color: #f5f3ee; text-align: right; }
  #wrrf-tt .tt-desc  { font-size: 12px; color: #ccc; margin-top: 5px; line-height: 1.45; }
  .diagram { position: relative; }
  .back-to-top {
    display: inline-block; margin: 6px 0 2px;
    padding: 4px 14px; font-size: 12px; font-family: 'Inter','Helvetica Neue',Arial,sans-serif;
    color: #555; text-decoration: none;
    border: 1px solid #d9d6cf; border-radius: 12px; background: white;
  }
  .back-to-top:hover { background: #1a1a1a; color: #fbfbf8; border-color: #1a1a1a; }
</style>
"""

# --------------------------------------------------------------------------- #
# 5. Embedded flow JSON (compact)
# --------------------------------------------------------------------------- #
FLOW_JSON = json.dumps(flow_data, separators=(',', ':'))

# --------------------------------------------------------------------------- #
# 6. Interaction JS
# --------------------------------------------------------------------------- #
INTERACT_JS = (
"""
<div id="wrrf-tt"><div class="tt-name" id="tt-nm"></div><div id="tt-bd"></div></div>
<script>
// ── simulation flow data ──────────────────────────────────────────────────────
const FD = """
+ FLOW_JSON +
""";

// ── SVG stream label text -> QSDsan stream ID ─────────────────────────────────
// 'Effluent' maps to 'SE' because the Dis. unit is not modeled — the SE stream
// from FC/MBR IS the final plant effluent in the simulation.
const LABEL_MAP = {
  'Raw WW'          : 'RWW',
  'PE'              : 'PE',
  'ML'              : 'treated',
  'SE'              : 'SE',
  'Effluent'        : 'SE',
  'RAS (return activated sludge)': 'RAS',
  'RAS'             : 'RAS',
  'MLR'             : 'MLR',
  'MLR (mixed liquor recycle)'   : 'MLR',
  'PS'              : 'PS',
  'WAS'             : 'WAS',
  'sludge cake'     : 'cake',
  'biogas'          : 'biogas',
  'reject water (centrate + thickener supernatant)': 'reject',
  'reject water'    : 'reject',
  'permeate'        : 'SE',
  'Permeate'        : 'SE',
  'FeCl3'           : 'FeCl3',
  'alum'            : 'alum',
  'settled effluent': 'settled_effluent',
};

// ── unit descriptions ─────────────────────────────────────────────────────────
const UNIT_DESC = {
  PC:  ['Primary Clarifier',          'Gravity settling removes ~60% of influent TSS. Produces primary effluent (PE) to bioreactor and primary sludge (PS) to solids train. Reject water from downstream returns here.'],
  SC:  ['Secondary Clarifier',         'Gravity settling of mixed liquor. Produces return activated sludge (RAS) recycled to bioreactor head, and waste activated sludge (WAS) to thickening/digestion.'],
  MBR: ['Membrane Bioreactor',        'Submerged membrane replaces the final clarifier for solid-liquid separation. Produces high-quality secondary effluent (SE); WAS wasted periodically.'],
  GT:  ['Gravity Thickener (PS)',     'Passive gravity thickening of primary sludge. Supernatant (reject water) returned to plant headworks. Underflow proceeds to digestion or dewatering.'],
  MT:  ['Mechanical Thickener (WAS)', 'Belt-press or centrifuge thickening of waste activated sludge. Filtrate returned as reject water. Underflow proceeds to digestion or dewatering.'],
  AD:  ['Anaerobic Digester',         'Mesophilic anaerobic digestion at 35 deg C. Stabilises combined thickened sludge, producing biogas and digested sludge for dewatering. Modeled with ADM1.'],
  AED: ['Aerobic Digester',           'Aerobic stabilisation of waste sludge at ambient temperature. No biogas production; reduces VSS before dewatering.'],
  DW:  ['Dewatering',                 'Centrifuge or belt-press dewatering. Produces sludge cake (~20-25% DS) for land application or disposal. High-strength centrate returned to plant headworks.'],
};
function unitDesc(id) {
  if (id in UNIT_DESC) return UNIT_DESC[id];
  if (/^O\\d/.test(id))  return ['Aerobic Zone ' + id,   'Aerobic compartment of the plug-flow bioreactor (DO setpoint ~2 mg/L). Nitrification depends on SRT and DO profile.'];
  if (/^AN\\d/.test(id)) return ['Anaerobic Zone ' + id, 'Anaerobic selector zone. Promotes phosphorus-accumulating organisms (PAOs) by creating feast conditions.'];
  if (/^A\\d/.test(id))  return ['Anoxic Zone ' + id,   'Anoxic denitrification zone. Uses nitrate recycled from downstream aerobic zones (MLR) as electron acceptor.'];
  return null;
}

// ── hover highlight colours ───────────────────────────────────────────────────
const STREAM_HOVER = {
  mainflow: { stroke: '#2055b8', sw: '2.2' },
  recycle:  { stroke: '#e03030', sw: '2'   },
  sludge:   { stroke: '#a07040', sw: '2'   },
  gas:      { stroke: '#8855cc', sw: '2'   },
};
const UNIT_SHAPE_RE = /clarifier|aerobic|anoxic|anaerobic|thickener|digester|mbr|chem|dw/;
function unitHoverStyle(cls) {
  if (cls.includes('clarifier'))                              return {stroke:'#5588cc', sw:'2.2', br:'1.06'};
  if (cls.includes('digester-an') || cls.includes('anaerobic')) return {stroke:'#2a1060', sw:'2.2', br:'1.08'};
  if (cls.includes('digester-ae') || cls.includes('aerobic'))   return {stroke:'#1a4a80', sw:'2',   br:'1.08'};
  if (cls.includes('anoxic'))                                 return {stroke:'#1f5b50', sw:'2',   br:'1.08'};
  if (cls.includes('thickener'))                              return {stroke:'#5a4a20', sw:'2',   br:'1.06'};
  if (cls.includes('dw'))                                     return {stroke:'#4a3a10', sw:'2',   br:'1.06'};
  if (cls.includes('mbr'))                                    return {stroke:'#1a4e75', sw:'2.2', br:'1.08'};
  if (cls.includes('chem'))                                   return {stroke:'#a85a30', sw:'2',   br:'1.06'};
  return {stroke:'#5588cc', sw:'2', br:'1.05'};
}

// ── tooltip DOM ───────────────────────────────────────────────────────────────
const tt   = document.getElementById('wrrf-tt');
const ttNm = document.getElementById('tt-nm');
const ttBd = document.getElementById('tt-bd');

function fmt(v, unit) {
  const n = parseFloat(v);
  if (isNaN(n)) return '-';
  return n.toLocaleString('en-US', {maximumFractionDigits: 0}) + ' ' + unit;
}
function showStream(e, label, sd) {
  tt.className = 'show tt-stream';
  ttNm.textContent = label;
  const q   = sd && sd.q_m3d   != null ? fmt(sd.q_m3d, 'm\\u00b3/d') : '-';
  const tss = sd && sd.tss_kgd != null ? fmt(sd.tss_kgd, 'kg/d')
              : (sd && sd.is_gas ? 'N/A (gas)' : '-');
  ttBd.innerHTML =
    '<div class="tt-row"><span class="tt-lbl">Volumetric flow</span><span class="tt-val">' + q   + '</span></div>' +
    '<div class="tt-row"><span class="tt-lbl">TSS load</span><span class="tt-val">'        + tss + '</span></div>';
  move(e);
}
function showUnit(e, name, desc) {
  tt.className = 'show tt-unit';
  ttNm.textContent = name;
  ttBd.innerHTML = '<div class="tt-desc">' + desc + '</div>';
  move(e);
}
function move(e) {
  const tw = tt.offsetWidth, th = tt.offsetHeight;
  const vw = window.innerWidth,  vh = window.innerHeight;
  tt.style.left = (e.clientX + 18 + tw > vw ? e.clientX - tw - 8 : e.clientX + 18) + 'px';
  tt.style.top  = (e.clientY + 18 + th > vh ? e.clientY - th - 8 : e.clientY + 18) + 'px';
}
function hide() { tt.classList.remove('show'); }

// ── legend detection: elements inside <g transform="..."> are legend items ────
function isInLegend(el, svg) {
  var p = el.parentElement;
  while (p && p !== svg) {
    if (p.tagName.toLowerCase() === 'g' && p.hasAttribute('transform')) return true;
    p = p.parentElement;
  }
  return false;
}

// ── geometry helpers ──────────────────────────────────────────────────────────
function midpointOf(el) {
  const tag = el.tagName.toLowerCase();
  if (tag === 'line') {
    return [(+el.getAttribute('x1') + +el.getAttribute('x2')) / 2,
            (+el.getAttribute('y1') + +el.getAttribute('y2')) / 2];
  }
  const nums = (el.getAttribute('d') || '').match(/[-\\d.]+/g);
  return nums && nums.length >= 2 ? [+nums[0], +nums[1]] : [0, 0];
}
function dist2(a, b) { return (a[0]-b[0])**2 + (a[1]-b[1])**2; }

// ── stream hover apply/reset ──────────────────────────────────────────────────
function applyStreamHover(lines, lblEl, on) {
  lines.forEach(function(item) {
    var h = STREAM_HOVER[item.fc] || STREAM_HOVER.mainflow;
    item.el.style.stroke      = on ? h.stroke : '';
    item.el.style.strokeWidth = on ? h.sw     : '';
  });
  if (lblEl) {
    lblEl.style.fill       = on ? '#2055b8' : '';
    lblEl.style.fontWeight = on ? '700'     : '';
  }
}

// ── unit hover apply/reset ────────────────────────────────────────────────────
function applyUnitHover(shapeEl, txtEl, on) {
  if (shapeEl) {
    var h = unitHoverStyle(shapeEl.getAttribute('class') || '');
    shapeEl.style.stroke      = on ? h.stroke                     : '';
    shapeEl.style.strokeWidth = on ? h.sw                         : '';
    shapeEl.style.filter      = on ? 'brightness(' + h.br + ')' : '';
  }
  if (txtEl) { txtEl.style.fontWeight = on ? '900' : ''; }
}

// ── per-diagram setup ─────────────────────────────────────────────────────────
const FLOW_CLASSES = ['mainflow', 'recycle', 'sludge', 'gas'];

function resolveStreamId(labelText, sysData) {
  if (labelText in sysData) return labelText;
  var mapped = LABEL_MAP[labelText];
  if (mapped && mapped in sysData) return mapped;
  var lo = labelText.toLowerCase();
  var keys = Object.keys(sysData);
  for (var i = 0; i < keys.length; i++) {
    if (keys[i].toLowerCase() === lo) return keys[i];
  }
  var entries = Object.entries(LABEL_MAP);
  for (var j = 0; j < entries.length; j++) {
    var lbl = entries[j][0], sid = entries[j][1];
    if ((labelText.startsWith(lbl) || lbl.startsWith(labelText)) && sid in sysData) return sid;
  }
  return null;
}

function setupDiagram(divEl) {
  var sysId   = divEl.id;
  var sysData = FD[sysId] || {};
  var svg     = divEl.querySelector('svg');
  if (!svg) return;
  var nsv = 'http://www.w3.org/2000/svg';

  // -- collect flow elements (exclude legend) ----------------------------------
  var flowEls = [];
  svg.querySelectorAll('line, path').forEach(function(el) {
    if (isInLegend(el, svg)) return;
    var cls = el.getAttribute('class') || '';
    var fc  = FLOW_CLASSES.find(function(c) { return cls.includes(c); });
    if (fc) flowEls.push({ el: el, fc: fc, mid: midpointOf(el) });
  });

  // -- collect stream labels (exclude legend) ----------------------------------
  var labels = [];
  svg.querySelectorAll('text.streamlabel').forEach(function(el) {
    if (isInLegend(el, svg)) return;
    labels.push({ el: el, text: el.textContent.trim(),
                  pos: [+el.getAttribute('x') || 0, +el.getAttribute('y') || 0] });
  });

  // -- match each label to nearest flow line(s) --------------------------------
  var labelGroups = labels.map(function(lbl) {
    if (!flowEls.length) return { lbl: lbl, lines: [] };
    var dists = flowEls.map(function(f) { return { f: f, d: dist2(f.mid, lbl.pos) }; });
    dists.sort(function(a, b) { return a.d - b.d; });
    var threshold = Math.max(dists[0].d * 4, 2500);
    return { lbl: lbl, lines: dists.filter(function(x) { return x.d <= threshold; }).map(function(x) { return x.f; }) };
  });

  // -- stream hit areas + hover highlights -------------------------------------
  labelGroups.forEach(function(group) {
    var lbl     = group.lbl;
    var lines   = group.lines;
    var sid     = resolveStreamId(lbl.text, sysData);
    var sd      = sid ? sysData[sid] : null;
    var ttLabel = lbl.text;

    var onEnter = function(e) { applyStreamHover(lines, lbl.el, true);  showStream(e, ttLabel, sd); };
    var onLeave = function()  { applyStreamHover(lines, lbl.el, false); hide(); };

    lines.forEach(function(item) {
      var hit = item.el.cloneNode(false);
      hit.setAttribute('stroke', 'transparent');
      hit.setAttribute('stroke-width', '18');
      hit.setAttribute('fill', 'none');
      hit.setAttribute('class', '');
      hit.removeAttribute('marker-end');
      hit.removeAttribute('marker-start');
      hit.removeAttribute('style');
      hit.style.cursor = 'crosshair';
      hit.style.pointerEvents = 'stroke';
      hit.addEventListener('mouseenter', onEnter);
      hit.addEventListener('mousemove',  move);
      hit.addEventListener('mouseleave', onLeave);
      svg.appendChild(hit);
    });

    lbl.el.style.cursor = 'crosshair';
    lbl.el.addEventListener('mouseenter', onEnter);
    lbl.el.addEventListener('mousemove',  move);
    lbl.el.addEventListener('mouseleave', onLeave);
  });

  // -- unit hit rects + hover highlights (exclude legend) ----------------------
  // DOM order in SVGs is always: shape -> zonelabel text -> optional air-bubbles
  // so previousElementSibling of the zonelabel IS the unit shape.
  svg.querySelectorAll('text.zonelabel, text.zonelabel-d').forEach(function(txtEl) {
    if (isInLegend(txtEl, svg)) return;
    var uid = txtEl.textContent.trim();
    var ud  = unitDesc(uid);
    if (!ud) return;
    var uName = ud[0], uDesc = ud[1];

    var prev    = txtEl.previousElementSibling;
    var prevCls = prev ? (prev.getAttribute('class') || '') : '';
    var shapeEl = (prev && UNIT_SHAPE_RE.test(prevCls)) ? prev : null;

    try {
      var bbox = txtEl.getBBox();
      var pad  = 30;
      var hitRect = document.createElementNS(nsv, 'rect');
      hitRect.setAttribute('x',      bbox.x - pad);
      hitRect.setAttribute('y',      bbox.y - pad);
      hitRect.setAttribute('width',  bbox.width  + 2 * pad);
      hitRect.setAttribute('height', bbox.height + 2 * pad);
      hitRect.setAttribute('fill',   'transparent');
      hitRect.setAttribute('stroke', 'none');
      hitRect.style.cursor = 'help';
      hitRect.style.pointerEvents = 'fill';
      hitRect.addEventListener('mouseenter', function(e) { applyUnitHover(shapeEl, txtEl, true);  showUnit(e, uName, uDesc); });
      hitRect.addEventListener('mousemove',  move);
      hitRect.addEventListener('mouseleave', function()  { applyUnitHover(shapeEl, txtEl, false); hide(); });
      svg.appendChild(hitRect);
    } catch(err) { /* getBBox may fail for off-screen SVG */ }
  });

  // -- back-to-top button -----------------------------------------------------
  var btn = document.createElement('a');
  btn.href = '#';
  btn.className = 'back-to-top';
  btn.textContent = '↑ Back to top';
  divEl.appendChild(btn);
}

document.querySelectorAll('div.diagram[id]').forEach(setupDiagram);
document.addEventListener('mousemove', function(e) { if (tt.classList.contains('show')) move(e); });
</script>
"""
)

# --------------------------------------------------------------------------- #
# 7. Assemble output
# --------------------------------------------------------------------------- #
html = html.replace('</head>', TOOLTIP_CSS + '</head>', 1)
html = html.replace(
    '<title>WERF benchmark WRRF configurations — process flow diagrams</title>',
    '<title>WERF benchmark WRRF configurations — interactive process flow diagrams</title>',
    1
)
html = html.replace('</body>', INTERACT_JS + '\n</body>', 1)

with open(OUT_HTML, 'w', encoding='utf-8') as f:
    f.write(html)

print('Written -> ' + OUT_HTML)
