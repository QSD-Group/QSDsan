.. _wrrf_interactive:

WERF Benchmark WRRF Configurations
====================================

Process flow diagrams for all 18 water resource recovery facility (WRRF)
configurations as described in In Zhang et al., 2026 [1]_. These configurations were based on the net-zero energy solutions for WRRFs report by the Water Environment Research Foundation (WERF, now a part of the Water Research Foundation, WRF) [2]_.

Source codes for these configurations can be found in the `werf EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/werf>`_.

To create any of the system configurations, use the following code, replacing ``B1`` with the desired configuration code (e.g., ``B2``, ``C3``, etc.). 

.. code-block:: python

   from exposan.werf import create_system

   B1 = create_system('B1')
   B1.simulate(method='BDF', t_span=(0, 400))
   B1.diagram(file='B1', format='html')


Valid configuration codes are: B1, B2, B3, C1, C2, C3, E2, E2P, F1, G1, G2, G3, H1, I1, I2, I3, N1, N2.

.. figure:: ../images/wrrf_configs_light.png
   :class: only-light

.. figure:: ../images/wrrf_configs_dark.png
   :class: only-dark

   Distinguishing features of benchmark WRRF configurations. A WRRF configuration is referred to as a unique combination of a liquid code (a) and a solid code (b) as defined in [2]_. Configurations followed by a star (\*) are implemented in ``EXPOsan``.


.. raw:: html

   <style>
   .wrrf-diagram {
     margin: 18px 0 28px 0; padding: 16px;
     background: white; border: 1px solid #d9d6cf; border-radius: 6px;
     box-shadow: 0 1px 2px rgba(0,0,0,0.03); position: relative;
   }
   .wrrf-diagram svg { width: 100%; height: auto; display: block; }

   /* Diagram SVG classes. Scoped to .wrrf-diagram so generic names
      (.title, .note, .section, ...) do not leak to the page or footer. */
   .wrrf-diagram .bg { fill: #fbfbf8; }
   .wrrf-diagram .title { font: 700 32px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #1a1a1a; }
   .wrrf-diagram .subtitle { font: 400 21px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #555; }
   .wrrf-diagram .smalllabel { font: 500 21px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #444; text-anchor: middle; }
   .wrrf-diagram .streamlabel { font: 500 21px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #444; }
   .wrrf-diagram .zonelabel { font: 700 21px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #fff; text-anchor: middle; }
   .wrrf-diagram .zonelabel-d { font: 700 21px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #1a1a1a; text-anchor: middle; }
   .wrrf-diagram .section { font: 700 18px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #888; text-anchor: middle; letter-spacing:2px; }
   .wrrf-diagram .note { font: 400 21px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #555; }
   .wrrf-diagram .dim { font: 400 18px 'Inter','Helvetica Neue',Arial,sans-serif; fill: #888; text-anchor: middle; }
   .wrrf-diagram .unmodeled { fill: #ececec; stroke: #aaa; stroke-width:1; stroke-dasharray:3,2; }
   .wrrf-diagram .anaerobic { fill: #6b4a8a; stroke: #4a3463; stroke-width:1; }
   .wrrf-diagram .anoxic { fill: #2f7d6e; stroke: #1f5b50; stroke-width:1; }
   .wrrf-diagram .aerobic { fill: #4a90c0; stroke: #2e6b95; stroke-width:1; }
   .wrrf-diagram .clarifier { fill: #e0d4b8; stroke: #8c7a55; stroke-width:1.2; }
   .wrrf-diagram .digester-an { fill: #6b4a8a; stroke: #4a3463; stroke-width:1.2; }
   .wrrf-diagram .digester-ae { fill: #4a90c0; stroke: #2e6b95; stroke-width:1.2; }
   .wrrf-diagram .thickener { fill: #c8bb9c; stroke: #8c7a55; stroke-width:1; }
   .wrrf-diagram .mbr { fill: #2e6b95; stroke: #1a4e75; stroke-width:1; }
   .wrrf-diagram .chem { fill: #d97a4a; stroke: #a85a30; stroke-width:1; }
   .wrrf-diagram .dw { fill: #b8a878; stroke: #7a6a40; stroke-width:1; }
   .wrrf-diagram .mainflow { stroke: #1a1a1a; stroke-width:1.6; fill:none; }
   .wrrf-diagram .recycle { stroke: #c2272d; stroke-width:1.4; fill:none; stroke-dasharray: 5,3; }
   .wrrf-diagram .sludge { stroke: #7a5a30; stroke-width:1.5; fill:none; }
   .wrrf-diagram .gas { stroke: #6b4a8a; stroke-width:1.4; fill:none; stroke-dasharray: 3,2; }
   .wrrf-diagram .dose { stroke: #d97a4a; stroke-width:1.4; fill:none; }

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

   .wrrf-nav {
     margin: 20px 0 28px 0; border: 1px solid #d9d6cf; border-radius: 6px;
     overflow: hidden; font: 13px/1.5 'Inter','Helvetica Neue',Arial,sans-serif;
   }
   .wrrf-nav-row {
     display: flex; align-items: center; gap: 14px; padding: 7px 16px;
     border-bottom: 1px solid #e8e5de; background: #faf9f6;
   }
   .wrrf-nav-row:last-child { border-bottom: none; }
   .wrrf-nav-row:nth-child(even) { background: #f5f3ee; }
   .wrrf-nav-series {
     min-width: 230px; color: #444; font-weight: 500;
     font-size: 13px; text-decoration: none;
   }
   .wrrf-nav-series:hover { color: #2a5a9a; text-decoration: underline; }
   .wrrf-nav-links { display: flex; gap: 6px; flex-wrap: wrap; }
   .wrrf-nav-links a {
     display: inline-block; padding: 2px 10px; background: #fff;
     border: 1px solid #c8c4bb; border-radius: 4px; color: #2a5a9a;
     font-weight: 600; font-size: 13px; text-decoration: none;
     transition: background .12s, border-color .12s;
   }
   .wrrf-nav-links a:hover { background: #e8f0fb; border-color: #2a5a9a; }
   </style>


Reading the Diagrams
--------------------

Every diagram uses the same convention: liquid train across the top with the
bioreactor zones in the middle, mixed-liquor recycles (MLR) above the zones,
return activated sludge (RAS) below them; solids train across the middle;
reject water (centrate + thickener supernatants) returned to the head of the
plant along the bottom.

``PC`` Primary clarifier · ``SC`` Secondary clarifier · ``DW`` Dewatering 

``GT`` Gravity thickener · ``MT`` Mechanical thickener · ``MBR`` Membrane bioreactor

``AD`` Anaerobic digester · ``AED`` Aerobic digester 

``RWW`` Raw wastewater · ``ML`` Mixed liquor

``PE`` Primary effluent · ``PS`` Primary sludge 

``SE`` Secondary effluent · ``WAS`` Waste activated sludge

The configurations are categorized based on their treatment goals:
  * rBOD for basic BOD removal (r for readily): B (with primary clarifier) and C (without primary clarifier) series.
  * NIT for nitrification: E (aerobic digestion) and F series (anaerobic digestion).
  * BNR for biological nutrient removal: G (Johannesburg), H (4-stage Bardenpho+chemicals), I (5-stage Bardenpho), N (5-stage Bardenpho+MBR) series.


.. raw:: html

   <div class="wrrf-nav">
     <div class="wrrf-nav-row">
       <a class="wrrf-nav-series" href="#b-series-rbod-with-primary-clarifier">B &nbsp;·&nbsp; rBOD with primary clarifier</a>
       <span class="wrrf-nav-links">
         <a href="#B1" title="rBOD · Primary clarifier · 6-zone aerobic ASP · Anaerobic digestion">B1</a>
         <a href="#B2" title="rBOD · Primary clarifier · 6-zone aerobic ASP · Aerobic digestion">B2</a>
         <a href="#B3" title="rBOD · Primary clarifier · 6-zone aerobic ASP · No digestion">B3</a>
       </span>
     </div>
     <div class="wrrf-nav-row">
       <a class="wrrf-nav-series" href="#c-series-rbod-without-primary-clarifier">C &nbsp;·&nbsp; rBOD, no primary clarifier</a>
       <span class="wrrf-nav-links">
         <a href="#C1" title="rBOD · No primary clarifier · 6-zone aerobic ASP · Anaerobic digestion">C1</a>
         <a href="#C2" title="rBOD · No primary clarifier · 6-zone aerobic ASP · Aerobic digestion">C2</a>
         <a href="#C3" title="rBOD · No primary clarifier · 6-zone aerobic ASP · No digestion">C3</a>
       </span>
     </div>
     <div class="wrrf-nav-row">
       <a class="wrrf-nav-series" href="#e-series-nit-with-aerobic-digestion">E &nbsp;·&nbsp; NIT, aerobic digestion</a>
       <span class="wrrf-nav-links">
         <a href="#E2"  title="NIT · No primary clarifier · Nitrifying ASP · Aerobic digestion">E2</a>
         <a href="#E2P" title="NIT · Primary clarifier (prime variant) · Nitrifying ASP · Aerobic digestion">E2P</a>
       </span>
     </div>
     <div class="wrrf-nav-row">
       <a class="wrrf-nav-series" href="#f-series-nit-with-anaerobic-digestion">F &nbsp;·&nbsp; NIT, anaerobic digestion</a>
       <span class="wrrf-nav-links">
         <a href="#F1" title="NIT · Primary clarifier · Nitrifying ASP · Anaerobic digestion">F1</a>
       </span>
     </div>
     <div class="wrrf-nav-row">
       <a class="wrrf-nav-series" href="#g-series-bnr-with-johannesburg-process">G &nbsp;·&nbsp; BNR, Johannesburg</a>
       <span class="wrrf-nav-links">
         <a href="#G1" title="BNR · Primary clarifier · Johannesburg step-fed · Anaerobic digestion">G1</a>
         <a href="#G2" title="BNR · Primary clarifier · Johannesburg step-fed · Aerobic digestion">G2</a>
         <a href="#G3" title="BNR · Primary clarifier · Johannesburg step-fed · No digestion">G3</a>
       </span>
     </div>
     <div class="wrrf-nav-row">
       <a class="wrrf-nav-series" href="#h-series-bnr-with-chemical-p-removal">H &nbsp;·&nbsp; BNR, chemical P removal</a>
       <span class="wrrf-nav-links">
         <a href="#H1" title="BNR · Primary clarifier + chemical P removal · Modified 5-stage Bardenpho · Anaerobic digestion">H1</a>
       </span>
     </div>
     <div class="wrrf-nav-row">
       <a class="wrrf-nav-series" href="#i-series-bnr-with-5-stage-bardenpho">I &nbsp;·&nbsp; BNR, 5-stage Bardenpho</a>
       <span class="wrrf-nav-links">
         <a href="#I1" title="BNR · No primary clarifier · 5-stage Bardenpho · Anaerobic digestion">I1</a>
         <a href="#I2" title="BNR · No primary clarifier · 5-stage Bardenpho · Aerobic digestion">I2</a>
         <a href="#I3" title="BNR · No primary clarifier · 5-stage Bardenpho · No digestion">I3</a>
       </span>
     </div>
     <div class="wrrf-nav-row">
       <a class="wrrf-nav-series" href="#n-series-bnr-with-mbr">N &nbsp;·&nbsp; BNR + MBR</a>
       <span class="wrrf-nav-links">
         <a href="#N1" title="BNR · Primary clarifier · 5-stage Bardenpho + MBR · Anaerobic digestion">N1</a>
         <a href="#N2" title="BNR · No primary clarifier · 5-stage Bardenpho + MBR · Aerobic digestion">N2</a>
       </span>
     </div>
   </div>


B Series: rBOD with Primary Clarifier
-------------------------------------

Conventional six-zone aerobic activated sludge with short SRT (~2 d) — only BOD
removal, nitrifiers wash out. The three B variants differ only in the solids
train (AD, AED, or no digestion).

.. raw:: html

   <div class="wrrf-diagram" id="B1">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">B1</text><text x="36" y="78" class="subtitle">rBOD · Primary clarifier · 6-zone aerobic ASP · Anaerobic digestion</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="118" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><polygon points="118,205 182,205 175,261 125,261" class="clarifier"/><text x="150.0" y="238.0" class="zonelabel-d">PC</text><line x1="182" y1="233.0" x2="210" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="196.0" y="225.0" class="streamlabel" text-anchor="middle">PE</text><rect x="210" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="238.0" y="237.0" class="zonelabel">O1</text><circle cx="230.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="238.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="246.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="266" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="294.0" y="237.0" class="zonelabel">O2</text><circle cx="286.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="294.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="302.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="322" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="350.0" y="237.0" class="zonelabel">O3</text><circle cx="342.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="350.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="358.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="378" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="406.0" y="237.0" class="zonelabel">O4</text><circle cx="398.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="406.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="414.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="434" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="462.0" y="237.0" class="zonelabel">O5</text><circle cx="454.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="462.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="470.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="490" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="518.0" y="237.0" class="zonelabel">O6</text><circle cx="510.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="518.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="526.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><line x1="546" y1="233.0" x2="586" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="566.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="586,205 650,205 643,261 593,261" class="clarifier"/><text x="618.0" y="238.0" class="zonelabel-d">SC</text><line x1="650" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="955.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 600 261 L 600 320 L 218 320 L 218 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="409.0" y="315" class="streamlabel" text-anchor="middle">RAS (return activated sludge)</text><line x1="150.0" y1="261" x2="150.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="142.0" y="275" class="streamlabel" text-anchor="end">PS</text><line x1="636" y1="261" x2="636" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="644" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="319.0,470 375.0,470 369.0,526 325.0,526" class="thickener"/><text x="347.0" y="503.0" class="zonelabel-d">GT</text><polygon points="485.0,470 541.0,470 535.0,526 491.0,526" class="thickener"/><text x="513.0" y="503.0" class="zonelabel-d">MT</text><path d="M651.0,480 Q651.0,470 683.0,470 Q715.0,470 715.0,480 L715.0,526 L651.0,526 z" class="digester-an"/><text x="683.0" y="503.0" class="zonelabel">AD</text><line x1="683.0" y1="468" x2="683.0" y2="448" class="gas" marker-end="url(#arr-gas)"/><text x="688.0" y="445" class="streamlabel" fill="#6b4a8a">biogas</text><rect x="825.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="853.0" y="503.0" class="zonelabel-d">DW</text><line x1="375.0" y1="498.0" x2="485.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="541.0" y1="498.0" x2="651.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="715.0" y1="498.0" x2="825.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="150.0" y1="420" x2="347.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="347.0" y1="420" x2="347.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="636" y1="420" x2="513.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="513.0" y1="420" x2="513.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="881.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="1016.5" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="347.0" y1="526" x2="347.0" y2="600" class="recycle"/><line x1="513.0" y1="526" x2="513.0" y2="600" class="recycle"/><line x1="853.0" y1="526" x2="853.0" y2="600" class="recycle"/><line x1="347.0" y1="600" x2="853.0" y2="600" class="recycle"/><line x1="347.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="116" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="112" y1="255.0" x2="112" y2="236.0" class="recycle"/><line x1="112" y1="236.0" x2="116" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
      <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
      <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
      <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
      <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
      <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
      <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
      <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
      <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
      <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
      <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
      <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
      <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
      <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
      <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
      <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
      </g></svg>
   </div>

   <div id="wrrf-tt"><div class="tt-name" id="tt-nm"></div><div id="tt-bd"></div></div>
   <script>
   const FD = {"B1":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":992.9,"tss_kgd":765.97,"is_gas":false,"phase":"l"},"PE":{"id":"PE","q_m3d":38597.18,"tss_kgd":2792.68,"is_gas":false,"phase":"l"},"PS":{"id":"PS","q_m3d":249.84,"tss_kgd":4234.26,"is_gas":false,"phase":"l"},"ws3":{"id":"ws3","q_m3d":162.77,"tss_kgd":275.86,"is_gas":false,"phase":"l"},"thickened_PS":{"id":"thickened_PS","q_m3d":87.06,"tss_kgd":3958.15,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":25362.26,"tss_kgd":122804.85,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":63959.44,"tss_kgd":127217.48,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37840.1,"tss_kgd":746.82,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":757.08,"tss_kgd":3665.8,"is_gas":false,"phase":"l"},"ws6":{"id":"ws6","q_m3d":679.1,"tss_kgd":164.41,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":77.98,"tss_kgd":3501.43,"is_gas":false,"phase":"l"},"ws8":{"id":"ws8","q_m3d":165.04,"tss_kgd":7459.57,"is_gas":false,"phase":"l"},"ws9":{"id":"ws9","q_m3d":165.04,"tss_kgd":6219.38,"is_gas":false,"phase":"l"},"biogas":{"id":"biogas","q_m3d":5.99,"tss_kgd":null,"is_gas":true,"phase":"g"},"digestate":{"id":"digestate","q_m3d":165.04,"tss_kgd":3675.76,"is_gas":false,"phase":"l"},"ws13":{"id":"ws13","q_m3d":165.04,"tss_kgd":3715.38,"is_gas":false,"phase":"l"},"ws14":{"id":"ws14","q_m3d":144.68,"tss_kgd":325.7,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":20.37,"tss_kgd":3390.51,"is_gas":false,"phase":"l"},"ws16":{"id":"ws16","q_m3d":986.55,"tss_kgd":765.97,"is_gas":false,"phase":"l"}},"B2":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":999.49,"tss_kgd":559.2,"is_gas":false,"phase":"l"},"PE":{"id":"PE","q_m3d":38603.77,"tss_kgd":2710.51,"is_gas":false,"phase":"l"},"PS":{"id":"PS","q_m3d":249.84,"tss_kgd":4109.66,"is_gas":false,"phase":"l"},"ws19":{"id":"ws19","q_m3d":162.77,"tss_kgd":267.74,"is_gas":false,"phase":"l"},"thickened_PS":{"id":"thickened_PS","q_m3d":87.06,"tss_kgd":3841.67,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":25362.26,"tss_kgd":118247.21,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":63966.03,"tss_kgd":122519.67,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37846.69,"tss_kgd":742.71,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":757.08,"tss_kgd":3529.76,"is_gas":false,"phase":"l"},"ws22":{"id":"ws22","q_m3d":681.37,"tss_kgd":158.84,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":75.71,"tss_kgd":3371.01,"is_gas":false,"phase":"l"},"digestate":{"id":"digestate","q_m3d":162.77,"tss_kgd":1442.19,"is_gas":false,"phase":"l"},"ws25":{"id":"ws25","q_m3d":149.68,"tss_kgd":132.62,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":13.1,"tss_kgd":1309.84,"is_gas":false,"phase":"l"},"ws27":{"id":"ws27","q_m3d":993.82,"tss_kgd":559.2,"is_gas":false,"phase":"l"}},"B3":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":972.02,"tss_kgd":1022.08,"is_gas":false,"phase":"l"},"PE":{"id":"PE","q_m3d":38576.3,"tss_kgd":2894.46,"is_gas":false,"phase":"l"},"PS":{"id":"PS","q_m3d":249.84,"tss_kgd":4388.6,"is_gas":false,"phase":"l"},"ws30":{"id":"ws30","q_m3d":162.77,"tss_kgd":285.92,"is_gas":false,"phase":"l"},"thickened_PS":{"id":"thickened_PS","q_m3d":87.06,"tss_kgd":4102.42,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":25362.26,"tss_kgd":123234.24,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":63938.56,"tss_kgd":127659.35,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37819.22,"tss_kgd":746.47,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":757.08,"tss_kgd":3678.62,"is_gas":false,"phase":"l"},"ws33":{"id":"ws33","q_m3d":682.89,"tss_kgd":165.91,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":74.19,"tss_kgd":3512.53,"is_gas":false,"phase":"l"},"ws35":{"id":"ws35","q_m3d":161.26,"tss_kgd":7615.42,"is_gas":false,"phase":"l"},"ws36":{"id":"ws36","q_m3d":120.75,"tss_kgd":570.24,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":40.5,"tss_kgd":7044.42,"is_gas":false,"phase":"l"},"ws38":{"id":"ws38","q_m3d":966.42,"tss_kgd":1022.09,"is_gas":false,"phase":"l"}},"C1":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":25362.26,"tss_kgd":149319.94,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":1188.66,"tss_kgd":627.83,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":64405.04,"tss_kgd":156760.05,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37831.45,"tss_kgd":308.41,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":1211.33,"tss_kgd":7131.69,"is_gas":false,"phase":"l"},"ws40":{"id":"ws40","q_m3d":1053.86,"tss_kgd":310.23,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":157.47,"tss_kgd":6821.33,"is_gas":false,"phase":"l"},"ws42":{"id":"ws42","q_m3d":157.47,"tss_kgd":5557.11,"is_gas":false,"phase":"l"},"biogas":{"id":"biogas","q_m3d":4.61,"tss_kgd":null,"is_gas":true,"phase":"g"},"digestate":{"id":"digestate","q_m3d":157.47,"tss_kgd":3711.87,"is_gas":false,"phase":"l"},"ws46":{"id":"ws46","q_m3d":157.47,"tss_kgd":3732.17,"is_gas":false,"phase":"l"},"ws47":{"id":"ws47","q_m3d":134.0,"tss_kgd":317.59,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":23.47,"tss_kgd":3414.71,"is_gas":false,"phase":"l"},"ws49":{"id":"ws49","q_m3d":1187.86,"tss_kgd":627.83,"is_gas":false,"phase":"l"}},"C2":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":25362.26,"tss_kgd":144142.1,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":1197.58,"tss_kgd":431.91,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":64413.96,"tss_kgd":151332.34,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37840.36,"tss_kgd":305.86,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":1211.33,"tss_kgd":6884.39,"is_gas":false,"phase":"l"},"ws51":{"id":"ws51","q_m3d":1056.13,"tss_kgd":300.12,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":155.2,"tss_kgd":6584.2,"is_gas":false,"phase":"l"},"digestate":{"id":"digestate","q_m3d":155.2,"tss_kgd":1447.45,"is_gas":false,"phase":"l"},"ws54":{"id":"ws54","q_m3d":141.31,"tss_kgd":131.79,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":13.89,"tss_kgd":1315.45,"is_gas":false,"phase":"l"},"ws56":{"id":"ws56","q_m3d":1197.44,"tss_kgd":431.91,"is_gas":false,"phase":"l"}},"C3":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":25362.26,"tss_kgd":140632.76,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":1169.83,"tss_kgd":759.52,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":64386.21,"tss_kgd":147653.2,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37812.62,"tss_kgd":303.65,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":1211.33,"tss_kgd":6716.78,"is_gas":false,"phase":"l"},"ws58":{"id":"ws58","q_m3d":1059.92,"tss_kgd":293.86,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":151.42,"tss_kgd":6423.08,"is_gas":false,"phase":"l"},"ws60":{"id":"ws60","q_m3d":109.78,"tss_kgd":465.68,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":41.64,"tss_kgd":5957.33,"is_gas":false,"phase":"l"},"ws62":{"id":"ws62","q_m3d":1169.69,"tss_kgd":759.52,"is_gas":false,"phase":"l"}},"E2":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":25362.26,"tss_kgd":114910.77,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":1236.84,"tss_kgd":361.44,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":64453.22,"tss_kgd":120860.82,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37841.77,"tss_kgd":290.29,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":1249.19,"tss_kgd":5659.8,"is_gas":false,"phase":"l"},"ws64":{"id":"ws64","q_m3d":1116.7,"tss_kgd":252.98,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":132.49,"tss_kgd":5406.83,"is_gas":false,"phase":"l"},"digestate":{"id":"digestate","q_m3d":132.49,"tss_kgd":1197.54,"is_gas":false,"phase":"l"},"ws67":{"id":"ws67","q_m3d":120.0,"tss_kgd":108.47,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":12.49,"tss_kgd":1088.91,"is_gas":false,"phase":"l"},"ws69":{"id":"ws69","q_m3d":1236.69,"tss_kgd":361.44,"is_gas":false,"phase":"l"}},"E2P":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":876.78,"tss_kgd":748.06,"is_gas":false,"phase":"l"},"PE":{"id":"PE","q_m3d":38450.78,"tss_kgd":2783.32,"is_gas":false,"phase":"l"},"PS":{"id":"PS","q_m3d":280.12,"tss_kgd":4225.66,"is_gas":false,"phase":"l"},"ws72":{"id":"ws72","q_m3d":193.06,"tss_kgd":291.23,"is_gas":false,"phase":"l"},"thickened_PS":{"id":"thickened_PS","q_m3d":87.06,"tss_kgd":3934.23,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":25362.26,"tss_kgd":132458.34,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":63813.04,"tss_kgd":135979.33,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37833.76,"tss_kgd":298.49,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":617.02,"tss_kgd":3222.48,"is_gas":false,"phase":"l"},"ws75":{"id":"ws75","q_m3d":545.1,"tss_kgd":142.34,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":71.92,"tss_kgd":3080.03,"is_gas":false,"phase":"l"},"digestate":{"id":"digestate","q_m3d":158.99,"tss_kgd":3763.18,"is_gas":false,"phase":"l"},"ws78":{"id":"ws78","q_m3d":132.87,"tss_kgd":314.49,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":26.12,"tss_kgd":3448.72,"is_gas":false,"phase":"l"},"ws80":{"id":"ws80","q_m3d":871.02,"tss_kgd":748.06,"is_gas":false,"phase":"l"}},"F1":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":881.66,"tss_kgd":756.02,"is_gas":false,"phase":"l"},"PE":{"id":"PE","q_m3d":38455.66,"tss_kgd":2786.48,"is_gas":false,"phase":"l"},"PS":{"id":"PS","q_m3d":280.12,"tss_kgd":4230.45,"is_gas":false,"phase":"l"},"ws83":{"id":"ws83","q_m3d":181.7,"tss_kgd":274.41,"is_gas":false,"phase":"l"},"thickened_PS":{"id":"thickened_PS","q_m3d":98.42,"tss_kgd":3956.02,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":25362.26,"tss_kgd":134505.87,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":63817.92,"tss_kgd":138081.12,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37838.63,"tss_kgd":302.93,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":617.02,"tss_kgd":3272.3,"is_gas":false,"phase":"l"},"ws86":{"id":"ws86","q_m3d":542.83,"tss_kgd":143.94,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":74.19,"tss_kgd":3128.19,"is_gas":false,"phase":"l"},"ws88":{"id":"ws88","q_m3d":172.61,"tss_kgd":7084.22,"is_gas":false,"phase":"l"},"ws89":{"id":"ws89","q_m3d":172.61,"tss_kgd":6037.28,"is_gas":false,"phase":"l"},"biogas":{"id":"biogas","q_m3d":5.14,"tss_kgd":null,"is_gas":true,"phase":"g"},"digestate":{"id":"digestate","q_m3d":172.61,"tss_kgd":3825.99,"is_gas":false,"phase":"l"},"ws93":{"id":"ws93","q_m3d":172.61,"tss_kgd":3865.77,"is_gas":false,"phase":"l"},"ws94":{"id":"ws94","q_m3d":150.77,"tss_kgd":337.66,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":21.84,"tss_kgd":3527.91,"is_gas":false,"phase":"l"},"ws96":{"id":"ws96","q_m3d":875.3,"tss_kgd":756.02,"is_gas":false,"phase":"l"}},"G1":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":770.5,"tss_kgd":966.32,"is_gas":false,"phase":"l"},"PE":{"id":"PE","q_m3d":38344.5,"tss_kgd":2869.93,"is_gas":false,"phase":"l"},"PS":{"id":"PS","q_m3d":280.12,"tss_kgd":4357.3,"is_gas":false,"phase":"l"},"ws99":{"id":"ws99","q_m3d":181.7,"tss_kgd":282.64,"is_gas":false,"phase":"l"},"thickened_PS":{"id":"thickened_PS","q_m3d":98.42,"tss_kgd":4074.64,"is_gas":false,"phase":"l"},"ws101":{"id":"ws101","q_m3d":30675.6,"tss_kgd":2295.94,"is_gas":false,"phase":"l"},"ws102":{"id":"ws102","q_m3d":7668.9,"tss_kgd":573.99,"is_gas":false,"phase":"l"},"carbon":{"id":"carbon","q_m3d":2.17,"tss_kgd":0.0,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":15141.65,"tss_kgd":158073.47,"is_gas":false,"phase":"l"},"ws103":{"id":"ws103","q_m3d":15143.82,"tss_kgd":157864.23,"is_gas":false,"phase":"l"},"ws105":{"id":"ws105","q_m3d":45819.42,"tss_kgd":156711.47,"is_gas":false,"phase":"l"},"intr":{"id":"intr","q_m3d":152823.77,"tss_kgd":467746.27,"is_gas":false,"phase":"l"},"ws107":{"id":"ws107","q_m3d":206312.09,"tss_kgd":628222.4,"is_gas":false,"phase":"l"},"ws109":{"id":"ws109","q_m3d":206312.09,"tss_kgd":629990.47,"is_gas":false,"phase":"l"},"ws111":{"id":"ws111","q_m3d":206312.09,"tss_kgd":632067.49,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":53488.32,"tss_kgd":163711.2,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37835.64,"tss_kgd":302.97,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":511.03,"tss_kgd":5334.97,"is_gas":false,"phase":"l"},"ws115":{"id":"ws115","q_m3d":419.8,"tss_kgd":219.13,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":91.23,"tss_kgd":5115.94,"is_gas":false,"phase":"l"},"ws117":{"id":"ws117","q_m3d":189.65,"tss_kgd":9190.56,"is_gas":false,"phase":"l"},"ws118":{"id":"ws118","q_m3d":189.65,"tss_kgd":8348.93,"is_gas":false,"phase":"l"},"biogas":{"id":"biogas","q_m3d":5.3,"tss_kgd":null,"is_gas":true,"phase":"g"},"digestate":{"id":"digestate","q_m3d":189.65,"tss_kgd":5537.81,"is_gas":false,"phase":"l"},"ws122":{"id":"ws122","q_m3d":189.65,"tss_kgd":5418.9,"is_gas":false,"phase":"l"},"ws123":{"id":"ws123","q_m3d":162.58,"tss_kgd":464.54,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":27.07,"tss_kgd":4955.11,"is_gas":false,"phase":"l"},"ws126":{"id":"ws126","q_m3d":764.09,"tss_kgd":966.32,"is_gas":false,"phase":"l"}},"G2":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":785.2,"tss_kgd":767.25,"is_gas":false,"phase":"l"},"PE":{"id":"PE","q_m3d":38359.2,"tss_kgd":2790.89,"is_gas":false,"phase":"l"},"PS":{"id":"PS","q_m3d":280.12,"tss_kgd":4237.28,"is_gas":false,"phase":"l"},"ws129":{"id":"ws129","q_m3d":181.7,"tss_kgd":274.85,"is_gas":false,"phase":"l"},"thickened_PS":{"id":"thickened_PS","q_m3d":98.42,"tss_kgd":3962.4,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":15141.65,"tss_kgd":154940.24,"is_gas":false,"phase":"l"},"carbon":{"id":"carbon","q_m3d":2.17,"tss_kgd":0.0,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":53503.02,"tss_kgd":160510.88,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37846.56,"tss_kgd":303.05,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":514.82,"tss_kgd":5268.01,"is_gas":false,"phase":"l"},"ws132":{"id":"ws132","q_m3d":425.1,"tss_kgd":217.5,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":89.71,"tss_kgd":5050.23,"is_gas":false,"phase":"l"},"digestate":{"id":"digestate","q_m3d":188.13,"tss_kgd":2993.51,"is_gas":false,"phase":"l"},"ws135":{"id":"ws135","q_m3d":172.77,"tss_kgd":274.91,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":15.37,"tss_kgd":2718.9,"is_gas":false,"phase":"l"},"ws137":{"id":"ws137","q_m3d":779.57,"tss_kgd":767.25,"is_gas":false,"phase":"l"}},"G3":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":776.68,"tss_kgd":1075.97,"is_gas":false,"phase":"l"},"PE":{"id":"PE","q_m3d":38350.67,"tss_kgd":2913.47,"is_gas":false,"phase":"l"},"PS":{"id":"PS","q_m3d":280.12,"tss_kgd":4423.4,"is_gas":false,"phase":"l"},"ws140":{"id":"ws140","q_m3d":181.7,"tss_kgd":286.92,"is_gas":false,"phase":"l"},"thickened_PS":{"id":"thickened_PS","q_m3d":98.42,"tss_kgd":4136.46,"is_gas":false,"phase":"l"},"ws142":{"id":"ws142","q_m3d":30680.54,"tss_kgd":2330.78,"is_gas":false,"phase":"l"},"ws143":{"id":"ws143","q_m3d":7670.13,"tss_kgd":582.69,"is_gas":false,"phase":"l"},"carbon":{"id":"carbon","q_m3d":0.77,"tss_kgd":0.0,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":15141.65,"tss_kgd":115489.65,"is_gas":false,"phase":"l"},"ws144":{"id":"ws144","q_m3d":15142.41,"tss_kgd":115691.94,"is_gas":false,"phase":"l"},"ws146":{"id":"ws146","q_m3d":45822.95,"tss_kgd":115076.36,"is_gas":false,"phase":"l"},"intr":{"id":"intr","q_m3d":152837.39,"tss_kgd":342384.41,"is_gas":false,"phase":"l"},"ws148":{"id":"ws148","q_m3d":206330.48,"tss_kgd":460708.77,"is_gas":false,"phase":"l"},"ws150":{"id":"ws150","q_m3d":206330.48,"tss_kgd":462096.51,"is_gas":false,"phase":"l"},"ws152":{"id":"ws152","q_m3d":206330.48,"tss_kgd":462798.01,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":53493.09,"tss_kgd":119834.55,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37821.48,"tss_kgd":302.76,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":529.96,"tss_kgd":4042.15,"is_gas":false,"phase":"l"},"ws156":{"id":"ws156","q_m3d":461.82,"tss_kgd":176.12,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":68.14,"tss_kgd":3866.16,"is_gas":false,"phase":"l"},"ws158":{"id":"ws158","q_m3d":166.56,"tss_kgd":8002.59,"is_gas":false,"phase":"l"},"ws159":{"id":"ws159","q_m3d":127.57,"tss_kgd":612.93,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":38.99,"tss_kgd":7389.63,"is_gas":false,"phase":"l"},"ws161":{"id":"ws161","q_m3d":771.09,"tss_kgd":1075.97,"is_gas":false,"phase":"l"}},"H1":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"ws162":{"id":"ws162","q_m3d":37854.12,"tss_kgd":7379.59,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":642.54,"tss_kgd":904.32,"is_gas":false,"phase":"l"},"PE":{"id":"PE","q_m3d":38216.53,"tss_kgd":3289.45,"is_gas":false,"phase":"l"},"PS":{"id":"PS","q_m3d":280.12,"tss_kgd":4994.44,"is_gas":false,"phase":"l"},"ws165":{"id":"ws165","q_m3d":181.7,"tss_kgd":323.96,"is_gas":false,"phase":"l"},"thickened_PS":{"id":"thickened_PS","q_m3d":98.42,"tss_kgd":4670.46,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":15141.65,"tss_kgd":151535.65,"is_gas":false,"phase":"l"},"carbon":{"id":"carbon","q_m3d":1.1,"tss_kgd":0.0,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":53359.28,"tss_kgd":155627.0,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37839.09,"tss_kgd":302.98,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":378.54,"tss_kgd":3788.38,"is_gas":false,"phase":"l"},"ws168":{"id":"ws168","q_m3d":312.41,"tss_kgd":156.33,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":66.13,"tss_kgd":3632.0,"is_gas":false,"phase":"l"},"ws170":{"id":"ws170","q_m3d":164.55,"tss_kgd":8302.46,"is_gas":false,"phase":"l"},"ws171":{"id":"ws171","q_m3d":164.55,"tss_kgd":7272.15,"is_gas":false,"phase":"l"},"biogas":{"id":"biogas","q_m3d":5.44,"tss_kgd":null,"is_gas":true,"phase":"g"},"digestate":{"id":"digestate","q_m3d":164.55,"tss_kgd":4934.94,"is_gas":false,"phase":"l"},"ws175":{"id":"ws175","q_m3d":164.55,"tss_kgd":4952.23,"is_gas":false,"phase":"l"},"ws176":{"id":"ws176","q_m3d":140.89,"tss_kgd":424.02,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":23.66,"tss_kgd":4528.48,"is_gas":false,"phase":"l"},"ws178":{"id":"ws178","q_m3d":635.0,"tss_kgd":904.31,"is_gas":false,"phase":"l"}},"I1":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":25362.26,"tss_kgd":143955.08,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":1046.53,"tss_kgd":582.62,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":64262.91,"tss_kgd":150308.74,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37835.06,"tss_kgd":305.4,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":1065.59,"tss_kgd":6048.24,"is_gas":false,"phase":"l"},"ws180":{"id":"ws180","q_m3d":966.04,"tss_kgd":274.16,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":99.56,"tss_kgd":5774.32,"is_gas":false,"phase":"l"},"ws182":{"id":"ws182","q_m3d":99.56,"tss_kgd":5168.75,"is_gas":false,"phase":"l"},"biogas":{"id":"biogas","q_m3d":2.71,"tss_kgd":null,"is_gas":true,"phase":"g"},"digestate":{"id":"digestate","q_m3d":99.56,"tss_kgd":3878.24,"is_gas":false,"phase":"l"},"ws186":{"id":"ws186","q_m3d":99.56,"tss_kgd":3841.35,"is_gas":false,"phase":"l"},"ws187":{"id":"ws187","q_m3d":79.95,"tss_kgd":308.47,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":19.61,"tss_kgd":3533.02,"is_gas":false,"phase":"l"},"ws189":{"id":"ws189","q_m3d":1045.99,"tss_kgd":582.63,"is_gas":false,"phase":"l"}},"I2":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":25362.26,"tss_kgd":138321.88,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":1056.18,"tss_kgd":393.6,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":64272.56,"tss_kgd":144436.06,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37844.7,"tss_kgd":302.6,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":1065.59,"tss_kgd":5811.56,"is_gas":false,"phase":"l"},"ws191":{"id":"ws191","q_m3d":954.3,"tss_kgd":260.23,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":111.29,"tss_kgd":5551.3,"is_gas":false,"phase":"l"},"digestate":{"id":"digestate","q_m3d":111.29,"tss_kgd":1458.66,"is_gas":false,"phase":"l"},"ws194":{"id":"ws194","q_m3d":101.75,"tss_kgd":133.36,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":9.54,"tss_kgd":1325.42,"is_gas":false,"phase":"l"},"ws196":{"id":"ws196","q_m3d":1056.05,"tss_kgd":393.59,"is_gas":false,"phase":"l"}},"I3":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":25362.26,"tss_kgd":143740.78,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":1034.07,"tss_kgd":688.92,"is_gas":false,"phase":"l"},"treated":{"id":"treated","q_m3d":64250.45,"tss_kgd":150085.17,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37822.6,"tss_kgd":305.12,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":1065.59,"tss_kgd":6039.24,"is_gas":false,"phase":"l"},"ws198":{"id":"ws198","q_m3d":949.38,"tss_kgd":269.03,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":116.21,"tss_kgd":5770.12,"is_gas":false,"phase":"l"},"ws200":{"id":"ws200","q_m3d":84.57,"tss_kgd":419.91,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":31.65,"tss_kgd":5351.0,"is_gas":false,"phase":"l"},"ws202":{"id":"ws202","q_m3d":1033.95,"tss_kgd":688.93,"is_gas":false,"phase":"l"}},"N1":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":906.9,"tss_kgd":1926.03,"is_gas":false,"phase":"l"},"PE":{"id":"PE","q_m3d":38480.9,"tss_kgd":3251.12,"is_gas":false,"phase":"l"},"PS":{"id":"PS","q_m3d":280.12,"tss_kgd":4935.83,"is_gas":false,"phase":"l"},"carbon":{"id":"carbon","q_m3d":6.39,"tss_kgd":0.0,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":151416.48,"tss_kgd":1611048.16,"is_gas":false,"phase":"l"},"ws205":{"id":"ws205","q_m3d":189903.77,"tss_kgd":1618097.36,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37843.77,"tss_kgd":40.27,"is_gas":false,"phase":"l"},"ws206":{"id":"ws206","q_m3d":152060.0,"tss_kgd":1617895.1,"is_gas":false,"phase":"l"},"ws208":{"id":"ws208","q_m3d":151416.48,"tss_kgd":1611048.14,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":643.52,"tss_kgd":6846.95,"is_gas":false,"phase":"l"},"ws209":{"id":"ws209","q_m3d":743.83,"tss_kgd":1423.35,"is_gas":false,"phase":"l"},"thickened_sludge":{"id":"thickened_sludge","q_m3d":179.81,"tss_kgd":10359.61,"is_gas":false,"phase":"l"},"ws211":{"id":"ws211","q_m3d":179.81,"tss_kgd":9307.8,"is_gas":false,"phase":"l"},"biogas":{"id":"biogas","q_m3d":6.47,"tss_kgd":null,"is_gas":true,"phase":"g"},"digestate":{"id":"digestate","q_m3d":179.81,"tss_kgd":5994.02,"is_gas":false,"phase":"l"},"ws215":{"id":"ws215","q_m3d":179.81,"tss_kgd":5881.3,"is_gas":false,"phase":"l"},"ws216":{"id":"ws216","q_m3d":153.69,"tss_kgd":502.7,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":26.12,"tss_kgd":5378.66,"is_gas":false,"phase":"l"},"ws218":{"id":"ws218","q_m3d":897.52,"tss_kgd":1926.04,"is_gas":false,"phase":"l"}},"N2":{"RWW":{"id":"RWW","q_m3d":37854.12,"tss_kgd":6260.92,"is_gas":false,"phase":"l"},"carbon":{"id":"carbon","q_m3d":5.11,"tss_kgd":0.0,"is_gas":false,"phase":"l"},"RAS":{"id":"RAS","q_m3d":151416.48,"tss_kgd":1287802.72,"is_gas":false,"phase":"l"},"reject":{"id":"reject","q_m3d":1128.93,"tss_kgd":685.53,"is_gas":false,"phase":"l"},"ws219":{"id":"ws219","q_m3d":190404.64,"tss_kgd":1297552.13,"is_gas":false,"phase":"l"},"SE":{"id":"SE","q_m3d":37852.53,"tss_kgd":32.19,"is_gas":false,"phase":"l"},"ws220":{"id":"ws220","q_m3d":152552.1,"tss_kgd":1297461.21,"is_gas":false,"phase":"l"},"ws222":{"id":"ws222","q_m3d":151416.48,"tss_kgd":1287802.72,"is_gas":false,"phase":"l"},"WAS":{"id":"WAS","q_m3d":1135.62,"tss_kgd":9658.49,"is_gas":false,"phase":"l"},"ws223":{"id":"ws223","q_m3d":992.16,"tss_kgd":421.92,"is_gas":false,"phase":"l"},"thickened_WAS":{"id":"thickened_WAS","q_m3d":143.47,"tss_kgd":9236.79,"is_gas":false,"phase":"l"},"digestate":{"id":"digestate","q_m3d":143.47,"tss_kgd":2897.68,"is_gas":false,"phase":"l"},"ws226":{"id":"ws226","q_m3d":130.52,"tss_kgd":263.61,"is_gas":false,"phase":"l"},"cake":{"id":"cake","q_m3d":12.95,"tss_kgd":2634.8,"is_gas":false,"phase":"l"},"ws228":{"id":"ws228","q_m3d":1122.68,"tss_kgd":685.53,"is_gas":false,"phase":"l"}}};

   const LABEL_MAP = {
     'Raw WW': 'RWW', 'PE': 'PE', 'ML': 'treated', 'SE': 'SE',
     'Effluent': 'SE', 'Permeate': 'SE',
     'RAS (return activated sludge)': 'RAS', 'RAS': 'RAS',
     'MLR': 'MLR', 'MLR (mixed liquor recycle)': 'MLR',
     'PS': 'PS', 'WAS': 'WAS', 'sludge cake': 'cake', 'biogas': 'biogas',
     'reject water (centrate + thickener supernatant)': 'reject',
     'reject water': 'reject', 'permeate': 'SE',
     'FeCl3': 'FeCl3', 'alum': 'alum', 'settled effluent': 'settled_effluent',
   };

   const STREAM_DISPLAY = {
     'WAS': 'Waste activated sludge',
     'RAS': 'Return activated sludge',
   };

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
     if (/^O\d/.test(id))  return ['Aerobic Zone ' + id,   'Aerobic compartment of the plug-flow bioreactor (DO setpoint ~2 mg/L). Nitrification depends on SRT and DO profile.'];
     if (/^AN\d/.test(id)) return ['Anaerobic Zone ' + id, 'Anaerobic selector zone. Promotes phosphorus-accumulating organisms (PAOs) by creating feast conditions.'];
     if (/^A\d/.test(id))  return ['Anoxic Zone ' + id,    'Anoxic denitrification zone. Uses nitrate recycled from downstream aerobic zones (MLR) as electron acceptor.'];
     return null;
   }

   const STREAM_HOVER = {
     mainflow: { stroke: '#2055b8', sw: '2.2' },
     recycle:  { stroke: '#e03030', sw: '2'   },
     sludge:   { stroke: '#a07040', sw: '2'   },
     gas:      { stroke: '#8855cc', sw: '2'   },
   };
   const UNIT_SHAPE_RE = /clarifier|aerobic|anoxic|anaerobic|thickener|digester|mbr|chem|dw/;
   function unitHoverStyle(cls) {
     if (cls.includes('clarifier'))                               return {stroke:'#5588cc', sw:'2.2', br:'1.06'};
     if (cls.includes('digester-an') || cls.includes('anaerobic')) return {stroke:'#2a1060', sw:'2.2', br:'1.08'};
     if (cls.includes('digester-ae') || cls.includes('aerobic'))   return {stroke:'#1a4a80', sw:'2',   br:'1.08'};
     if (cls.includes('anoxic'))                                  return {stroke:'#1f5b50', sw:'2',   br:'1.08'};
     if (cls.includes('thickener'))                               return {stroke:'#5a4a20', sw:'2',   br:'1.06'};
     if (cls.includes('dw'))                                      return {stroke:'#4a3a10', sw:'2',   br:'1.06'};
     if (cls.includes('mbr'))                                     return {stroke:'#1a4e75', sw:'2.2', br:'1.08'};
     if (cls.includes('chem'))                                    return {stroke:'#a85a30', sw:'2',   br:'1.06'};
     return {stroke:'#5588cc', sw:'2', br:'1.05'};
   }

   const tt   = document.getElementById('wrrf-tt');
   const ttNm = document.getElementById('tt-nm');
   const ttBd = document.getElementById('tt-bd');

   function fmt(v, unit) {
     var n = parseFloat(v);
     if (isNaN(n)) return '-';
     return n.toLocaleString('en-US', {maximumFractionDigits: 0}) + ' ' + unit;
   }
   function showStream(e, label, sd) {
     tt.className = 'show tt-stream';
     ttNm.textContent = label;
     var q   = sd && sd.q_m3d   != null ? fmt(sd.q_m3d,   'm³/d') : '-';
     var tss = sd && sd.tss_kgd != null ? fmt(sd.tss_kgd, 'kg/d')
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
     var tw = tt.offsetWidth, th = tt.offsetHeight;
     var vw = window.innerWidth,  vh = window.innerHeight;
     tt.style.left = (e.clientX + 18 + tw > vw ? e.clientX - tw - 8 : e.clientX + 18) + 'px';
     tt.style.top  = (e.clientY + 18 + th > vh ? e.clientY - th - 8 : e.clientY + 18) + 'px';
   }
   function hide() { tt.classList.remove('show'); }

   function isInLegend(el, svg) {
     var p = el.parentElement;
     while (p && p !== svg) {
       if (p.tagName.toLowerCase() === 'g' && p.hasAttribute('transform')) return true;
       p = p.parentElement;
     }
     return false;
   }

   function midpointOf(el) {
     var tag = el.tagName.toLowerCase();
     if (tag === 'line') {
       return [(+el.getAttribute('x1') + +el.getAttribute('x2')) / 2,
               (+el.getAttribute('y1') + +el.getAttribute('y2')) / 2];
     }
     var nums = (el.getAttribute('d') || '').match(/[-\d.]+/g);
     return nums && nums.length >= 2 ? [+nums[0], +nums[1]] : [0, 0];
   }
   function dist2(a, b) { return (a[0]-b[0])**2 + (a[1]-b[1])**2; }

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

   function applyUnitHover(shapeEl, txtEl, on) {
     if (shapeEl) {
       var h = unitHoverStyle(shapeEl.getAttribute('class') || '');
       shapeEl.style.stroke      = on ? h.stroke                   : '';
       shapeEl.style.strokeWidth = on ? h.sw                       : '';
       shapeEl.style.filter      = on ? 'brightness(' + h.br + ')' : '';
     }
     if (txtEl) { txtEl.style.fontWeight = on ? '900' : ''; }
   }

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

     var flowEls = [];
     svg.querySelectorAll('line, path').forEach(function(el) {
       if (isInLegend(el, svg)) return;
       var cls = el.getAttribute('class') || '';
       var fc  = FLOW_CLASSES.find(function(c) { return cls.includes(c); });
       if (fc) flowEls.push({ el: el, fc: fc, mid: midpointOf(el) });
     });

     var labels = [];
     svg.querySelectorAll('text.streamlabel').forEach(function(el) {
       if (isInLegend(el, svg)) return;
       labels.push({ el: el, text: el.textContent.trim(),
                     pos: [+el.getAttribute('x') || 0, +el.getAttribute('y') || 0] });
     });

     var labelGroups = labels.map(function(lbl) {
       if (!flowEls.length) return { lbl: lbl, lines: [] };
       var dists = flowEls.map(function(f) { return { f: f, d: dist2(f.mid, lbl.pos) }; });
       dists.sort(function(a, b) { return a.d - b.d; });
       var threshold = Math.max(dists[0].d * 4, 2500);
       return { lbl: lbl, lines: dists.filter(function(x) { return x.d <= threshold; }).map(function(x) { return x.f; }) };
     });

     labelGroups.forEach(function(group) {
       var lbl   = group.lbl, lines = group.lines;
       var sid   = resolveStreamId(lbl.text, sysData);
       var sd    = sid ? sysData[sid] : null;
       var ttLabel = STREAM_DISPLAY[lbl.text] || lbl.text;
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
       } catch(err) {}
     });


   }

   document.addEventListener('DOMContentLoaded', function() {
     document.querySelectorAll('div.wrrf-diagram[id]').forEach(setupDiagram);
     document.addEventListener('mousemove', function(e) { if (tt.classList.contains('show')) move(e); });
   });
   </script>


.. raw:: html


   <div class="wrrf-diagram" id="B2">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">B2</text><text x="36" y="78" class="subtitle">rBOD · Primary clarifier · 6-zone aerobic ASP · Aerobic digestion</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="118" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><polygon points="118,205 182,205 175,261 125,261" class="clarifier"/><text x="150.0" y="238.0" class="zonelabel-d">PC</text><line x1="182" y1="233.0" x2="210" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="196.0" y="225.0" class="streamlabel" text-anchor="middle">PE</text><rect x="210" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="238.0" y="237.0" class="zonelabel">O1</text><circle cx="230.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="238.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="246.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="266" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="294.0" y="237.0" class="zonelabel">O2</text><circle cx="286.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="294.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="302.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="322" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="350.0" y="237.0" class="zonelabel">O3</text><circle cx="342.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="350.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="358.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="378" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="406.0" y="237.0" class="zonelabel">O4</text><circle cx="398.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="406.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="414.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="434" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="462.0" y="237.0" class="zonelabel">O5</text><circle cx="454.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="462.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="470.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="490" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="518.0" y="237.0" class="zonelabel">O6</text><circle cx="510.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="518.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="526.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><line x1="546" y1="233.0" x2="586" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="566.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="586,205 650,205 643,261 593,261" class="clarifier"/><text x="618.0" y="238.0" class="zonelabel-d">SC</text><line x1="650" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="955.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 600 261 L 600 320 L 218 320 L 218 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="409.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="150.0" y1="261" x2="150.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="142.0" y="275" class="streamlabel" text-anchor="end">PS</text><line x1="636" y1="261" x2="636" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="644" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="319.0,470 375.0,470 369.0,526 325.0,526" class="thickener"/><text x="347.0" y="503.0" class="zonelabel-d">GT</text><polygon points="485.0,470 541.0,470 535.0,526 491.0,526" class="thickener"/><text x="513.0" y="503.0" class="zonelabel-d">MT</text><rect x="651.0" y="470" width="64" height="56" class="digester-ae" rx="3"/><text x="683.0" y="503.0" class="zonelabel">AED</text><circle cx="671.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="679.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="687.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="695.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><rect x="825.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="853.0" y="503.0" class="zonelabel-d">DW</text><line x1="541.0" y1="498.0" x2="651.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><path d="M 375 498 L 375 456 L 651 456 L 651 492" class="sludge" marker-end="url(#arr-sludge)"/><line x1="715.0" y1="498.0" x2="825.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="150.0" y1="420" x2="347.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="347.0" y1="420" x2="347.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="636" y1="420" x2="513.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="513.0" y1="420" x2="513.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="881.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="1016.5" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="347.0" y1="526" x2="347.0" y2="600" class="recycle"/><line x1="513.0" y1="526" x2="513.0" y2="600" class="recycle"/><line x1="853.0" y1="526" x2="853.0" y2="600" class="recycle"/><line x1="347.0" y1="600" x2="853.0" y2="600" class="recycle"/><line x1="347.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="116" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="112" y1="255.0" x2="112" y2="236.0" class="recycle"/><line x1="112" y1="236.0" x2="116" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>

.. raw:: html


   <div class="wrrf-diagram" id="B3">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">B3</text><text x="36" y="78" class="subtitle">rBOD · Primary clarifier · 6-zone aerobic ASP · No digestion (thickening + dewatering only)</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="118" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><polygon points="118,205 182,205 175,261 125,261" class="clarifier"/><text x="150.0" y="238.0" class="zonelabel-d">PC</text><line x1="182" y1="233.0" x2="210" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="196.0" y="225.0" class="streamlabel" text-anchor="middle">PE</text><rect x="210" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="238.0" y="237.0" class="zonelabel">O1</text><circle cx="230.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="238.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="246.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="266" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="294.0" y="237.0" class="zonelabel">O2</text><circle cx="286.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="294.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="302.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="322" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="350.0" y="237.0" class="zonelabel">O3</text><circle cx="342.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="350.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="358.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="378" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="406.0" y="237.0" class="zonelabel">O4</text><circle cx="398.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="406.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="414.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="434" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="462.0" y="237.0" class="zonelabel">O5</text><circle cx="454.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="462.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="470.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="490" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="518.0" y="237.0" class="zonelabel">O6</text><circle cx="510.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="518.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="526.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><line x1="546" y1="233.0" x2="586" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="566.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="586,205 650,205 643,261 593,261" class="clarifier"/><text x="618.0" y="238.0" class="zonelabel-d">SC</text><line x1="650" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="955.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 600 261 L 600 320 L 218 320 L 218 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="409.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="150.0" y1="261" x2="150.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="142.0" y="275" class="streamlabel" text-anchor="end">PS</text><line x1="636" y1="261" x2="636" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="644" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="406.0,470 462.0,470 456.0,526 412.0,526" class="thickener"/><text x="434.0" y="503.0" class="zonelabel-d">GT</text><polygon points="572.0,470 628.0,470 622.0,526 578.0,526" class="thickener"/><text x="600.0" y="503.0" class="zonelabel-d">MT</text><rect x="738.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="766.0" y="503.0" class="zonelabel-d">DW</text><line x1="628.0" y1="498.0" x2="738.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><path d="M 462 498 L 462 456 L 738 456 L 738 492" class="sludge" marker-end="url(#arr-sludge)"/><line x1="150.0" y1="420" x2="434.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="434.0" y1="420" x2="434.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="636" y1="420" x2="600.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="600.0" y1="420" x2="600.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="794.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="973.0" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="434.0" y1="526" x2="434.0" y2="600" class="recycle"/><line x1="600.0" y1="526" x2="600.0" y2="600" class="recycle"/><line x1="766.0" y1="526" x2="766.0" y2="600" class="recycle"/><line x1="434.0" y1="600" x2="766.0" y2="600" class="recycle"/><line x1="434.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="116" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="112" y1="255.0" x2="112" y2="236.0" class="recycle"/><line x1="112" y1="236.0" x2="116" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>


C Series: rBOD without Primary Clarifier
----------------------------------------
Same six-zone aerobic activated sludge as the B series but without primary clarification — all influent organics enter the bioreactor directly, increasing aeration demand. The three C variants differ only in the solids train (AD, AED, or no digestion).

.. raw:: html


   <div class="wrrf-diagram" id="C1">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">C1</text><text x="36" y="78" class="subtitle">rBOD · No primary clarifier · 6-zone aerobic ASP · Anaerobic digestion</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="178" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><rect x="178" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="206.0" y="237.0" class="zonelabel">O1</text><circle cx="198.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="206.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="214.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="234" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="262.0" y="237.0" class="zonelabel">O2</text><circle cx="254.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="262.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="270.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="290" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="318.0" y="237.0" class="zonelabel">O3</text><circle cx="310.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="318.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="326.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="346" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="374.0" y="237.0" class="zonelabel">O4</text><circle cx="366.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="374.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="382.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="402" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="430.0" y="237.0" class="zonelabel">O5</text><circle cx="422.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="430.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="438.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="458" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="486.0" y="237.0" class="zonelabel">O6</text><circle cx="478.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="486.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="494.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><line x1="514" y1="233.0" x2="554" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="534.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="554,205 618,205 611,261 561,261" class="clarifier"/><text x="586.0" y="238.0" class="zonelabel-d">SC</text><line x1="618" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="939.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 568 261 L 568 320 L 186 320 L 186 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="377.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="604" y1="261" x2="604" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="612" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="402.0,470 458.0,470 452.0,526 408.0,526" class="thickener"/><text x="430.0" y="503.0" class="zonelabel-d">MT</text><path d="M568.0,480 Q568.0,470 600.0,470 Q632.0,470 632.0,480 L632.0,526 L568.0,526 z" class="digester-an"/><text x="600.0" y="503.0" class="zonelabel">AD</text><line x1="600.0" y1="468" x2="600.0" y2="448" class="gas" marker-end="url(#arr-gas)"/><text x="605.0" y="445" class="streamlabel" fill="#6b4a8a">biogas</text><rect x="742.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="770.0" y="503.0" class="zonelabel-d">DW</text><line x1="458.0" y1="498.0" x2="568.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="632.0" y1="498.0" x2="742.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="604" y1="420" x2="430.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="430.0" y1="420" x2="430.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="798.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="975.0" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="430.0" y1="526" x2="430.0" y2="600" class="recycle"/><line x1="770.0" y1="526" x2="770.0" y2="600" class="recycle"/><line x1="430.0" y1="600" x2="770.0" y2="600" class="recycle"/><line x1="430.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="176" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="172" y1="255.0" x2="172" y2="236.0" class="recycle"/><line x1="172" y1="236.0" x2="176" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>

.. raw:: html


   <div class="wrrf-diagram" id="C2">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">C2</text><text x="36" y="78" class="subtitle">rBOD · No primary clarifier · 6-zone aerobic ASP · Aerobic digestion</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="178" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><rect x="178" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="206.0" y="237.0" class="zonelabel">O1</text><circle cx="198.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="206.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="214.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="234" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="262.0" y="237.0" class="zonelabel">O2</text><circle cx="254.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="262.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="270.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="290" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="318.0" y="237.0" class="zonelabel">O3</text><circle cx="310.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="318.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="326.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="346" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="374.0" y="237.0" class="zonelabel">O4</text><circle cx="366.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="374.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="382.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="402" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="430.0" y="237.0" class="zonelabel">O5</text><circle cx="422.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="430.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="438.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="458" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="486.0" y="237.0" class="zonelabel">O6</text><circle cx="478.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="486.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="494.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><line x1="514" y1="233.0" x2="554" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="534.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="554,205 618,205 611,261 561,261" class="clarifier"/><text x="586.0" y="238.0" class="zonelabel-d">SC</text><line x1="618" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="939.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 568 261 L 568 320 L 186 320 L 186 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="377.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="604" y1="261" x2="604" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="612" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="402.0,470 458.0,470 452.0,526 408.0,526" class="thickener"/><text x="430.0" y="503.0" class="zonelabel-d">MT</text><rect x="568.0" y="470" width="64" height="56" class="digester-ae" rx="3"/><text x="600.0" y="503.0" class="zonelabel">AED</text><circle cx="588.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="596.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="604.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="612.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><rect x="742.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="770.0" y="503.0" class="zonelabel-d">DW</text><line x1="458.0" y1="498.0" x2="568.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="632.0" y1="498.0" x2="742.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="604" y1="420" x2="430.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="430.0" y1="420" x2="430.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="798.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="975.0" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="430.0" y1="526" x2="430.0" y2="600" class="recycle"/><line x1="770.0" y1="526" x2="770.0" y2="600" class="recycle"/><line x1="430.0" y1="600" x2="770.0" y2="600" class="recycle"/><line x1="430.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="176" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="172" y1="255.0" x2="172" y2="236.0" class="recycle"/><line x1="172" y1="236.0" x2="176" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>

.. raw:: html


   <div class="wrrf-diagram" id="C3">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">C3</text><text x="36" y="78" class="subtitle">rBOD · No primary clarifier · 6-zone aerobic ASP · No digestion (thickening + dewatering only)</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="178" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><rect x="178" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="206.0" y="237.0" class="zonelabel">O1</text><circle cx="198.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="206.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="214.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="234" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="262.0" y="237.0" class="zonelabel">O2</text><circle cx="254.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="262.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="270.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="290" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="318.0" y="237.0" class="zonelabel">O3</text><circle cx="310.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="318.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="326.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="346" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="374.0" y="237.0" class="zonelabel">O4</text><circle cx="366.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="374.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="382.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="402" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="430.0" y="237.0" class="zonelabel">O5</text><circle cx="422.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="430.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="438.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="458" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="486.0" y="237.0" class="zonelabel">O6</text><circle cx="478.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="486.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="494.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><line x1="514" y1="233.0" x2="554" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="534.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="554,205 618,205 611,261 561,261" class="clarifier"/><text x="586.0" y="238.0" class="zonelabel-d">SC</text><line x1="618" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="939.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 568 261 L 568 320 L 186 320 L 186 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="377.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="604" y1="261" x2="604" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="612" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="489.0,470 545.0,470 539.0,526 495.0,526" class="thickener"/><text x="517.0" y="503.0" class="zonelabel-d">MT</text><rect x="655.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="683.0" y="503.0" class="zonelabel-d">DW</text><line x1="545.0" y1="498.0" x2="655.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="604" y1="420" x2="517.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="517.0" y1="420" x2="517.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="711.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="931.5" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="517.0" y1="526" x2="517.0" y2="600" class="recycle"/><line x1="683.0" y1="526" x2="683.0" y2="600" class="recycle"/><line x1="517.0" y1="600" x2="683.0" y2="600" class="recycle"/><line x1="517.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="176" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="172" y1="255.0" x2="172" y2="236.0" class="recycle"/><line x1="172" y1="236.0" x2="176" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>


E Series: NIT with Aerobic Digestion
------------------------------------
Longer SRT (~6 d) retains nitrifiers for NH₄⁺ oxidation. E2-P adds primary clarification to the E2 liquid train.

.. raw:: html


   <div class="wrrf-diagram" id="E2">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">E2</text><text x="36" y="78" class="subtitle">NIT · No primary clarifier · Nitrifying ASP · Aerobic digestion</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="178" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><rect x="178" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="206.0" y="237.0" class="zonelabel">O1</text><circle cx="198.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="206.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="214.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="234" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="262.0" y="237.0" class="zonelabel">O2</text><circle cx="254.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="262.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="270.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="290" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="318.0" y="237.0" class="zonelabel">O3</text><circle cx="310.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="318.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="326.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="346" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="374.0" y="237.0" class="zonelabel">O4</text><circle cx="366.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="374.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="382.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="402" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="430.0" y="237.0" class="zonelabel">O5</text><circle cx="422.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="430.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="438.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="458" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="486.0" y="237.0" class="zonelabel">O6</text><circle cx="478.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="486.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="494.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><line x1="514" y1="233.0" x2="554" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="534.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="554,205 618,205 611,261 561,261" class="clarifier"/><text x="586.0" y="238.0" class="zonelabel-d">SC</text><line x1="618" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="939.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 568 261 L 568 320 L 186 320 L 186 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="377.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="604" y1="261" x2="604" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="612" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="402.0,470 458.0,470 452.0,526 408.0,526" class="thickener"/><text x="430.0" y="503.0" class="zonelabel-d">MT</text><rect x="568.0" y="470" width="64" height="56" class="digester-ae" rx="3"/><text x="600.0" y="503.0" class="zonelabel">AED</text><circle cx="588.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="596.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="604.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="612.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><rect x="742.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="770.0" y="503.0" class="zonelabel-d">DW</text><line x1="458.0" y1="498.0" x2="568.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="632.0" y1="498.0" x2="742.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="604" y1="420" x2="430.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="430.0" y1="420" x2="430.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="798.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="975.0" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="430.0" y1="526" x2="430.0" y2="600" class="recycle"/><line x1="770.0" y1="526" x2="770.0" y2="600" class="recycle"/><line x1="430.0" y1="600" x2="770.0" y2="600" class="recycle"/><line x1="430.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="176" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="172" y1="255.0" x2="172" y2="236.0" class="recycle"/><line x1="172" y1="236.0" x2="176" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>

.. raw:: html


   <div class="wrrf-diagram" id="E2P">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">E2-P</text><text x="36" y="78" class="subtitle">NIT · Primary clarifier ("prime" variant) · Nitrifying ASP · Aerobic digestion</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="118" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><polygon points="118,205 182,205 175,261 125,261" class="clarifier"/><text x="150.0" y="238.0" class="zonelabel-d">PC</text><line x1="182" y1="233.0" x2="210" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="196.0" y="225.0" class="streamlabel" text-anchor="middle">PE</text><rect x="210" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="238.0" y="237.0" class="zonelabel">O1</text><circle cx="230.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="238.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="246.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="266" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="294.0" y="237.0" class="zonelabel">O2</text><circle cx="286.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="294.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="302.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="322" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="350.0" y="237.0" class="zonelabel">O3</text><circle cx="342.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="350.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="358.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="378" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="406.0" y="237.0" class="zonelabel">O4</text><circle cx="398.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="406.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="414.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="434" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="462.0" y="237.0" class="zonelabel">O5</text><circle cx="454.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="462.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="470.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="490" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="518.0" y="237.0" class="zonelabel">O6</text><circle cx="510.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="518.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="526.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><line x1="546" y1="233.0" x2="586" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="566.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="586,205 650,205 643,261 593,261" class="clarifier"/><text x="618.0" y="238.0" class="zonelabel-d">SC</text><line x1="650" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="955.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 600 261 L 600 320 L 218 320 L 218 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="409.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="150.0" y1="261" x2="150.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="142.0" y="275" class="streamlabel" text-anchor="end">PS</text><line x1="636" y1="261" x2="636" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="644" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="319.0,470 375.0,470 369.0,526 325.0,526" class="thickener"/><text x="347.0" y="503.0" class="zonelabel-d">GT</text><polygon points="485.0,470 541.0,470 535.0,526 491.0,526" class="thickener"/><text x="513.0" y="503.0" class="zonelabel-d">MT</text><rect x="651.0" y="470" width="64" height="56" class="digester-ae" rx="3"/><text x="683.0" y="503.0" class="zonelabel">AED</text><circle cx="671.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="679.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="687.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="695.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><rect x="825.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="853.0" y="503.0" class="zonelabel-d">DW</text><line x1="541.0" y1="498.0" x2="651.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><path d="M 375 498 L 375 456 L 651 456 L 651 492" class="sludge" marker-end="url(#arr-sludge)"/><line x1="715.0" y1="498.0" x2="825.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="150.0" y1="420" x2="347.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="347.0" y1="420" x2="347.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="636" y1="420" x2="513.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="513.0" y1="420" x2="513.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="881.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="1016.5" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="347.0" y1="526" x2="347.0" y2="600" class="recycle"/><line x1="513.0" y1="526" x2="513.0" y2="600" class="recycle"/><line x1="853.0" y1="526" x2="853.0" y2="600" class="recycle"/><line x1="347.0" y1="600" x2="853.0" y2="600" class="recycle"/><line x1="347.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="116" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="112" y1="255.0" x2="112" y2="236.0" class="recycle"/><line x1="112" y1="236.0" x2="116" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>


F Series: NIT with Anaerobic Digestion
--------------------------------------
Compared to E2-P, F1 differs in solids treatment by using anaerobic digestion.

.. raw:: html


   <div class="wrrf-diagram" id="F1">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">F1</text><text x="36" y="78" class="subtitle">NIT · Primary clarifier · Nitrifying ASP · Anaerobic digestion</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="118" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><polygon points="118,205 182,205 175,261 125,261" class="clarifier"/><text x="150.0" y="238.0" class="zonelabel-d">PC</text><line x1="182" y1="233.0" x2="210" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="196.0" y="225.0" class="streamlabel" text-anchor="middle">PE</text><rect x="210" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="238.0" y="237.0" class="zonelabel">O1</text><circle cx="230.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="238.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="246.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="266" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="294.0" y="237.0" class="zonelabel">O2</text><circle cx="286.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="294.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="302.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="322" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="350.0" y="237.0" class="zonelabel">O3</text><circle cx="342.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="350.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="358.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="378" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="406.0" y="237.0" class="zonelabel">O4</text><circle cx="398.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="406.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="414.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="434" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="462.0" y="237.0" class="zonelabel">O5</text><circle cx="454.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="462.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="470.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="490" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="518.0" y="237.0" class="zonelabel">O6</text><circle cx="510.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="518.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="526.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><line x1="546" y1="233.0" x2="586" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="566.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="586,205 650,205 643,261 593,261" class="clarifier"/><text x="618.0" y="238.0" class="zonelabel-d">SC</text><line x1="650" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="955.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 600 261 L 600 320 L 218 320 L 218 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="409.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="150.0" y1="261" x2="150.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="142.0" y="275" class="streamlabel" text-anchor="end">PS</text><line x1="636" y1="261" x2="636" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="644" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="319.0,470 375.0,470 369.0,526 325.0,526" class="thickener"/><text x="347.0" y="503.0" class="zonelabel-d">GT</text><polygon points="485.0,470 541.0,470 535.0,526 491.0,526" class="thickener"/><text x="513.0" y="503.0" class="zonelabel-d">MT</text><path d="M651.0,480 Q651.0,470 683.0,470 Q715.0,470 715.0,480 L715.0,526 L651.0,526 z" class="digester-an"/><text x="683.0" y="503.0" class="zonelabel">AD</text><line x1="683.0" y1="468" x2="683.0" y2="448" class="gas" marker-end="url(#arr-gas)"/><text x="688.0" y="445" class="streamlabel" fill="#6b4a8a">biogas</text><rect x="825.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="853.0" y="503.0" class="zonelabel-d">DW</text><line x1="541.0" y1="498.0" x2="651.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><path d="M 375 498 L 375 456 L 651 456 L 651 492" class="sludge" marker-end="url(#arr-sludge)"/><line x1="715.0" y1="498.0" x2="825.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="150.0" y1="420" x2="347.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="347.0" y1="420" x2="347.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="636" y1="420" x2="513.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="513.0" y1="420" x2="513.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="881.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="1016.5" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="347.0" y1="526" x2="347.0" y2="600" class="recycle"/><line x1="513.0" y1="526" x2="513.0" y2="600" class="recycle"/><line x1="853.0" y1="526" x2="853.0" y2="600" class="recycle"/><line x1="347.0" y1="600" x2="853.0" y2="600" class="recycle"/><line x1="347.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="116" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="112" y1="255.0" x2="112" y2="236.0" class="recycle"/><line x1="112" y1="236.0" x2="116" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>


G Series: BNR with Johannesburg Process
---------------------------------------
BNR configurations with primary clarification. G1–G3 use the Johannesburg step-fed process with external carbon for EBPR; 

.. raw:: html


   <div class="wrrf-diagram" id="G1">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">G1</text><text x="36" y="78" class="subtitle">BNR · Primary clarifier · Johannesburg process (step-fed) · Anaerobic digestion</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="118" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><polygon points="118,205 182,205 175,261 125,261" class="clarifier"/><text x="150.0" y="238.0" class="zonelabel-d">PC</text><line x1="182" y1="233.0" x2="202" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><rect x="210" y="205" width="56" height="56" class="anaerobic" rx="2"/><text x="238.0" y="237.0" class="zonelabel">A1</text><rect x="266" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="294.0" y="237.0" class="zonelabel">A2</text><rect x="322" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="350.0" y="237.0" class="zonelabel">A3</text><rect x="378" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="406.0" y="237.0" class="zonelabel">A4</text><rect x="434" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="462.0" y="237.0" class="zonelabel">O5</text><circle cx="454.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="462.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="470.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="490" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="518.0" y="237.0" class="zonelabel">O6</text><circle cx="510.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="518.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="526.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="202" cy="233.0" r="3.5" fill="#1a1a1a"/><path d="M 202 233.0 L 202 294 L 294.0 294 L 294.0 261" class="mainflow" marker-end="url(#arr-main)"/><text x="248.0" y="289" class="streamlabel" text-anchor="middle">80%</text><path d="M 202 233.0 L 202 178 L 350.0 178 L 350.0 205" class="mainflow" marker-end="url(#arr-main)"/><text x="276.0" y="173" class="streamlabel" text-anchor="middle">20%</text><rect x="178.0" y="128" width="120" height="30" class="chem" rx="15"/><line x1="238.0" y1="158" x2="238.0" y2="204" class="dose" marker-end="url(#arr-dose)"/><text x="238" y="143" class="zonelabel" dominant-baseline="central">External C</text><path d="M 532.0 205 L 532.0 145 L 336.0 145 L 336.0 205" class="recycle" marker-end="url(#arr-recycle)"/><text x="434.0" y="140" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text><line x1="546" y1="233.0" x2="586" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="566.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="586,205 650,205 643,261 593,261" class="clarifier"/><text x="618.0" y="238.0" class="zonelabel-d">SC</text><line x1="650" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="955.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 600 261 L 600 320 L 218 320 L 218 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="409.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="150.0" y1="261" x2="150.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="142.0" y="275" class="streamlabel" text-anchor="end">PS</text><line x1="636" y1="261" x2="636" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="644" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="319.0,470 375.0,470 369.0,526 325.0,526" class="thickener"/><text x="347.0" y="503.0" class="zonelabel-d">GT</text><polygon points="485.0,470 541.0,470 535.0,526 491.0,526" class="thickener"/><text x="513.0" y="503.0" class="zonelabel-d">MT</text><path d="M651.0,480 Q651.0,470 683.0,470 Q715.0,470 715.0,480 L715.0,526 L651.0,526 z" class="digester-an"/><text x="683.0" y="503.0" class="zonelabel">AD</text><line x1="683.0" y1="468" x2="683.0" y2="448" class="gas" marker-end="url(#arr-gas)"/><text x="688.0" y="445" class="streamlabel" fill="#6b4a8a">biogas</text><rect x="825.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="853.0" y="503.0" class="zonelabel-d">DW</text><line x1="541.0" y1="498.0" x2="651.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><path d="M 375 498 L 375 456 L 651 456 L 651 492" class="sludge" marker-end="url(#arr-sludge)"/><line x1="715.0" y1="498.0" x2="825.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="150.0" y1="420" x2="347.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="347.0" y1="420" x2="347.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="636" y1="420" x2="513.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="513.0" y1="420" x2="513.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="881.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="1016.5" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="347.0" y1="526" x2="347.0" y2="600" class="recycle"/><line x1="513.0" y1="526" x2="513.0" y2="600" class="recycle"/><line x1="853.0" y1="526" x2="853.0" y2="600" class="recycle"/><line x1="347.0" y1="600" x2="853.0" y2="600" class="recycle"/><line x1="347.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="116" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="112" y1="255.0" x2="112" y2="236.0" class="recycle"/><line x1="112" y1="236.0" x2="116" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>

.. raw:: html


   <div class="wrrf-diagram" id="G2">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">G2</text><text x="36" y="78" class="subtitle">BNR · Primary clarifier · Johannesburg process (step-fed) · Aerobic digestion</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="118" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><polygon points="118,205 182,205 175,261 125,261" class="clarifier"/><text x="150.0" y="238.0" class="zonelabel-d">PC</text><line x1="182" y1="233.0" x2="202" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><rect x="210" y="205" width="56" height="56" class="anaerobic" rx="2"/><text x="238.0" y="237.0" class="zonelabel">A1</text><rect x="266" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="294.0" y="237.0" class="zonelabel">A2</text><rect x="322" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="350.0" y="237.0" class="zonelabel">A3</text><rect x="378" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="406.0" y="237.0" class="zonelabel">A4</text><rect x="434" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="462.0" y="237.0" class="zonelabel">O5</text><circle cx="454.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="462.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="470.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="490" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="518.0" y="237.0" class="zonelabel">O6</text><circle cx="510.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="518.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="526.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="202" cy="233.0" r="3.5" fill="#1a1a1a"/><path d="M 202 233.0 L 202 294 L 294.0 294 L 294.0 261" class="mainflow" marker-end="url(#arr-main)"/><text x="248.0" y="289" class="streamlabel" text-anchor="middle">80%</text><path d="M 202 233.0 L 202 178 L 350.0 178 L 350.0 205" class="mainflow" marker-end="url(#arr-main)"/><text x="276.0" y="173" class="streamlabel" text-anchor="middle">20%</text><rect x="178.0" y="128" width="120" height="30" class="chem" rx="15"/><line x1="238.0" y1="158" x2="238.0" y2="204" class="dose" marker-end="url(#arr-dose)"/><text x="238" y="143" class="zonelabel" dominant-baseline="central">External C</text><path d="M 532.0 205 L 532.0 145 L 336.0 145 L 336.0 205" class="recycle" marker-end="url(#arr-recycle)"/><text x="434.0" y="140" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text><line x1="546" y1="233.0" x2="586" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="566.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="586,205 650,205 643,261 593,261" class="clarifier"/><text x="618.0" y="238.0" class="zonelabel-d">SC</text><line x1="650" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="955.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 600 261 L 600 320 L 218 320 L 218 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="409.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="150.0" y1="261" x2="150.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="142.0" y="275" class="streamlabel" text-anchor="end">PS</text><line x1="636" y1="261" x2="636" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="644" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="319.0,470 375.0,470 369.0,526 325.0,526" class="thickener"/><text x="347.0" y="503.0" class="zonelabel-d">GT</text><polygon points="485.0,470 541.0,470 535.0,526 491.0,526" class="thickener"/><text x="513.0" y="503.0" class="zonelabel-d">MT</text><rect x="651.0" y="470" width="64" height="56" class="digester-ae" rx="3"/><text x="683.0" y="503.0" class="zonelabel">AED</text><circle cx="671.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="679.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="687.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="695.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><rect x="825.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="853.0" y="503.0" class="zonelabel-d">DW</text><line x1="541.0" y1="498.0" x2="651.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><path d="M 375 498 L 375 456 L 651 456 L 651 492" class="sludge" marker-end="url(#arr-sludge)"/><line x1="715.0" y1="498.0" x2="825.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="150.0" y1="420" x2="347.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="347.0" y1="420" x2="347.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="636" y1="420" x2="513.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="513.0" y1="420" x2="513.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="881.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="1016.5" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="347.0" y1="526" x2="347.0" y2="600" class="recycle"/><line x1="513.0" y1="526" x2="513.0" y2="600" class="recycle"/><line x1="853.0" y1="526" x2="853.0" y2="600" class="recycle"/><line x1="347.0" y1="600" x2="853.0" y2="600" class="recycle"/><line x1="347.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="116" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="112" y1="255.0" x2="112" y2="236.0" class="recycle"/><line x1="112" y1="236.0" x2="116" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>

.. raw:: html


   <div class="wrrf-diagram" id="G3">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">G3</text><text x="36" y="78" class="subtitle">BNR · Primary clarifier · Johannesburg process (step-fed) · No digestion (thickening + dewatering only)</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="118" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><polygon points="118,205 182,205 175,261 125,261" class="clarifier"/><text x="150.0" y="238.0" class="zonelabel-d">PC</text><line x1="182" y1="233.0" x2="202" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><rect x="210" y="205" width="56" height="56" class="anaerobic" rx="2"/><text x="238.0" y="237.0" class="zonelabel">A1</text><rect x="266" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="294.0" y="237.0" class="zonelabel">A2</text><rect x="322" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="350.0" y="237.0" class="zonelabel">A3</text><rect x="378" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="406.0" y="237.0" class="zonelabel">A4</text><rect x="434" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="462.0" y="237.0" class="zonelabel">O5</text><circle cx="454.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="462.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="470.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="490" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="518.0" y="237.0" class="zonelabel">O6</text><circle cx="510.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="518.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="526.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="202" cy="233.0" r="3.5" fill="#1a1a1a"/><path d="M 202 233.0 L 202 294 L 294.0 294 L 294.0 261" class="mainflow" marker-end="url(#arr-main)"/><text x="248.0" y="289" class="streamlabel" text-anchor="middle">80%</text><path d="M 202 233.0 L 202 178 L 350.0 178 L 350.0 205" class="mainflow" marker-end="url(#arr-main)"/><text x="276.0" y="173" class="streamlabel" text-anchor="middle">20%</text><rect x="178.0" y="128" width="120" height="30" class="chem" rx="15"/><line x1="238.0" y1="158" x2="238.0" y2="204" class="dose" marker-end="url(#arr-dose)"/><text x="238" y="143" class="zonelabel" dominant-baseline="central">External C</text><path d="M 532.0 205 L 532.0 145 L 336.0 145 L 336.0 205" class="recycle" marker-end="url(#arr-recycle)"/><text x="434.0" y="140" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text><line x1="546" y1="233.0" x2="586" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="566.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="586,205 650,205 643,261 593,261" class="clarifier"/><text x="618.0" y="238.0" class="zonelabel-d">SC</text><line x1="650" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="955.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 600 261 L 600 320 L 218 320 L 218 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="409.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="150.0" y1="261" x2="150.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="142.0" y="275" class="streamlabel" text-anchor="end">PS</text><line x1="636" y1="261" x2="636" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="644" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="406.0,470 462.0,470 456.0,526 412.0,526" class="thickener"/><text x="434.0" y="503.0" class="zonelabel-d">GT</text><polygon points="572.0,470 628.0,470 622.0,526 578.0,526" class="thickener"/><text x="600.0" y="503.0" class="zonelabel-d">MT</text><rect x="738.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="766.0" y="503.0" class="zonelabel-d">DW</text><line x1="628.0" y1="498.0" x2="738.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><path d="M 462 498 L 462 456 L 738 456 L 738 492" class="sludge" marker-end="url(#arr-sludge)"/><line x1="150.0" y1="420" x2="434.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="434.0" y1="420" x2="434.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="636" y1="420" x2="600.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="600.0" y1="420" x2="600.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="794.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="973.0" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="434.0" y1="526" x2="434.0" y2="600" class="recycle"/><line x1="600.0" y1="526" x2="600.0" y2="600" class="recycle"/><line x1="766.0" y1="526" x2="766.0" y2="600" class="recycle"/><line x1="434.0" y1="600" x2="766.0" y2="600" class="recycle"/><line x1="434.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="116" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="112" y1="255.0" x2="112" y2="236.0" class="recycle"/><line x1="112" y1="236.0" x2="116" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>


H Series: BNR with Chemical P Removal
-------------------------------------
H1 combines chemical phosphorus removal (Fe:sup:`3+`/Al:sup:`3+`) with a 4-stage Bardenpho bioreactor and external carbon.

.. raw:: html


   <div class="wrrf-diagram" id="H1">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">H1</text><text x="36" y="78" class="subtitle">BNR · Primary clarifier + chemical P removal · Modified 5-stage Bardenpho · Anaerobic digestion</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="118" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><polygon points="118,205 182,205 175,261 125,261" class="clarifier"/><text x="150.0" y="238.0" class="zonelabel-d">PC</text><rect x="29.0" y="113" width="120" height="30" class="chem" rx="15"/><line x1="89.0" y1="143" x2="89.0" y2="232.0" class="dose" marker-end="url(#arr-dose)"/><text x="89" y="128" class="zonelabel" dominant-baseline="central">Coagulant</text><line x1="182" y1="233.0" x2="210" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="196.0" y="225.0" class="streamlabel" text-anchor="middle">PE</text><rect x="210" y="205" width="56" height="56" class="anaerobic" rx="2"/><text x="238.0" y="237.0" class="zonelabel">Anae</text><rect x="266" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="294.0" y="237.0" class="zonelabel">Anox</text><rect x="322" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="350.0" y="237.0" class="zonelabel">Aer</text><circle cx="342.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="350.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="358.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="378" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="406.0" y="237.0" class="zonelabel">Aer</text><circle cx="398.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="406.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="414.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="434" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="462.0" y="237.0" class="zonelabel">Anox</text><rect x="490" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="518.0" y="237.0" class="zonelabel">Aer</text><circle cx="510.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="518.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="526.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="178.0" y="113" width="120" height="30" class="chem" rx="15"/><line x1="238.0" y1="143" x2="238.0" y2="204" class="dose" marker-end="url(#arr-dose)"/><text x="238" y="128" class="zonelabel" dominant-baseline="central">External C</text><path d="M 420.0 205 L 420.0 145 L 224.0 145 L 224.0 205" class="recycle" marker-end="url(#arr-recycle)"/><text x="348.0" y="170" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text><line x1="546" y1="233.0" x2="586" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="566.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="586,205 650,205 643,261 593,261" class="clarifier"/><text x="618.0" y="238.0" class="zonelabel-d">SC</text><line x1="650" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="955.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 600 261 L 600 320 L 218 320 L 218 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="409.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="150.0" y1="261" x2="150.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="142.0" y="275" class="streamlabel" text-anchor="end">PS</text><line x1="636" y1="261" x2="636" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="644" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="319.0,470 375.0,470 369.0,526 325.0,526" class="thickener"/><text x="347.0" y="503.0" class="zonelabel-d">GT</text><polygon points="485.0,470 541.0,470 535.0,526 491.0,526" class="thickener"/><text x="513.0" y="503.0" class="zonelabel-d">MT</text><path d="M651.0,480 Q651.0,470 683.0,470 Q715.0,470 715.0,480 L715.0,526 L651.0,526 z" class="digester-an"/><text x="683.0" y="503.0" class="zonelabel">AD</text><line x1="683.0" y1="468" x2="683.0" y2="448" class="gas" marker-end="url(#arr-gas)"/><text x="688.0" y="445" class="streamlabel" fill="#6b4a8a">biogas</text><rect x="825.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="853.0" y="503.0" class="zonelabel-d">DW</text><line x1="541.0" y1="498.0" x2="651.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><path d="M 375 498 L 375 456 L 651 456 L 651 492" class="sludge" marker-end="url(#arr-sludge)"/><line x1="715.0" y1="498.0" x2="825.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="150.0" y1="420" x2="347.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="347.0" y1="420" x2="347.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="636" y1="420" x2="513.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="513.0" y1="420" x2="513.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="881.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="1016.5" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="347.0" y1="526" x2="347.0" y2="600" class="recycle"/><line x1="513.0" y1="526" x2="513.0" y2="600" class="recycle"/><line x1="853.0" y1="526" x2="853.0" y2="600" class="recycle"/><line x1="347.0" y1="600" x2="853.0" y2="600" class="recycle"/><line x1="347.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="116" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="112" y1="255.0" x2="112" y2="236.0" class="recycle"/><line x1="112" y1="236.0" x2="116" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>


I Series: BNR with 5-stage Bardenpho
------------------------------------
Five-stage Bardenpho without primary clarification — raw wastewater provides carbon for denitrification and EBPR without external dosing. The three I variants differ only in the solids train (AD, AED, or no digestion).

.. raw:: html


   <div class="wrrf-diagram" id="I1">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">I1</text><text x="36" y="78" class="subtitle">BNR · No primary clarifier · 5-stage Bardenpho · Anaerobic digestion</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="178" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><rect x="178" y="205" width="56" height="56" class="anaerobic" rx="2"/><text x="206.0" y="237.0" class="zonelabel">Anae</text><rect x="234" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="262.0" y="237.0" class="zonelabel">Anox</text><rect x="290" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="318.0" y="237.0" class="zonelabel">Aer</text><circle cx="310.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="318.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="326.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="346" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="374.0" y="237.0" class="zonelabel">Aer</text><circle cx="366.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="374.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="382.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="402" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="430.0" y="237.0" class="zonelabel">Anox</text><rect x="458" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="486.0" y="237.0" class="zonelabel">Aer</text><circle cx="478.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="486.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="494.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><path d="M 388.0 205 L 388.0 145 L 248.0 145 L 248.0 205" class="recycle" marker-end="url(#arr-recycle)"/><text x="318.0" y="140" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text><line x1="514" y1="233.0" x2="554" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="534.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="554,205 618,205 611,261 561,261" class="clarifier"/><text x="586.0" y="238.0" class="zonelabel-d">SC</text><line x1="618" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="939.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 568 261 L 568 320 L 186 320 L 186 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="377.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="604" y1="261" x2="604" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="612" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="402.0,470 458.0,470 452.0,526 408.0,526" class="thickener"/><text x="430.0" y="503.0" class="zonelabel-d">MT</text><path d="M568.0,480 Q568.0,470 600.0,470 Q632.0,470 632.0,480 L632.0,526 L568.0,526 z" class="digester-an"/><text x="600.0" y="503.0" class="zonelabel">AD</text><line x1="600.0" y1="468" x2="600.0" y2="448" class="gas" marker-end="url(#arr-gas)"/><text x="605.0" y="445" class="streamlabel" fill="#6b4a8a">biogas</text><rect x="742.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="770.0" y="503.0" class="zonelabel-d">DW</text><line x1="458.0" y1="498.0" x2="568.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="632.0" y1="498.0" x2="742.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="604" y1="420" x2="430.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="430.0" y1="420" x2="430.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="798.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="975.0" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="430.0" y1="526" x2="430.0" y2="600" class="recycle"/><line x1="770.0" y1="526" x2="770.0" y2="600" class="recycle"/><line x1="430.0" y1="600" x2="770.0" y2="600" class="recycle"/><line x1="430.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="176" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="172" y1="255.0" x2="172" y2="236.0" class="recycle"/><line x1="172" y1="236.0" x2="176" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>

.. raw:: html


   <div class="wrrf-diagram" id="I2">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">I2</text><text x="36" y="78" class="subtitle">BNR · No primary clarifier · 5-stage Bardenpho · Aerobic digestion</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="178" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><rect x="178" y="205" width="56" height="56" class="anaerobic" rx="2"/><text x="206.0" y="237.0" class="zonelabel">Anae</text><rect x="234" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="262.0" y="237.0" class="zonelabel">Anox</text><rect x="290" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="318.0" y="237.0" class="zonelabel">Aer</text><circle cx="310.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="318.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="326.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="346" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="374.0" y="237.0" class="zonelabel">Aer</text><circle cx="366.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="374.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="382.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="402" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="430.0" y="237.0" class="zonelabel">Anox</text><rect x="458" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="486.0" y="237.0" class="zonelabel">Aer</text><circle cx="478.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="486.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="494.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><path d="M 388.0 205 L 388.0 145 L 248.0 145 L 248.0 205" class="recycle" marker-end="url(#arr-recycle)"/><text x="318.0" y="140" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text><line x1="514" y1="233.0" x2="554" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="534.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="554,205 618,205 611,261 561,261" class="clarifier"/><text x="586.0" y="238.0" class="zonelabel-d">SC</text><line x1="618" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="939.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 568 261 L 568 320 L 186 320 L 186 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="377.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="604" y1="261" x2="604" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="612" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="402.0,470 458.0,470 452.0,526 408.0,526" class="thickener"/><text x="430.0" y="503.0" class="zonelabel-d">MT</text><rect x="568.0" y="470" width="64" height="56" class="digester-ae" rx="3"/><text x="600.0" y="503.0" class="zonelabel">AED</text><circle cx="588.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="596.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="604.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="612.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><rect x="742.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="770.0" y="503.0" class="zonelabel-d">DW</text><line x1="458.0" y1="498.0" x2="568.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="632.0" y1="498.0" x2="742.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="604" y1="420" x2="430.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="430.0" y1="420" x2="430.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="798.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="975.0" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="430.0" y1="526" x2="430.0" y2="600" class="recycle"/><line x1="770.0" y1="526" x2="770.0" y2="600" class="recycle"/><line x1="430.0" y1="600" x2="770.0" y2="600" class="recycle"/><line x1="430.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="176" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="172" y1="255.0" x2="172" y2="236.0" class="recycle"/><line x1="172" y1="236.0" x2="176" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>

.. raw:: html


   <div class="wrrf-diagram" id="I3">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">I3</text><text x="36" y="78" class="subtitle">BNR · No primary clarifier · 5-stage Bardenpho · No digestion (thickening + dewatering only)</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="178" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><rect x="178" y="205" width="56" height="56" class="anaerobic" rx="2"/><text x="206.0" y="237.0" class="zonelabel">Anae</text><rect x="234" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="262.0" y="237.0" class="zonelabel">Anox</text><rect x="290" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="318.0" y="237.0" class="zonelabel">Aer</text><circle cx="310.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="318.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="326.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="346" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="374.0" y="237.0" class="zonelabel">Aer</text><circle cx="366.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="374.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="382.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="402" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="430.0" y="237.0" class="zonelabel">Anox</text><rect x="458" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="486.0" y="237.0" class="zonelabel">Aer</text><circle cx="478.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="486.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="494.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><path d="M 388.0 205 L 388.0 145 L 248.0 145 L 248.0 205" class="recycle" marker-end="url(#arr-recycle)"/><text x="318.0" y="140" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text><line x1="514" y1="233.0" x2="554" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="534.0" y="225.0" class="streamlabel" text-anchor="middle">ML</text><polygon points="554,205 618,205 611,261 561,261" class="clarifier"/><text x="586.0" y="238.0" class="zonelabel-d">SC</text><line x1="618" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="939.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 568 261 L 568 320 L 186 320 L 186 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="377.0" y="315" class="streamlabel" text-anchor="middle">RAS</text><line x1="604" y1="261" x2="604" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="612" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="489.0,470 545.0,470 539.0,526 495.0,526" class="thickener"/><text x="517.0" y="503.0" class="zonelabel-d">MT</text><rect x="655.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="683.0" y="503.0" class="zonelabel-d">DW</text><line x1="545.0" y1="498.0" x2="655.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="604" y1="420" x2="517.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="517.0" y1="420" x2="517.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="711.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="931.5" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="517.0" y1="526" x2="517.0" y2="600" class="recycle"/><line x1="683.0" y1="526" x2="683.0" y2="600" class="recycle"/><line x1="517.0" y1="600" x2="683.0" y2="600" class="recycle"/><line x1="517.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="176" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="172" y1="255.0" x2="172" y2="236.0" class="recycle"/><line x1="172" y1="236.0" x2="176" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>


N Series: BNR with MBR
----------------------
Five-stage Bardenpho liquid train with a completely-mixed membrane bioreactor (MBR) replacing the conventional secondary clarifier. N1 includes primary clarification and anaerobic digestion; N2 omits primary clarification and uses aerobic digestion.

.. raw:: html


   <div class="wrrf-diagram" id="N1">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">N1</text><text x="36" y="78" class="subtitle">BNR · Primary clarifier · 5-stage Bardenpho with membrane bioreactor · Anaerobic digestion</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="118" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><polygon points="118,205 182,205 175,261 125,261" class="clarifier"/><text x="150.0" y="238.0" class="zonelabel-d">PC</text><line x1="182" y1="233.0" x2="210" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="196.0" y="225.0" class="streamlabel" text-anchor="middle">PE</text><rect x="210" y="205" width="56" height="56" class="anaerobic" rx="2"/><text x="238.0" y="237.0" class="zonelabel">Anae</text><rect x="266" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="294.0" y="237.0" class="zonelabel">Anox</text><rect x="322" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="350.0" y="237.0" class="zonelabel">Aer</text><circle cx="342.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="350.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="358.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="378" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="406.0" y="237.0" class="zonelabel">Aer</text><circle cx="398.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="406.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="414.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="434" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="462.0" y="237.0" class="zonelabel">Anox</text><rect x="490" y="205" width="66" height="56" class="mbr" rx="2"/><line x1="499" y1="210" x2="499" y2="256" stroke="#fff" stroke-width="0.9" opacity="0.7"/><line x1="507" y1="210" x2="507" y2="256" stroke="#fff" stroke-width="0.9" opacity="0.7"/><line x1="515" y1="210" x2="515" y2="256" stroke="#fff" stroke-width="0.9" opacity="0.7"/><line x1="523" y1="210" x2="523" y2="256" stroke="#fff" stroke-width="0.9" opacity="0.7"/><line x1="531" y1="210" x2="531" y2="256" stroke="#fff" stroke-width="0.9" opacity="0.7"/><line x1="539" y1="210" x2="539" y2="256" stroke="#fff" stroke-width="0.9" opacity="0.7"/><text x="523.0" y="238.0" class="zonelabel">MBR</text><rect x="178.0" y="113" width="120" height="30" class="chem" rx="15"/><line x1="238.0" y1="143" x2="238.0" y2="204" class="dose" marker-end="url(#arr-dose)"/><text x="238" y="128" class="zonelabel" dominant-baseline="central">External C</text><path d="M 420.0 205 L 420.0 145 L 280.0 145 L 280.0 205" class="recycle" marker-end="url(#arr-recycle)"/><text x="350.0" y="170" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text><line x1="556" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="918.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 550 261 L 550 320 L 350.0 320 L 350.0 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="450.0" y="315" class="streamlabel" text-anchor="middle">MBR mixed-liquor recycle</text><line x1="150.0" y1="261" x2="150.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="142.0" y="275" class="streamlabel" text-anchor="end">PS</text><line x1="522" y1="261" x2="522" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="530" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="402.0,470 458.0,470 452.0,526 408.0,526" class="thickener"/><text x="430.0" y="503.0" class="zonelabel-d">MT</text><path d="M568.0,480 Q568.0,470 600.0,470 Q632.0,470 632.0,480 L632.0,526 L568.0,526 z" class="digester-an"/><text x="600.0" y="503.0" class="zonelabel">AD</text><line x1="600.0" y1="468" x2="600.0" y2="448" class="gas" marker-end="url(#arr-gas)"/><text x="605.0" y="445" class="streamlabel" fill="#6b4a8a">biogas</text><rect x="742.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="770.0" y="503.0" class="zonelabel-d">DW</text><line x1="458.0" y1="498.0" x2="568.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="632.0" y1="498.0" x2="742.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="150.0" y1="420" x2="430.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="430.0" y1="420" x2="430.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="522" y1="420" x2="430.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="798.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="975.0" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="430.0" y1="526" x2="430.0" y2="600" class="recycle"/><line x1="770.0" y1="526" x2="770.0" y2="600" class="recycle"/><line x1="430.0" y1="600" x2="770.0" y2="600" class="recycle"/><line x1="430.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="116" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="112" y1="255.0" x2="112" y2="236.0" class="recycle"/><line x1="112" y1="236.0" x2="116" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>

.. raw:: html


   <div class="wrrf-diagram" id="N2">
   <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" width="1200" height="800">

   <defs>
     <marker id="arr-main" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#1a1a1a"/></marker>
     <marker id="arr-recycle" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#c2272d"/></marker>
     <marker id="arr-sludge" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#7a5a30"/></marker>
     <marker id="arr-gas" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#6b4a8a"/></marker>
     <marker id="arr-dose" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
       <path d="M0,0 L10,5 L0,10 z" fill="#d97a4a"/></marker>
   </defs>
   <rect class="bg" width="1200" height="800"/><text x="36" y="38" class="title">N2</text><text x="36" y="78" class="subtitle">BNR · No primary clarifier · 5-stage Bardenpho with membrane bioreactor · Aerobic digestion</text><text x="600.0" y="190" class="section">LIQUID TREATMENT</text><line x1="60" y1="233.0" x2="178" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="60" y="224.0" class="streamlabel" text-anchor="start">RWW</text><rect x="178" y="205" width="56" height="56" class="anaerobic" rx="2"/><text x="206.0" y="237.0" class="zonelabel">Anae</text><rect x="234" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="262.0" y="237.0" class="zonelabel">Anox</text><rect x="290" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="318.0" y="237.0" class="zonelabel">Aer</text><circle cx="310.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="318.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="326.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="346" y="205" width="56" height="56" class="aerobic" rx="2"/><text x="374.0" y="237.0" class="zonelabel">Aer</text><circle cx="366.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="374.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><circle cx="382.0" cy="253" r="1.5" fill="#fff" opacity="0.85"/><rect x="402" y="205" width="56" height="56" class="anoxic" rx="2"/><text x="430.0" y="237.0" class="zonelabel">Anox</text><rect x="458" y="205" width="66" height="56" class="mbr" rx="2"/><line x1="467" y1="210" x2="467" y2="256" stroke="#fff" stroke-width="0.9" opacity="0.7"/><line x1="475" y1="210" x2="475" y2="256" stroke="#fff" stroke-width="0.9" opacity="0.7"/><line x1="483" y1="210" x2="483" y2="256" stroke="#fff" stroke-width="0.9" opacity="0.7"/><line x1="491" y1="210" x2="491" y2="256" stroke="#fff" stroke-width="0.9" opacity="0.7"/><line x1="499" y1="210" x2="499" y2="256" stroke="#fff" stroke-width="0.9" opacity="0.7"/><line x1="507" y1="210" x2="507" y2="256" stroke="#fff" stroke-width="0.9" opacity="0.7"/><text x="491.0" y="238.0" class="zonelabel">MBR</text><rect x="146.0" y="113" width="120" height="30" class="chem" rx="15"/><line x1="206.0" y1="143" x2="206.0" y2="204" class="dose" marker-end="url(#arr-dose)"/><text x="206" y="128" class="zonelabel" dominant-baseline="central">External C</text><path d="M 388.0 205 L 388.0 145 L 248.0 145 L 248.0 205" class="recycle" marker-end="url(#arr-recycle)"/><text x="318.0" y="170" class="streamlabel" text-anchor="middle">MLR (internal recycle)</text><line x1="524" y1="233.0" x2="1164" y2="233.0" class="mainflow" marker-end="url(#arr-main)"/><text x="902.0" y="225.0" class="streamlabel" text-anchor="middle">Effluent</text><path d="M 518 261 L 518 320 L 318.0 320 L 318.0 261" class="recycle" marker-end="url(#arr-recycle)"/><text x="418.0" y="315" class="streamlabel" text-anchor="middle">MBR mixed-liquor recycle</text><line x1="490" y1="261" x2="490" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><text x="498" y="275" class="streamlabel" text-anchor="start">WAS</text><text x="600.0" y="392" class="section">SOLIDS TREATMENT</text><polygon points="402.0,470 458.0,470 452.0,526 408.0,526" class="thickener"/><text x="430.0" y="503.0" class="zonelabel-d">MT</text><rect x="568.0" y="470" width="64" height="56" class="digester-ae" rx="3"/><text x="600.0" y="503.0" class="zonelabel">AED</text><circle cx="588.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="596.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="604.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><circle cx="612.0" cy="517" r="1.5" fill="#fff" opacity="0.85"/><rect x="742.0" y="470" width="56" height="56" class="dw" rx="2"/><text x="770.0" y="503.0" class="zonelabel-d">DW</text><line x1="458.0" y1="498.0" x2="568.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="632.0" y1="498.0" x2="742.0" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><line x1="490" y1="420" x2="430.0" y2="420" class="sludge" marker-end="url(#arr-sludge)"/><line x1="430.0" y1="420" x2="430.0" y2="470" class="sludge" marker-end="url(#arr-sludge)"/><line x1="798.0" y1="498.0" x2="1152" y2="498.0" class="sludge" marker-end="url(#arr-sludge)"/><text x="975.0" y="491.0" class="streamlabel" text-anchor="middle">sludge cake</text><line x1="430.0" y1="526" x2="430.0" y2="600" class="recycle"/><line x1="770.0" y1="526" x2="770.0" y2="600" class="recycle"/><line x1="430.0" y1="600" x2="770.0" y2="600" class="recycle"/><line x1="430.0" y1="600" x2="50" y2="600" class="recycle"/><text x="600.0" y="594" class="streamlabel" text-anchor="middle">reject water (centrate + thickener supernatant)</text><line x1="50" y1="600" x2="50" y2="255.0" class="recycle"/><line x1="50" y1="255.0" x2="176" y2="255.0" class="recycle" marker-end="url(#arr-recycle)"/><line x1="172" y1="255.0" x2="172" y2="236.0" class="recycle"/><line x1="172" y1="236.0" x2="176" y2="236.0" class="recycle" marker-end="url(#arr-recycle)"/><g transform="translate(36,630)">
   <rect x="0" y="0" width="14" height="10" class="anaerobic" rx="1.5"/><text x="19" y="9" class="streamlabel">Anaerobic</text>
   <rect x="149" y="0" width="14" height="10" class="anoxic" rx="1.5"/><text x="168" y="9" class="streamlabel">Anoxic</text>
   <rect x="260" y="0" width="14" height="10" class="aerobic" rx="1.5"/><text x="279" y="9" class="streamlabel">Aerobic</text>
   <rect x="371" y="0" width="14" height="10" class="clarifier" rx="1.5"/><text x="390" y="9" class="streamlabel">Clarifier</text>
   <rect x="515" y="0" width="14" height="10" class="digester-an" rx="1.5"/><text x="534" y="9" class="streamlabel">Anaerobic dig. (AD)</text>
   <rect x="759" y="0" width="14" height="10" class="digester-ae" rx="1.5"/><text x="778" y="9" class="streamlabel">Aerobic dig. (AED)</text>
   <rect x="0" y="30" width="14" height="10" class="thickener" rx="1.5"/><text x="19" y="39" class="streamlabel">Thickener</text>
   <rect x="149" y="30" width="14" height="10" class="dw" rx="1.5"/><text x="168" y="39" class="streamlabel">Dewatering</text>
   <rect x="306" y="30" width="14" height="10" class="mbr" rx="1.5"/><text x="325" y="39" class="streamlabel">Membrane (MBR)</text>
   <rect x="500" y="30" width="14" height="10" class="chem" rx="1.5"/><text x="519" y="39" class="streamlabel">Chemical dose</text>
   <line x1="0" y1="70" x2="24" y2="70" class="mainflow" marker-end="url(#arr-main)"/><text x="30" y="74" class="streamlabel">Main flow</text>
   <line x1="145" y1="70" x2="169" y2="70" class="recycle" marker-end="url(#arr-recycle)"/><text x="175" y="74" class="streamlabel">Recycle</text>
   <line x1="273" y1="70" x2="297" y2="70" class="sludge" marker-end="url(#arr-sludge)"/><text x="303" y="74" class="streamlabel">Sludge</text>
   <line x1="395" y1="70" x2="419" y2="70" class="gas" marker-end="url(#arr-gas)"/><text x="425" y="74" class="streamlabel">Biogas</text>
   <line x1="517" y1="70" x2="541" y2="70" class="dose" marker-end="url(#arr-dose)"/><text x="547" y="74" class="streamlabel">Chemical addition</text>
   </g></svg>
   </div>

.. References
.. [1] Zhang, X.; Rai, S.; Wang, Z.; Li, Y.; Guest, J. S. An Agile Benchmarking Framework for Wastewater Resource Recovery Technologies. npj Clean Water 2025, 9 (1), 4. https://doi.org/10.1038/s41545-025-00537-4.

.. [2] Tarallo, S., Shaw, A., Kohl, P. & Eschborn, R. A Guide to Net-Zero Energy Solutions for Water Resource Recovery Facilities. https://iwaponline.com/ebooks/book/293/ (2015).