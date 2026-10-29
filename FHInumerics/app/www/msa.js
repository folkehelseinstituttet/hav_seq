// Alignment viewer (Alignment page, decision 0011). Draws the cell classes computed in
// R/msa.R; it compares no bases itself. Classes:
//   0 same as the anchor  1 SNP  2 ambiguity code that fits  3 ambiguity code that cannot
//   match  4 indel  5 N (no data)  6 missing start/end (no data)  7 the anchor has no base
//   8 anchor base (anchor row only)
// In the window, R sends 'msa_data' (names, letters, metadata) when the alignment changes
// and 'msa_view' (anchor, classes, row statistics, column masks) when the anchor changes;
// the page sends back input$msa_anchor_pick (right-click "Use as anchor") and
// input$msa_window (what is on screen, for the PDF and FASTA downloads). msa_html() in
// R/msa.R puts this file in a standalone page and calls fhiMsa.load() instead.
(function() {
  var FILL = ['', '#e6550d', '#c6dbef', '#2166ac', '#7b3294', '#e6e9ed', '#e6e9ed', '#cfd4da', '#fff3cd'];
  var INK = ['#9aa4b0', '#ffffff', '#08306b', '#ffffff', '#ffffff', '#8c959f', '#000000', '#333333', '#1f2328'];
  var RANK = [7, 0, 3, 1, 2, 5, 6, 4, 8];   // lower = drawn first when a pixel covers several columns
  var WHAT = ['same as the anchor', 'SNP (counted)', 'ambiguity code that fits (not counted)',
              'ambiguity code that cannot match (counted as a SNP)', 'indel (not counted as a SNP)',
              'N: no data', 'missing start or end: no data', 'nothing to compare: the anchor has no base here'];
  var LABEL_W = 250, ROW_H = 14, AXIS_H = 22, ANCHOR_H = 18, OV_H = 46, MAX_PX = 16, LETTER_PX = 8,
      MAX_FULL_ROWS = 300, ROW_BG = '#f3f5f8', HOT = '#cf222e', LINE = '#9aa4b0', ACCENT = '#0969da',
      MUTED = '#656d76';

  var S = {data: null, view: null, opts: {}, shown: [], order: [], px: 1, start: 0, top: 0, selected: -1,
           varOnly: false, hideGap: false, mark: -1, inserts: [], sort: null, fitted: true};
  var el = {};
  var arr = function(x) { return x == null ? [] : (Array.isArray(x) ? x : [x]); };
  var plural = function(n, word) { return n + ' ' + word + (n === 1 ? '' : 's'); };

  function init() {
    var root = document.getElementById('msa');
    if (!root || el.root) return;
    el.root = root;
    var standalone = !!S.opts.standalone;
    root.innerHTML =
      '<div class="msa-bar">' +
      '<label>Sort rows <select class="msa-sort"></select></label>' +
      (standalone ? '<label><input type="checkbox" data-opt="variable"> Variable columns only</label>' +
                    '<label><input type="checkbox" data-opt="hide_gaps"> Hide columns that are a gap in the anchor</label>' : '') +
      '<label>Go to position <input type="number" min="1" class="msa-goto"></label>' +
      '<button type="button" data-act="in">Zoom in</button><button type="button" data-act="out">Zoom out</button>' +
      '<button type="button" data-act="fit">Whole alignment</button><button type="button" data-act="letters">Letters</button>' +
      '<button type="button" data-act="prev">&#9664; Previous difference</button>' +
      '<button type="button" data-act="next">Next difference &#9654;</button>' +
      '<span class="msa-range"></span></div>' +
      '<div class="msa-hint">Mouse wheel = zoom &middot; drag = move &middot; double-click = whole alignment &middot; ' +
      'point = details &middot; click a name = choose that row' +
      (standalone ? '' : ' &middot; right-click a name = use it as anchor') + '</div>' +
      '<div class="msa-message"></div>' +
      '<canvas class="msa-overview"></canvas><div class="msa-legend"></div>' +
      '<div class="msa-frame"><canvas class="msa-canvas"></canvas><div class="msa-vscroll"><div></div></div></div>' +
      '<div class="msa-tip fhi-tip"></div><div class="msa-menu"></div>';
    el.canvas = root.querySelector('.msa-canvas');
    el.ov = root.querySelector('.msa-overview');
    el.vscroll = root.querySelector('.msa-vscroll');
    el.tip = root.querySelector('.msa-tip');
    el.menu = root.querySelector('.msa-menu');
    el.msg = root.querySelector('.msa-message');
    el.range = root.querySelector('.msa-range');
    el.legend = root.querySelector('.msa-legend');
    el.sort = root.querySelector('.msa-sort');
    el.gotoInput = root.querySelector('.msa-goto');
    root.querySelector('.msa-bar').addEventListener('click', function(e) {
      var b = e.target.closest('button');
      if (b) act(b.getAttribute('data-act'));
    });
    root.querySelectorAll('input[data-opt]').forEach(function(box) {
      box.addEventListener('change', function() {
        setOptions({variable: root.querySelector('input[data-opt=variable]').checked,
                    hide_gaps: root.querySelector('input[data-opt=hide_gaps]').checked});
      });
    });
    el.sort.addEventListener('change', function() { S.sort = el.sort.value; computeOrder(); layout(); draw(); });
    el.gotoInput.addEventListener('keydown', function(e) { if (e.key === 'Enter') act('goto'); });
    el.vscroll.addEventListener('scroll', function() { S.top = el.vscroll.scrollTop / ROW_H; draw(); });
    el.canvas.addEventListener('wheel', wheel, {passive: false});
    el.canvas.addEventListener('mousemove', hover);
    el.canvas.addEventListener('mouseleave', function() { el.tip.style.display = 'none'; });
    el.canvas.addEventListener('pointerdown', pointerdown);
    el.canvas.addEventListener('dblclick', function() { if (S.view) { fit(); draw(); } });
    el.canvas.addEventListener('contextmenu', contextmenu);
    el.ov.addEventListener('pointerdown', function(e) { el.ov.setPointerCapture(e.pointerId); S.ovDrag = true; overviewTo(e); });
    el.ov.addEventListener('pointermove', function(e) { if (S.ovDrag) overviewTo(e); });
    el.ov.addEventListener('pointerup', function() { S.ovDrag = false; });
    window.addEventListener('pointermove', dragmove);
    window.addEventListener('pointerup', dragend);
    document.addEventListener('click', function(e) { if (!el.menu.contains(e.target)) el.menu.style.display = 'none'; });
    new ResizeObserver(function() {
      layout();
      if (S.fitted && S.view) fit();   // keep showing the whole alignment when the width changes
      if (S.pending && S.width) { var v = S.pending; S.pending = null; setView(v); } else draw();
    }).observe(root);
  }

  // ---- Data -------------------------------------------------------------------------
  function setData(d) {
    init();
    S.data = d;
    S.view = null;
    S.selected = -1;
    S.mark = -1;
    if (d) {
      ['names', 'chars', 'cols', 'query', 'dropped', 'revcomp', 'meta'].forEach(function(k) { d[k] = arr(d[k]); });
      if (d.dates) { d.dates.text = arr(d.dates.text); d.dates.sort = arr(d.dates.sort); }
      var opts = [['tree', 'Tree order (as the SNP table)'], ['fewest', 'Fewest SNPs first'], ['most', 'Most SNPs first'],
                  ['name', 'Name'], ['date', 'Sample date']].filter(function(o) {
        return (o[0] !== 'tree' || d.tree_order) && (o[0] !== 'date' || d.dates);
      });
      el.sort.innerHTML = '';
      opts.forEach(function(o) { var x = document.createElement('option'); x.value = o[0]; x.textContent = o[1]; el.sort.appendChild(x); });
      var valid = function(s) { return opts.some(function(o) { return o[0] === s; }); };
      if (!valid(S.sort)) S.sort = valid(S.opts.sort) ? S.opts.sort : (d.tree_order ? 'tree' : 'fewest');
      el.sort.value = S.sort;
    }
    el.msg.textContent = '';
  }

  function setView(v) {
    init();
    if (!S.data) return;
    ['rows', 'classes', 'diffs', 'ndiff'].forEach(function(k) { v[k] = arr(v[k]); });
    ['indel_events', 'indel_cols', 'ambiguity', 'nodata'].forEach(function(k) { v.stats[k] = arr(v.stats[k]); });
    var keep = S.view && S.shown.length ? S.shown[Math.floor(S.start)] : null;
    S.view = v;
    buildShown();
    computeOrder();
    layout();
    if (!S.width) { S.pending = v; return; }
    if (v.select != null) S.selected = S.data.names.indexOf(v.select);
    if (v.goto != null) goColumn(v.goto - 1, true);
    else if (keep != null && !S.fitted) S.start = nearestShown(keep);
    else fit();
    drawLegend();
    draw();
  }

  function setOptions(o) {
    S.varOnly = !!o.variable;
    S.hideGap = !!o.hide_gaps;
    if (!S.view) return;
    var keep = S.shown[Math.floor(S.start)];
    buildShown();
    if (S.fitted) fit(); else S.start = nearestShown(keep || 0);
    if (S.px < minPx()) S.px = minPx();
    draw();
  }

  // Row order below the anchor: R's order (tree, or fewest SNPs), or the chosen sort.
  function computeOrder() {
    if (!S.view) return;
    var v = S.view, d = S.data, rows = v.rows.slice(), pos = {};
    rows.forEach(function(r, i) { pos[r] = i; });
    var snps = function(r) { return v.diffs[r] == null ? Infinity : v.diffs[r]; };
    var by = {
      tree: null,
      fewest: function(a, b) { return snps(a) - snps(b); },
      most: function(a, b) { return (v.diffs[b] == null ? -1 : v.diffs[b]) - (v.diffs[a] == null ? -1 : v.diffs[a]); },
      name: function(a, b) { return d.names[a].localeCompare(d.names[b], undefined, {numeric: true}); },
      date: function(a, b) {
        var x = d.dates.sort[a] || '9999', y = d.dates.sort[b] || '9999';
        return x < y ? -1 : x > y ? 1 : 0;
      }
    }[S.sort];
    if (by) rows.sort(function(a, b) { return by(a, b) || pos[a] - pos[b]; });
    S.order = rows;
    S.top = 0;
  }

  // Columns shown (indices into the alignment), after the two toggles. Also notes where
  // anchor-gap columns were hidden, for the marks on the anchor row.
  function buildShown() {
    var v = S.view, n = S.data.ncol, hidden = 0;
    S.shown = [];
    S.inserts = [];
    for (var c = 0; c < n; c++) {
      var gap = v.anchor_gap.charCodeAt(c) === 49;
      if ((S.varOnly && v.variable.charCodeAt(c) !== 49) || (S.hideGap && gap)) {
        if (S.hideGap && gap) hidden++;
        continue;
      }
      if (hidden) { S.inserts.push({at: S.shown.length, n: hidden}); hidden = 0; }
      S.shown.push(c);
    }
    if (hidden) S.inserts.push({at: S.shown.length, n: hidden});
  }
  function nearestShown(c) {
    for (var k = 0; k < S.shown.length; k++) if (S.shown[k] >= c) return k;
    return Math.max(0, S.shown.length - 1);
  }

  // ---- Geometry ---------------------------------------------------------------------
  function mainW() { return Math.max(10, S.width - LABEL_W - 4); }
  function minPx() { return Math.min(MAX_PX, mainW() / Math.max(1, S.shown.length)); }
  function nRows() { return S.view ? S.order.length : 0; }
  function rowsTop() { return AXIS_H + ANCHOR_H; }
  function scrolling() { return nRows() > MAX_FULL_ROWS; }
  function visibleRows() { return Math.max(1, Math.floor((S.height - rowsTop() - 4) / ROW_H)); }

  // All rows are drawn and the page scrolls; above MAX_FULL_ROWS the rows get their own
  // scroll bar (one canvas of thousands of rows would be too large).
  function layout() {
    if (!el.root) return;
    var w = el.root.clientWidth - (scrolling() ? 18 : 0);
    if (w <= 0) return;
    S.width = w;
    var rows = nRows();
    S.height = scrolling() ? rowsTop() + 4 + ROW_H * Math.max(10, Math.floor((window.innerHeight - 260) / ROW_H))
                           : rowsTop() + 4 + rows * ROW_H;
    var dpr = window.devicePixelRatio || 1;
    [[el.canvas, S.height], [el.ov, OV_H]].forEach(function(p) {
      p[0].width = Math.round(w * dpr); p[0].height = Math.round(p[1] * dpr);
      p[0].style.width = w + 'px'; p[0].style.height = p[1] + 'px';
    });
    el.vscroll.style.display = scrolling() ? 'block' : 'none';
    el.vscroll.style.height = S.height + 'px';
    el.vscroll.firstChild.style.height = (S.height + Math.max(0, rows - visibleRows()) * ROW_H) + 'px';
    if (S.px < minPx()) S.px = minPx();
  }

  function clamp() {
    S.start = Math.min(Math.max(0, S.start), Math.max(0, S.shown.length - mainW() / S.px));
    S.top = scrolling() ? Math.min(Math.max(0, S.top), Math.max(0, nRows() - visibleRows())) : 0;
  }
  function zoom(factor, centreK) {
    var mw = mainW(), centre = centreK != null ? centreK : S.start + mw / S.px / 2, rel = (centre - S.start) * S.px / mw;
    S.px = Math.min(MAX_PX, Math.max(minPx(), S.px * factor));
    S.fitted = S.px <= minPx() + 1e-9;
    S.start = centre - rel * mw / S.px;
    draw();
  }
  function fit() { S.px = minPx(); S.start = 0; S.fitted = true; }

  // Show column c (index into the alignment) in the middle, at letter size.
  function goColumn(c, mark) {
    var k = nearestShown(c);
    if (S.px < LETTER_PX * 1.5) S.px = Math.min(MAX_PX, 12);
    S.fitted = false;
    S.start = k - mainW() / S.px / 2 + 0.5;
    if (mark) S.mark = S.shown[k];
    var i = S.order.indexOf(S.selected);
    if (scrolling() && i >= 0 && (i < S.top || i >= S.top + visibleRows())) S.top = i - Math.floor(visibleRows() / 2);
  }

  // ---- Actions ----------------------------------------------------------------------
  function act(what) {
    if (!S.view) return;
    el.msg.textContent = '';
    if (what === 'in') zoom(2);
    if (what === 'out') zoom(0.5);
    if (what === 'fit') { fit(); draw(); }
    if (what === 'letters') zoom(12 / S.px);
    if (what === 'goto') gotoPosition(parseInt(el.gotoInput.value, 10));
    if (what === 'next' || what === 'prev') nextDifference(what === 'next' ? 1 : -1);
  }

  function gotoPosition(p) {
    var a = S.view.anchor_chars, count = 0;
    if (!(p > 0)) { el.msg.textContent = 'Type a position (a base number in ' + S.view.anchor_name + ').'; return; }
    for (var c = 0; c < a.length; c++) {
      if (a[c] !== '-') count++;
      if (a[c] !== '-' && count === p) { goColumn(c, true); draw(); return; }
    }
    el.msg.textContent = S.view.anchor_name + ' has only ' + count + ' bases.';
  }

  // The next column (from the marked one, or the middle of the window) where the chosen
  // row - or any row - differs from the anchor (SNP or indel).
  function nextDifference(dir) {
    var v = S.view, from = S.mark >= 0 ? S.mark : (S.shown[Math.floor(S.start + mainW() / S.px / 2)] || 0);
    var shown = {};
    S.shown.forEach(function(c) { shown[c] = true; });
    for (var c = from + dir; c >= 0 && c < S.data.ncol; c += dir) {
      var hit = S.selected >= 0 ? [1, 3, 4].indexOf(v.classes[S.selected].charCodeAt(c) - 48) >= 0 : v.ndiff[c] > 0;
      if (hit && shown[c]) { goColumn(c, true); draw(); return; }
    }
    el.msg.textContent = 'No more differences ' + (dir > 0 ? 'after' : 'before') + ' this column' +
      (S.selected >= 0 ? ' in ' + S.data.names[S.selected] : '') + '.';
  }

  // ---- Drawing ----------------------------------------------------------------------
  function context(cv) {
    var ctx = cv.getContext('2d'), dpr = window.devicePixelRatio || 1;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.textBaseline = 'middle';
    return ctx;
  }
  function X(k) { return LABEL_W + (k - S.start) * S.px; }

  function draw() {
    if (!el.root || !S.view || !S.width) return;
    clamp();
    if (scrolling()) el.vscroll.scrollTop = S.top * ROW_H;
    drawOverview();
    var ctx = context(el.canvas);
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, S.width, S.height);
    drawAxis(ctx);
    var first = Math.floor(S.top), n = scrolling() ? Math.min(visibleRows(), nRows() - first) : nRows();
    drawRow(ctx, -1, AXIS_H, -1);
    for (var i = 0; i < n; i++) drawRow(ctx, S.order[first + i], rowsTop() + i * ROW_H, i);
    ctx.save(); ctx.beginPath(); ctx.rect(LABEL_W, AXIS_H, mainW(), S.height - AXIS_H); ctx.clip();
    drawCompared(ctx, AXIS_H, ANCHOR_H + n * ROW_H);
    drawInserts(ctx);
    if (S.mark >= 0) {
      var k = S.shown.indexOf(S.mark);
      if (k >= 0) {
        ctx.strokeStyle = '#ff7f00'; ctx.lineWidth = 2;
        ctx.strokeRect(X(k), AXIS_H, Math.max(S.px, 2), ANCHOR_H + n * ROW_H + 2);
      }
    }
    ctx.restore();
    var kFirst = Math.floor(S.start), kLast = Math.min(S.shown.length, Math.ceil(S.start + mainW() / S.px));
    el.range.textContent = S.shown.length ? 'columns ' + S.data.cols[S.shown[kFirst]] + '–' +
      S.data.cols[S.shown[kLast - 1]] + ' of ' + S.data.cols[S.data.ncol - 1] +
      (S.shown.length < S.data.ncol ? ' (' + S.shown.length + ' shown)' : '') : 'no columns to show';
    reportWindow(first, n);
  }

  // The class string, letters and name of a row (-1 = the anchor row).
  function rowInfo(r) {
    if (r < 0) return {cls: S.view.anchor_classes, chars: S.view.anchor_chars, name: S.view.anchor_name};
    return {cls: S.view.classes[r], chars: S.data.chars[r], name: S.data.names[r]};
  }
  function isQuery(r) { return r >= 0 ? S.data.query[r] : S.view.anchor != null && S.data.query[S.view.anchor]; }

  function drawRow(ctx, r, y, i) {
    var info = rowInfo(r), mw = mainW(), px = S.px, h = r < 0 ? ANCHOR_H - 2 : ROW_H;
    // Background: anchor pale yellow, chosen row light blue, every other row grey.
    var bg = r < 0 ? '#fff8e1' : (r === S.selected ? '#dbe9fd' : (i % 2 === 0 ? ROW_BG : ''));
    if (bg) { ctx.fillStyle = bg; ctx.fillRect(0, y, S.width, h); }
    if (isQuery(r)) { ctx.fillStyle = '#fde0dd'; ctx.fillRect(0, y, LABEL_W - 6, h); }
    // Label: name, date, SNPs to the anchor.
    var dropped = r >= 0 && S.data.dropped[r], dates = S.data.dates;
    ctx.font = (dropped ? 'italic ' : '') + (isQuery(r) || r < 0 ? 'bold ' : '') + '11px system-ui, sans-serif';
    ctx.fillStyle = dropped ? '#8c959f' : (isQuery(r) ? '#d7301f' : '#1f2328');
    ctx.textAlign = 'left';
    var nameW = LABEL_W - (dates && r >= 0 ? 112 : 50);
    ctx.save(); ctx.beginPath(); ctx.rect(0, y, nameW, h); ctx.clip();
    ctx.fillText((r < 0 ? '▶ ' : '') + info.name, 4, y + h / 2 + 1);
    ctx.restore();
    ctx.font = '11px system-ui, sans-serif';
    ctx.fillStyle = MUTED;
    if (dates && r >= 0) ctx.fillText(dates.text[r] || '', nameW + 4, y + h / 2 + 1);
    ctx.textAlign = 'right';
    ctx.fillText(r < 0 ? 'anchor' : (S.view.diffs[r] == null ? '–' : String(S.view.diffs[r])), LABEL_W - 10, y + h / 2 + 1);

    // Cells: only what differs is drawn; no data as a light grey band.
    var cls = info.cls, chars = info.chars;
    ctx.save(); ctx.beginPath(); ctx.rect(LABEL_W, y, mw, h); ctx.clip();
    if (px >= 1) {
      var kFirst = Math.floor(S.start), kLast = Math.min(S.shown.length - 1, Math.ceil(S.start + mw / px));
      var letters = px >= LETTER_PX;
      ctx.font = Math.min(12, Math.floor(px)) + 'px ui-monospace, Consolas, monospace';
      ctx.textAlign = 'center';
      for (var k = kFirst; k <= kLast; k++) {
        var c = S.shown[k], kc = cls.charCodeAt(c) - 48, x = X(k);
        if (FILL[kc]) { ctx.fillStyle = FILL[kc]; ctx.fillRect(x, y + 1, Math.max(1, px - (px > 3 ? 1 : 0)), h - 2); }
        if (letters && kc !== 6) {
          ctx.fillStyle = INK[kc];
          ctx.fillText(kc === 0 ? '·' : chars[c], x + px / 2, y + h / 2 + 1);
        }
      }
    } else {
      // Each pixel covers several columns: the most important class among them.
      for (var p = 0; p < mw; p++) {
        var a = Math.floor(S.start + p / px), b = Math.min(S.shown.length, Math.floor(S.start + (p + 1) / px));
        if (a >= S.shown.length) break;
        var best = 0, bestRank = 99;
        for (var j = a; j < Math.max(b, a + 1); j++) {
          var kk = cls.charCodeAt(S.shown[j]) - 48;
          if (RANK[kk] < bestRank) { bestRank = RANK[kk]; best = kk; }
        }
        if (FILL[best]) { ctx.fillStyle = FILL[best]; ctx.fillRect(LABEL_W + p, y + 1, 1, h - 2); }
      }
    }
    ctx.restore();
    if (r < 0) { ctx.fillStyle = '#d0d7de'; ctx.fillRect(0, y + h + 1, S.width, 1); }
  }

  // Column numbers of the alignment. Where columns are skipped (Variable columns only /
  // hidden anchor gaps), a red tick marks the jump.
  function drawAxis(ctx) {
    var mw = mainW(), px = S.px, cols = S.data.cols;
    ctx.font = '11px system-ui, sans-serif';
    ctx.fillStyle = MUTED; ctx.textAlign = 'left';
    ctx.fillText('alignment column', 4, AXIS_H / 2 - 2);
    var steps = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000], step = steps[steps.length - 1];
    for (var i = 0; i < steps.length; i++) if (steps[i] * px >= 80) { step = steps[i]; break; }
    var kFirst = Math.floor(S.start), kLast = Math.min(S.shown.length - 1, Math.ceil(S.start + mw / px));
    var skipped = S.varOnly || S.hideGap, lastLabel = -1e9;
    ctx.textAlign = 'center';
    for (var k = Math.max(0, kFirst); k <= kLast; k++) {
      var x = X(k + 0.5), col = cols[S.shown[k]];
      if (x < LABEL_W) continue;
      if (k > 0 && S.shown[k] !== S.shown[k - 1] + 1) { ctx.fillStyle = HOT; ctx.fillRect(X(k) - 1, AXIS_H - 8, 2, 8); }
      if ((skipped ? k % step === 0 : col % step === 0) && x - lastLabel > 60) {
        ctx.fillStyle = MUTED; ctx.fillRect(x - 0.5, AXIS_H - 6, 1, 5);
        ctx.fillText(String(col), x, AXIS_H / 2 - 2);
        lastLabel = x;
      }
    }
  }

  // Grey wash over columns outside the compared region (run results).
  function drawCompared(ctx, y, h) {
    var cmp = S.data.compared;
    if (!cmp) return;
    ctx.fillStyle = 'rgba(110,110,110,0.30)';
    var mw = mainW(), step = Math.max(1, S.px);
    for (var x = 0; x < mw; x += step) {
      var k = Math.floor(S.start + x / S.px);
      if (k >= S.shown.length) break;
      if (cmp.charCodeAt(S.shown[k]) === 48) ctx.fillRect(S.px >= 1 ? X(k) : LABEL_W + x, y, Math.ceil(step), h);
    }
  }

  // Marks on the anchor row where columns that are a gap in the anchor were hidden.
  function drawInserts(ctx) {
    ctx.fillStyle = '#7b3294';
    S.inserts.forEach(function(ins) {
      var x = X(ins.at);
      ctx.beginPath(); ctx.moveTo(x - 4, AXIS_H + 1); ctx.lineTo(x + 4, AXIS_H + 1); ctx.lineTo(x, AXIS_H + 8); ctx.fill();
      ctx.fillRect(x - 0.5, AXIS_H + 1, 1, ANCHOR_H);
    });
  }

  // Overview: how many sequences differ from the anchor (SNP or indel) along the whole
  // alignment; red where 3 or more do. The blue box is the part on screen.
  function drawOverview() {
    var ctx = context(el.ov), mw = mainW(), n = S.shown.length, per = n / mw, nd = S.view.ndiff, cmp = S.data.compared;
    ctx.clearRect(0, 0, S.width, OV_H);
    ctx.font = '11px system-ui, sans-serif'; ctx.fillStyle = MUTED; ctx.textAlign = 'left';
    ctx.fillText('Whole alignment: sequences differing', 4, 15);
    ctx.fillText('(click or drag to move)', 4, 30);
    var max = 1;
    S.shown.forEach(function(c) { if (nd[c] > max) max = nd[c]; });
    for (var p = 0; p < mw; p++) {
      var a = Math.floor(p * per), b = Math.max(a + 1, Math.floor((p + 1) * per)), m = 0, outside = !!cmp;
      for (var j = a; j < b && j < n; j++) {
        if (nd[S.shown[j]] > m) m = nd[S.shown[j]];
        if (cmp && cmp.charCodeAt(S.shown[j]) === 49) outside = false;
      }
      if (outside) { ctx.fillStyle = '#eef0f2'; ctx.fillRect(LABEL_W + p, 4, 1, 36); }
      if (m) { var hh = Math.max(2, 36 * m / max); ctx.fillStyle = m >= 3 ? HOT : LINE; ctx.fillRect(LABEL_W + p, 40 - hh, 1, hh); }
    }
    ctx.fillStyle = '#d0d7de'; ctx.fillRect(LABEL_W, 40, mw, 1);
    var x0 = LABEL_W + S.start / per, ww = Math.max(3, mw / S.px / per);
    ctx.globalAlpha = 0.08; ctx.fillStyle = ACCENT; ctx.fillRect(x0, 2, ww, 38); ctx.globalAlpha = 1;
    ctx.strokeStyle = ACCENT; ctx.lineWidth = 2; ctx.strokeRect(x0, 2, ww, 38);
  }

  function drawLegend() {
    var sw = function(fill, ink, text, label) {
      return '<span><i class="sw" style="background:' + fill + ';color:' + ink + '">' + text + '</i>' + label + '</span>';
    };
    el.legend.innerHTML =
      sw('#ffffff', INK[0], '·', 'same as the anchor') + sw(FILL[1], INK[1], 'T', 'SNP') +
      sw(FILL[3], INK[3], 'K', 'ambiguity code, cannot match (counted)') +
      sw(FILL[2], INK[2], 'R', 'ambiguity code, fits (not counted)') + sw(FILL[4], INK[4], '-', 'indel') +
      sw(FILL[5], INK[5], 'N', 'no data (N, missing start or end)') + sw(FILL[7], INK[7], 'A', 'anchor has no base here') +
      (S.data.compared ? sw('rgba(110,110,110,0.30)', '#333', '', 'outside the compared region') : '') +
      '<span><i class="sw" style="background:' + HOT + '"></i>overview: 3 or more sequences differ</span>';
  }

  // ---- Pointer ----------------------------------------------------------------------
  // What is under the pointer: {row (-1 anchor), k (index into shown), part}.
  function at(e) {
    var rect = el.canvas.getBoundingClientRect(), x = e.clientX - rect.left, y = e.clientY - rect.top;
    var res = {x: x, y: y, part: null};
    if (y < AXIS_H) return res;
    if (y < rowsTop()) res.row = -1;
    else {
      var idx = Math.floor(S.top) + Math.floor((y - rowsTop()) / ROW_H);
      if (idx >= nRows()) return res;
      res.row = S.order[idx];
    }
    if (x < LABEL_W) { res.part = 'name'; return res; }
    res.part = 'cell';
    res.k = Math.floor(S.start + (x - LABEL_W) / S.px);
    if (res.k >= S.shown.length) res.part = null;
    return res;
  }

  function wheel(e) {
    if (!S.view) return;
    var rect = el.canvas.getBoundingClientRect(), x = e.clientX - rect.left;
    if (x < LABEL_W) {   // over the names: scroll the rows (or the page)
      if (scrolling()) { e.preventDefault(); S.top += e.deltaY / ROW_H / 2; draw(); }
      return;
    }
    e.preventDefault();
    if (e.shiftKey || Math.abs(e.deltaX) > Math.abs(e.deltaY)) { S.start += (e.deltaX || e.deltaY) / S.px; S.fitted = false; draw(); return; }
    zoom(Math.exp(-e.deltaY * 0.0015), S.start + (x - LABEL_W) / S.px);
  }

  function pointerdown(e) {
    if (!S.view || e.button !== 0) return;
    S.drag = {x: e.clientX, y: e.clientY, start: S.start, top: S.top, moved: false, ev: e};
  }
  function dragmove(e) {
    if (!S.drag) return;
    var dx = e.clientX - S.drag.x, dy = e.clientY - S.drag.y;
    if (Math.abs(dx) + Math.abs(dy) > 3) S.drag.moved = true;
    if (!S.drag.moved) return;
    el.canvas.classList.add('dragging');
    el.tip.style.display = 'none';
    S.start = S.drag.start - dx / S.px;
    S.fitted = false;
    if (scrolling()) S.top = S.drag.top - dy / ROW_H;
    draw();
  }
  function dragend() {
    if (S.drag && !S.drag.moved) click(S.drag.ev);
    S.drag = null;
    el.canvas.classList.remove('dragging');
  }
  function overviewTo(e) {
    var rect = el.ov.getBoundingClientRect(), per = S.shown.length / mainW();
    S.start = (e.clientX - rect.left - LABEL_W) * per - mainW() / S.px / 2;
    S.fitted = false;
    draw();
  }

  // A click on a row chooses it (Previous/Next difference then look only at that row).
  function click(e) {
    var p = at(e);
    if (p.row == null || p.row < 0 || (p.part !== 'name' && p.part !== 'cell')) return;
    S.selected = S.selected === p.row ? -1 : p.row;
    if (p.part === 'cell') S.mark = S.shown[p.k];
    draw();
  }

  function contextmenu(e) {
    var p = at(e);
    if (S.opts.standalone || p.row == null || (p.part !== 'name' && p.part !== 'cell')) return;
    e.preventDefault();
    el.menu.innerHTML = '';
    var add = function(text, name) {
      var item = document.createElement('div');
      item.className = 'msa-menu-item';
      item.textContent = text;
      item.addEventListener('click', function() {
        el.menu.style.display = 'none';
        Shiny.setInputValue('msa_anchor_pick', {name: name}, {priority: 'event'});
      });
      el.menu.appendChild(item);
    };
    if (p.row >= 0) add('Use ' + S.data.names[p.row] + ' as anchor', S.data.names[p.row]);
    if (S.view.anchor != null) add('Compare with the majority consensus', '');
    if (!el.menu.childNodes.length) return;
    var box = el.root.getBoundingClientRect();
    el.menu.style.left = (e.clientX - box.left) + 'px';
    el.menu.style.top = (e.clientY - box.top) + 'px';
    el.menu.style.display = 'block';
  }

  // Base number of column c in a sequence: the count of non-gap letters up to it (as
  // msa_positions() in R/msa.R). Gaps give "gap after N", missing ends say so.
  function position(chars, c) {
    var n = 0, i;
    for (i = 0; i <= c; i++) if (chars[i] !== '-') n++;
    if (chars[c] !== '-') return String(n);
    if (n === 0) return 'missing (before its start)';
    for (i = c + 1; i < chars.length; i++) if (chars[i] !== '-') return 'gap after ' + n;
    return 'missing (after its end)';
  }

  function hover(e) {
    if (!S.view || (S.drag && S.drag.moved)) return;
    var p = at(e);
    if (p.part !== 'cell' && p.part !== 'name') { el.tip.style.display = 'none'; return; }
    var info = rowInfo(p.row), r = p.row, v = S.view, d = S.data, tip = el.tip;
    var metaRow = r >= 0 ? r : v.anchor;
    tip.innerHTML = '';
    var line = function(text, cls) {
      var div = document.createElement('div');
      if (cls) div.className = cls;
      div.textContent = text;
      tip.appendChild(div);
    };
    var head = document.createElement('div'), b = document.createElement('b');
    b.textContent = info.name;
    head.appendChild(b);
    var extra = [];
    if (r < 0) extra.push('anchor');
    if (metaRow != null && d.dates && d.dates.text[metaRow]) extra.push(d.dates.text[metaRow]);
    if (metaRow != null && d.query[metaRow]) extra.push('query');
    if (extra.length) head.appendChild(document.createTextNode(' · ' + extra.join(' · ')));
    tip.appendChild(head);
    if (r >= 0 && d.dropped[r]) line('Left out of the SNP table and tree (covers < 50% of the compared region)', 'msa-tip-warn');
    if (metaRow != null && d.meta[metaRow]) {
      var m = d.meta[metaRow], parts = Object.keys(m).filter(function(k) { return m[k] !== '?' && m[k] !== ''; })
        .map(function(k) { return k + ' ' + m[k]; });
      if (parts.length) line(parts.join(' · '), 'msa-tip-muted');
    }
    if (r >= 0) {
      var st = v.stats, nd = st.n_counted ? (100 * st.nodata[r] / st.n_counted).toFixed(1) : '0';
      line((v.diffs[r] == null ? 'not in the SNP table' : plural(v.diffs[r], 'SNP')) + ' · ' +
           plural(st.indel_events[r], 'indel event') + ' (' + plural(st.indel_cols[r], 'column') + ') · ' +
           plural(st.ambiguity[r], 'ambiguity code') + ' · ' + nd + '% no data' +
           (d.compared ? ' (compared region)' : ''));
    }
    if (p.part === 'cell') {
      var c = S.shown[p.k], kc = info.cls.charCodeAt(c) - 48;
      tip.appendChild(document.createElement('hr'));
      var ins = S.hideGap && S.inserts.filter(function(x) { return Math.abs(X(x.at) - p.x) < 5; })[0];
      if (r < 0 && ins) line(plural(ins.n, 'column') + ' hidden here: a gap in the anchor (extra bases in other sequences)');
      var rc = metaRow != null && d.revcomp[metaRow] ? ' (on the reverse complement)' : '';
      if (r >= 0) {
        var arc = v.anchor != null && d.revcomp[v.anchor] ? ' (on the reverse complement)' : '';
        line('Column ' + d.cols[c] + ' · position ' + position(info.chars, c) + rc +
             ' · in the anchor ' + position(v.anchor_chars, c) + arc);
        line(info.chars[c] + ' / ' + v.anchor_chars[c] + ' (this / anchor): ' + WHAT[kc]);
        if (d.compared && d.compared.charCodeAt(c) === 48) line('Outside the compared region (not in the SNP table)', 'msa-tip-muted');
      } else {
        line('Column ' + d.cols[c] + ' · position ' + position(info.chars, c) + rc + ' · base ' + info.chars[c]);
      }
      line(plural(v.ndiff[c], 'sequence') + ' of ' + v.rows.length + (v.ndiff[c] === 1 ? ' differs' : ' differ') +
           ' from the anchor here', 'msa-tip-muted');
    }
    var right = e.clientX > window.innerWidth / 2, below = e.clientY < window.innerHeight / 2;
    tip.style.left = right ? '' : (e.clientX + 14) + 'px';
    tip.style.right = right ? (window.innerWidth - e.clientX + 14) + 'px' : '';
    tip.style.top = below ? (e.clientY + 14) + 'px' : '';
    tip.style.bottom = below ? '' : (window.innerHeight - e.clientY + 14) + 'px';
    tip.style.display = 'block';
  }

  // What is on screen (1-based, as R counts), for the PDF and FASTA downloads.
  var reportTimer = null;
  function reportWindow(first, n) {
    if (S.opts.standalone || !window.Shiny || !Shiny.setInputValue) return;
    clearTimeout(reportTimer);
    reportTimer = setTimeout(function() {
      var kFirst = Math.floor(S.start), kLast = Math.min(S.shown.length - 1, Math.floor(S.start + mainW() / S.px - 1e-9));
      Shiny.setInputValue('msa_window', {
        cols: S.shown.slice(kFirst, kLast + 1).map(function(c) { return c + 1; }),
        rows: S.order.slice(first, first + n).map(function(r) { return r + 1; }),
        order: S.order.map(function(r) { return r + 1; }), sort: S.sort
      });
    }, 250);
  }

  function savePicture(name) {
    if (!S.view) return;
    el.canvas.toBlob(function(blob) {
      var a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = name;
      document.body.appendChild(a);
      a.click();
      setTimeout(function() { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
    });
  }

  if (window.jQuery) {
    jQuery(document).on('shiny:connected', function() {
      Shiny.addCustomMessageHandler('msa_data', function(d) { setData(d); });
      Shiny.addCustomMessageHandler('msa_view', function(v) { setView(v); });
      Shiny.addCustomMessageHandler('msa_options', setOptions);
      Shiny.addCustomMessageHandler('msa_png', savePicture);
    });
  }
  window.fhiMsa = {
    state: S, draw: draw, act: act,
    // Standalone page (msa_html() in R/msa.R): data and classes are in the file.
    load: function(d, v, opts) { S.opts = opts || {}; init(); setData(d); setView(v); }
  };
})();
