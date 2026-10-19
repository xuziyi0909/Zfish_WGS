"""Build the tutorial slide deck as one self-contained web page for GitHub Pages.

Input : tutorial/deck.json (title and slide order) and tutorial/slides/<id>.html,
        one <section> per slide on a 1920x1080 canvas with inline styles.
Output: docs/slides/index.html, served at https://xuziyi0909.github.io/Zfish_WGS/slides/

Usage (from the pipeline folder):  python3 tutorial/build_slides.py

The slides were written for Claude's Slides format, so the few custom elements in it
(<x-icon>, <x-shape>) are converted to inline SVG here, and uploaded images are mapped
to the copies in docs/assets/.
"""

import html
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent / "docs" / "slides" / "index.html"

# uploaded image -> file served by the guide (relative to docs/slides/)
IMAGES = {
    "/_blob/4f61db685ca08aa4ca0d5182233b7278": "../assets/variant_density.png",
    "/_blob/0cf4213a8e603b7452fb33398465eb01": "../assets/homozygosity.png",
}

# 24x24 line icons, drawn in the icon's text colour
ICONS = {
    "Check": '<polyline points="20 6 9 17 4 12"/>',
    "Database": '<ellipse cx="12" cy="5" rx="8" ry="3"/><path d="M4 5v14c0 1.7 3.6 3 8 3s8-1.3 8-3V5"/>'
                '<path d="M4 12c0 1.7 3.6 3 8 3s8-1.3 8-3"/>',
    "Key": '<circle cx="7.5" cy="15.5" r="4.5"/><path d="M10.7 12.3 21 2"/><path d="m16 7 3 3"/>'
           '<path d="m18.5 4.5 2 2"/>',
    "Code": '<polyline points="16 18 22 12 16 6"/><polyline points="8 6 2 12 8 18"/>',
}


def style_value(style, prop):
    m = re.search(rf"(?:^|;)\s*{prop}\s*:\s*([^;]+)", style)
    return m.group(1).strip() if m else ""


def px(value):
    return float(value.replace("px", "")) if value else 0.0


def icon_svg(match):
    name, style = match.group(1), match.group(2)
    if name not in ICONS:
        raise SystemExit(f"No SVG for <x-icon name=\"{name}\">: add it to ICONS in {__file__}")
    return (f'<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" '
            f'stroke-linecap="round" stroke-linejoin="round" style="flex:none;{style}" '
            f'aria-hidden="true">{ICONS[name]}</svg>')


def arrow_svg(match):
    kind, style = match.group(1), match.group(2)
    if kind != "arrow-right":
        raise SystemExit(f"<x-shape kind=\"{kind}\"> is not supported by {__file__}")
    w, h = px(style_value(style, "width")), px(style_value(style, "height"))
    head = w * 0.4  # the head is the last 40% of the length, as in the Slides format
    pts = f"0,{h * .28:g} {w - head:g},{h * .28:g} {w - head:g},0 {w:g},{h / 2:g} " \
          f"{w - head:g},{h:g} {w - head:g},{h * .72:g} 0,{h * .72:g}"
    return (f'<svg width="{w:g}" height="{h:g}" viewBox="0 0 {w:g} {h:g}" style="flex:none" '
            f'aria-hidden="true"><polygon points="{pts}" fill="{style_value(style, "background")}"/></svg>')


def convert(section):
    section = re.sub(r'<x-icon name="(\w+)" style="([^"]*)"></x-icon>', icon_svg, section)
    section = re.sub(r'<x-shape kind="([\w-]+)" style="([^"]*)"></x-shape>', arrow_svg, section)
    for blob, path in IMAGES.items():
        section = section.replace(blob, path)
    if re.search(r"<x-|/_blob/", section):
        raise SystemExit(f"Unconverted element left in:\n{section[:200]}")
    return section.replace("<section ", '<section class="slide" ', 1)


CSS = """
html,body{margin:0;height:100%;background:#0e141b;overflow:hidden}
#stage{position:absolute;left:50%;top:50%;width:1920px;height:1080px;transform-origin:center center}
#stage *{box-sizing:border-box;margin:0}
#stage section.slide{position:absolute;inset:0;width:1920px;height:1080px;overflow:hidden}
#stage section.slide:not(.active){display:none!important}
#stage section.slide.active{animation:fade .3s ease-out}
@keyframes fade{from{opacity:0}}
#stage section div{display:flex;flex-direction:column;min-width:0}
#stage h1{font-size:96px;font-weight:600;line-height:1.1}
#stage h2{font-size:64px;font-weight:600;line-height:1.15}
#stage h3{font-size:44px;font-weight:600;line-height:1.2}
#stage p{font-size:32px;line-height:1.4}
#stage p[style*="JetBrains"]{white-space:pre-wrap}
#stage ul,#stage ol{padding-left:1.3em}
#stage table{border-collapse:collapse;width:100%}
#stage th,#stage td{padding:.35em .6em;border-bottom:1px solid rgba(22,32,43,.14);text-align:left;vertical-align:top}
#stage th{font-weight:600}
#stage hr{border:0;height:0;align-self:stretch}
#stage aside{display:none}
#bar{position:fixed;right:16px;bottom:12px;display:flex;gap:8px;align-items:center;
  font:14px system-ui,sans-serif;color:#c3c2b7;opacity:.35;transition:opacity .2s}
#bar:hover,#bar:focus-within{opacity:1}
#bar button,#bar a{background:#232d38;color:#f2f0ea;border:1px solid #3a4756;border-radius:6px;
  padding:6px 10px;font:inherit;cursor:pointer;text-decoration:none}
#notes{position:fixed;left:0;right:0;bottom:0;max-height:38vh;overflow:auto;display:none;
  background:rgba(14,20,27,.94);color:#f2f0ea;font:18px/1.5 system-ui,sans-serif;padding:16px 24px 56px}
#notes.show{display:block}
@media print{
  @page{size:1920px 1080px;margin:0}
  html,body{height:auto;overflow:visible;background:none}
  #stage{position:static;transform:none!important;width:auto;height:auto}
  #stage section.slide{position:relative;display:flex!important;animation:none;break-after:page}
  #bar,#notes{display:none!important}
}
"""

JS = """
const slides = [...document.querySelectorAll('section.slide')];
const stage = document.getElementById('stage');
const notes = document.getElementById('notes');
const counter = document.getElementById('counter');
let current = 0;

function fit() {
  const s = Math.min(innerWidth / 1920, innerHeight / 1080);
  stage.style.transform = `translate(-50%, -50%) scale(${s})`;
}
function show(n) {
  current = Math.max(0, Math.min(slides.length - 1, n));
  slides.forEach((s, k) => s.classList.toggle('active', k === current));
  counter.textContent = `${current + 1} / ${slides.length}`;
  const aside = slides[current].querySelector('aside');
  notes.textContent = aside ? aside.textContent : 'No speaker notes for this slide.';
  history.replaceState(null, '', '#' + (current + 1));
}
const next = () => show(current + 1), prev = () => show(current - 1);

document.addEventListener('keydown', e => {
  if (e.altKey || e.ctrlKey || e.metaKey) return;
  if (['ArrowRight', 'ArrowDown', 'PageDown', ' ', 'Enter'].includes(e.key)) { next(); e.preventDefault(); }
  else if (['ArrowLeft', 'ArrowUp', 'PageUp', 'Backspace'].includes(e.key)) { prev(); e.preventDefault(); }
  else if (e.key === 'Home') show(0);
  else if (e.key === 'End') show(slides.length - 1);
  else if (e.key === 'n' || e.key === 'N') notes.classList.toggle('show');
  else if (e.key === 'f' || e.key === 'F') document.fullscreenElement ? document.exitFullscreen() : document.documentElement.requestFullscreen();
});
stage.addEventListener('click', e => { if (!e.target.closest('a') && !getSelection().toString()) next(); });
let touchX = null;
addEventListener('touchstart', e => { touchX = e.touches[0].clientX; }, {passive: true});
addEventListener('touchend', e => {
  if (touchX === null) return;
  const dx = e.changedTouches[0].clientX - touchX;
  if (Math.abs(dx) > 40) (dx < 0 ? next : prev)();
  touchX = null;
});
document.getElementById('prev').onclick = prev;
document.getElementById('next').onclick = next;
document.getElementById('notesbtn').onclick = () => notes.classList.toggle('show');
addEventListener('resize', fit);
fit();
show((parseInt(location.hash.slice(1), 10) || 1) - 1);
"""

FONTS = ("https://fonts.googleapis.com/css2?family=Source+Serif+4:wght@400..700"
         "&family=IBM+Plex+Sans:wght@400;600&family=JetBrains+Mono:wght@400;600&display=swap")


def main():
    deck = json.loads((ROOT / "deck.json").read_text())
    sections = [convert((ROOT / "slides" / f"{sid}.html").read_text().strip()) for sid in deck["order"]]
    title = html.escape(deck["title"])
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="{FONTS}">
<style>{CSS}</style>
</head>
<body>
<!-- Generated by tutorial/build_slides.py from tutorial/slides/: edit those files, not this one. -->
<main id="stage" aria-label="{title}">
{chr(10).join(sections)}
</main>
<div id="notes" aria-live="polite"></div>
<nav id="bar" aria-label="Slide controls">
<a href="../">Guide</a>
<button id="prev" aria-label="Previous slide">&#8592;</button>
<span id="counter"></span>
<button id="next" aria-label="Next slide">&#8594;</button>
<button id="notesbtn" title="Speaker notes (N)">Notes</button>
</nav>
<script>{JS}</script>
</body>
</html>
"""
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(page)
    print(f"wrote {OUT} ({len(sections)} slides)")


if __name__ == "__main__":
    main()
