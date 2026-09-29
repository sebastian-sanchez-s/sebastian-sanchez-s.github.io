import os
from nbconvert import HTMLExporter
import nbformat

# Palette: background #000 | cyan #00ffff | darkcyan #008b8b | text cream #f5f4ee
BRAND_CSS = """
<style>
 
  /* Page */
  html, body {
    background-color: #000000 !important;
    color: #f5f4ee !important;
  }
  body {
    margin: 0 10%;
    font-family: 'Styrene A', 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
    font-size: 18px;
  }
 
  /* Paper-page layout (~US Letter). Margins live on the page, not between cells */
  .jp-Notebook, #notebook-container, .container {
    box-sizing: border-box;
    background-color: #000000 !important;
    border: 1px solid rgba(0, 139, 139, 0.45);
    box-shadow: none !important;
  }
  .jp-Cell, .cell { margin-bottom: 0.35rem !important; }
  .jp-Cell:last-child, .cell:last-child { margin-bottom: 0 !important; }
 
  /* Markdown */
  .jp-RenderedMarkdown, .jp-RenderedHTMLCommon, .text_cell_render,
  .jp-RenderedMarkdown p, .jp-RenderedMarkdown li {
    font-family: 'Anthropic Serif', 'Tiempos Text', Georgia, 'Times New Roman', serif;
    color: #f5f4ee !important;
    line-height: 1.7;
  }
 
  h1, h2, h3 {
    font-family: 'Tiempos Headline', Georgia, 'Times New Roman', serif;
    font-weight: 500;
  }
  h1 {
    font-size: 2rem;
    color: #00ffff !important;
    border-bottom: 1px solid #008b8b;
    padding-bottom: 0.4em;
    margin-top: 2rem;
    margin-bottom: 1rem;
  }
  h2 {
    font-size: 1.5rem;
    color: #00ffff !important;
    margin-top: 1.75rem;
    margin-bottom: 0.75rem;
  }
  h3 {
    font-size: 1.2rem;
    color: #008b8b !important;
    filter: brightness(1.4);
    margin-top: 1.5rem;
  }
 
  .jp-RenderedMarkdown a, .text_cell_render a {
    color: #00ffff !important;
    text-decoration: none;
    border-bottom: 1px solid rgba(0, 255, 255, 0.4);
  }
  .jp-RenderedMarkdown code, .text_cell_render code {
    background-color: #0f1a1a !important;
    color: #00ffff !important;
    padding: 0.15em 0.45em;
    border-radius: 4px;
    font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace;
    font-size: 0.85em;
  }
  blockquote {
    color: #d8d5c8 !important;
    background-color: #0a0a0a;
    border-left: 4px solid #008b8b;
    padding: 0.75em 1.25em;
    margin: 0.75rem 0;
  }
  ul li::marker, ol li::marker { color: #008b8b; }
  hr {
    border: none;
    border-top: 1px solid rgba(0, 139, 139, 0.6);
    margin: 1.5rem 0;
  }
 
  /* Code cells */
  .jp-CodeCell .jp-Cell-inputWrapper, div.input_area {
    background-color: #0a0a0a !important;
    border: 1px solid #008b8b !important;
    border-radius: 6px !important;
    overflow: hidden;
  }
  .jp-InputArea-editor, div.input_area pre, .highlight {
    background-color: #0a0a0a !important;
    color: #f5f4ee !important;
    font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace !important;
    font-size: 0.88em;
  }
 
  /* Pygments syntax highlighting */
  .highlight .k, .highlight .kn, .highlight .kd, .highlight .kc, .highlight .kr { color: #00ffff !important; }
  .highlight .nf, .highlight .fm, .highlight .nc { color: #00ffff !important; }
  .highlight .nb, .highlight .bp { color: #20b2aa !important; }
  .highlight .n, .highlight .nn, .highlight .p { color: #f5f4ee !important; }
  .highlight .s, .highlight .s1, .highlight .s2, .highlight .sa { color: #e8dcc0 !important; }
  .highlight .mi, .highlight .mf, .highlight .mh { color: #40e0d0 !important; }
  .highlight .o, .highlight .ow { color: #f5f4ee !important; }
  .highlight .c, .highlight .c1, .highlight .cm { color: #7d7b72 !important; font-style: italic; }
 
  /* Prompts: In / Out */
  .jp-InputPrompt, div.input_prompt { color: #00ffff !important; }
  .jp-OutputPrompt, div.output_prompt { color: #008b8b !important; }
  .jp-InputPrompt, .jp-OutputPrompt, div.prompt {
    font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace !important;
    font-size: 0.8rem;
  }
 
  /* Outputs: same as canvas so plots blend into the page */
  .jp-OutputArea-output, div.output_area, .jp-RenderedText {
    background-color: #000000 !important;
    color: #f5f4ee !important;
  }
  .jp-OutputArea-output pre, div.output_subarea pre {
    color: #d8d5c8 !important;
    font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace;
  }
  .jp-OutputArea-output img, div.output_area img {
    display: block;
    margin: 1.25rem auto;
    max-width: 88%;
    background-color: transparent;
    border: none;
    box-shadow: none;
  }
 
  /* Tables */
  table.dataframe, .jp-RenderedHTMLCommon table {
    background-color: #000000 !important;
    color: #f5f4ee !important;
    border-collapse: collapse;
    margin: 1rem auto;
  }
  table.dataframe th, .jp-RenderedHTMLCommon th {
    background-color: #0a0a0a !important;
    color: #00ffff !important;
    border-bottom: 2px solid #008b8b !important;
    text-transform: uppercase;
    font-size: 0.75rem;
    letter-spacing: 0.05em;
    font-weight: 500;
    padding: 0.6em 1em;
  }
  table.dataframe td, .jp-RenderedHTMLCommon td {
    border: none;
    border-bottom: 1px solid rgba(0, 139, 139, 0.35) !important;
    padding: 0.5em 1em;
  }
  table.dataframe tr:hover td { background-color: rgba(0, 139, 139, 0.15); }
 
  /* Error tracebacks */
  .jp-OutputArea-output.jp-mod-error, div.output_area.output_error {
    background-color: rgba(224, 119, 111, 0.08) !important;
    border-left: 4px solid #e0776f;
    padding: 1rem 1.25rem;
  }
 
  /* Scrollbars */
  ::-webkit-scrollbar { width: 10px; height: 10px; }
  ::-webkit-scrollbar-track { background: #000000; }
  ::-webkit-scrollbar-thumb { background: #008b8b; border-radius: 6px; }
  ::-webkit-scrollbar-thumb:hover { background: #00ffff; }
 
</style>
"""
DARK_CSS = """
<style>
  @import url('https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;700&family=Roboto+Mono:wght@400;500&display=swap');

  :root {
    /* Material Design dark theme baseline */
    --md-bg:         #000000;   /* true black canvas — lets images blend seamlessly */
    --md-surface:    #121212;   /* elevation 0dp surface */
    --md-surface1:   #1d1d1d;   /* +1dp overlay (~5%) */
    --md-surface2:   #212121;   /* +2dp overlay (~7%) */
    --md-surface4:   #272727;   /* +4dp overlay (~9%) */
    --md-surface8:   #2c2c2c;   /* +8dp overlay (~12%) */
    --md-divider:    rgba(255,255,255,0.12);
    --md-divider-strong: rgba(255,255,255,0.20);

    /* On-dark emphasis levels (Material spec) */
    --md-on-high:    rgba(255,255,255,0.87);
    --md-on-medium:  rgba(255,255,255,0.60);
    --md-on-disabled:rgba(255,255,255,0.38);

    /* Material accent roles */
    --md-primary:    #BB86FC;
    --md-primary-var:#3700B3;
    --md-secondary:  #03DAC6;
    --md-error:      #CF6679;
    --md-amber:      #FFD54F;
    --md-orange:     #FFB74D;
    --md-blue:       #82B1FF;
    --md-green:      #A5D6A7;

    /* Map onto JupyterLab / nbconvert vars */
    --jp-layout-color0: var(--md-bg) !important;
    --jp-layout-color1: var(--md-surface) !important;
    --jp-layout-color2: var(--md-surface1) !important;
    --jp-layout-color3: var(--md-surface2) !important;
    --jp-layout-color4: var(--md-surface4) !important;

    --jp-ui-font-color0: var(--md-on-high) !important;
    --jp-ui-font-color1: var(--md-on-medium) !important;
    --jp-ui-font-color2: var(--md-on-disabled) !important;

    --jp-content-font-color0: var(--md-on-high) !important;
    --jp-content-font-color1: var(--md-on-high) !important;
    --jp-content-link-color: var(--md-secondary) !important;

    --jp-border-color0: var(--md-divider) !important;
    --jp-border-color1: var(--md-divider) !important;
    --jp-border-color2: var(--md-divider-strong) !important;
    --jp-border-color3: var(--md-divider-strong) !important;

    --jp-cell-editor-background: var(--md-surface1) !important;
    --jp-cell-editor-border-color: var(--md-divider) !important;
    --jp-cell-prompt-not-active-font-color: var(--md-on-disabled) !important;

    --jp-rendermime-error-background: rgba(207, 102, 121, 0.12) !important;
    --jp-warn-color0: var(--md-amber) !important;
    --jp-error-color0: var(--md-error) !important;

    --jp-brand-color1: var(--md-primary) !important;
    --jp-accent-color1: var(--md-secondary) !important;
  }

  html, body {
    background-color: var(--md-bg) !important;
    color: var(--md-on-high) !important;
  }

  body {
    margin: 0 10%;
  }

  body, .jp-Notebook, #notebook-container, .container {
    background-color: var(--md-bg) !important;
    color: var(--md-on-high) !important;
    font-family: "Roboto", "Segoe UI", Arial, sans-serif !important;
    font-size: 16px;
    letter-spacing: 0.01em;
  }

  /* Layout: fixed "paper page" proportions (~US Letter), margins on the page itself, not between cells */
  #notebook-container, .container {
    max-width: 816px;           /* ~8.5in at 96dpi */
    min-height: 1056px;         /* ~11in at 96dpi, so short notebooks still read as a page */
    margin: 3rem auto;
    padding: 96px 72px;         /* ~1in top/bottom, 0.75in sides — like real page margins */
    background-color: var(--md-bg) !important;
    box-sizing: border-box;
  }

  /* Cells: tight, document-like flow — no artificial gaps beyond nbconvert's own spacing */
  .jp-Cell, .cell {
    margin-bottom: 0.35rem !important;
  }
  .jp-Cell:last-child, .cell:last-child {
    margin-bottom: 0 !important;
  }

  /* Typography scale (Material type system) */
  .jp-RenderedMarkdown, .text_cell_render {
    color: var(--md-on-high) !important;
    line-height: 1.7;
    font-weight: 400;
  }
  .jp-RenderedMarkdown h1, .text_cell_render h1 {
    font-size: 2.0rem;
    font-weight: 500;
    color: var(--md-on-high) !important;
    border-bottom: 1px solid var(--md-divider);
    padding-bottom: 0.5em;
    margin-top: 2rem;
    margin-bottom: 1rem;
  }
  .jp-RenderedMarkdown h2, .text_cell_render h2 {
    font-size: 1.5rem;
    font-weight: 500;
    color: var(--md-primary) !important;
    margin-top: 1.75rem;
    margin-bottom: 0.75rem;
  }
  .jp-RenderedMarkdown h3, .text_cell_render h3 {
    font-size: 1.2rem;
    font-weight: 500;
    color: var(--md-secondary) !important;
    margin-top: 1.5rem;
  }
  .jp-RenderedMarkdown a, .text_cell_render a {
    color: var(--md-secondary) !important;
    text-decoration: none;
    border-bottom: 1px solid rgba(3, 218, 198, 0.4);
  }
  .jp-RenderedMarkdown a:hover, .text_cell_render a:hover {
    border-bottom-color: var(--md-secondary);
  }
  .jp-RenderedMarkdown code, .text_cell_render code {
    background-color: var(--md-surface2) !important;
    color: var(--md-amber) !important;
    padding: 0.15em 0.45em;
    border-radius: 4px;
    font-family: "Roboto Mono", "Consolas", monospace;
    font-size: 0.85em;
  }
  .jp-RenderedMarkdown blockquote, .text_cell_render blockquote {
    border-left: 3px solid var(--md-primary);
    background-color: var(--md-surface1);
    color: var(--md-on-medium) !important;
    padding: 0.75em 1.25em;
    border-radius: 0 8px 8px 0;
    margin: 0.75rem 0;
  }
  .jp-RenderedMarkdown strong, .text_cell_render strong { color: var(--md-on-high); font-weight: 700; }
  .jp-RenderedMarkdown em, .text_cell_render em { color: var(--md-blue); }
  .jp-RenderedMarkdown hr, .text_cell_render hr {
    border: none;
    border-top: 1px solid var(--md-divider);
    margin: 1.5rem 0;
  }

  /* Code cells — Material "card" surface with elevation, rounded corners */
  .jp-CodeCell .jp-Cell-inputWrapper,
  div.input_area {
    background-color: var(--md-surface1) !important;
    border: none !important;
    border-radius: 8px !important;
    box-shadow: 0 1px 3px rgba(0,0,0,0.5), 0 1px 2px rgba(0,0,0,0.4);
    overflow: hidden;
    margin: 0.25rem 0;
  }
  .jp-InputArea-editor, div.input_area pre {
    background-color: var(--md-surface1) !important;
    font-family: "Roboto Mono", "Consolas", "Menlo", monospace !important;
    font-size: 0.88em;
    padding: 0.5rem 0.25rem;
  }

  /* Pygments syntax highlighting — Material accent mapping */
  .highlight .k,  .highlight .kn, .highlight .kd { color: var(--md-primary) !important; }     /* keywords */
  .highlight .n  { color: var(--md-on-high) !important; }                                     /* names */
  .highlight .s,  .highlight .s1, .highlight .s2 { color: var(--md-green) !important; }        /* strings */
  .highlight .mi, .highlight .mf { color: var(--md-orange) !important; }                       /* numbers */
  .highlight .c,  .highlight .c1, .highlight .cm { color: var(--md-on-disabled) !important; font-style: italic; }
  .highlight .o,  .highlight .ow { color: var(--md-secondary) !important; }                    /* operators */
  .highlight .nf, .highlight .fm { color: var(--md-blue) !important; }                         /* functions */
  .highlight .nb  { color: var(--md-amber) !important; }                                       /* builtins */
  .highlight .bp  { color: var(--md-error) !important; }                                       /* self/cls */
  .highlight { background-color: var(--md-surface1) !important; }

  /* Output areas — same as canvas so plots melt into the page */
  .jp-OutputArea-output, div.output_area {
    background-color: var(--md-bg) !important;
    color: var(--md-on-high) !important;
  }
  .jp-OutputArea-output pre, div.output_subarea pre {
    color: var(--md-on-medium) !important;
    font-family: "Roboto Mono", monospace;
  }

  /* Figure output — no border/background/shadow so black-bg plots blend into the page */
  .jp-OutputArea-output img, div.output_area img {
    display: block;
    margin: 1.25rem auto;
    max-width: 88%;
    background-color: transparent;
    border: none;
    box-shadow: none;
  }

  /* Tables — Material data-table styling */
  table.dataframe, .jp-RenderedHTMLCommon table {
    background-color: var(--md-surface) !important;
    color: var(--md-on-high) !important;
    border-collapse: collapse;
    border-radius: 8px;
    overflow: hidden;
    margin: 1rem auto;
  }
  table.dataframe th, .jp-RenderedHTMLCommon th {
    background-color: var(--md-surface2) !important;
    color: var(--md-on-medium) !important;
    border-bottom: 2px solid var(--md-divider-strong) !important;
    text-transform: uppercase;
    font-size: 0.75rem;
    letter-spacing: 0.05em;
    font-weight: 500;
    padding: 0.6em 1em;
  }
  table.dataframe td, .jp-RenderedHTMLCommon td {
    border: none;
    border-bottom: 1px solid var(--md-divider) !important;
    padding: 0.5em 1em;
  }
  table.dataframe tr:hover td {
    background-color: rgba(255,255,255,0.04);
  }

  /* Error tracebacks — Material error card */
  .jp-OutputArea-output.jp-mod-error, div.output_area.output_error {
    background-color: rgba(207, 102, 121, 0.08) !important;
    border-left: 4px solid var(--md-error);
    border-radius: 0 8px 8px 0;
    padding: 1rem 1.25rem;
    font-family: "Roboto Mono", monospace;
  }

  /* Cell execution prompts */
  .jp-InputPrompt, .jp-OutputPrompt, div.prompt {
    color: var(--md-on-disabled) !important;
    font-family: "Roboto Mono", monospace !important;
    font-size: 0.8rem;
  }

  /* Scrollbars */
  ::-webkit-scrollbar { width: 10px; height: 10px; }
  ::-webkit-scrollbar-track { background: var(--md-bg); }
  ::-webkit-scrollbar-thumb { background: var(--md-surface4); border-radius: 6px; }
  ::-webkit-scrollbar-thumb:hover { background: var(--md-surface8); }
</style>
"""


def export_notebook_with_catppuccin(notebook_path, output_path=None):
    if not output_path:
        output_path = notebook_path.replace(".ipynb", ".html")

    print(f"Reading {notebook_path}...")
    with open(notebook_path, "r", encoding="utf-8") as f:
        nb = nbformat.read(f, as_version=4)

    html_exporter = HTMLExporter()
    html_exporter.theme = "dark"  # base dark structure; our CSS overrides colors on top

    print("Generating base HTML...")
    (body, resources) = html_exporter.from_notebook_node(nb)

    CSS = BRAND_CSS

    print("Injecting Catppuccin styling...")
    if "</head>" in body:
        modified_body = body.replace("</head>", f"{CSS}\n</head>")
    else:
        modified_body = f"{CSS}\n{body}"

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(modified_body)

    print(f"Successfully created: {output_path}")

if __name__ == "__main__":
    import sys
    NOTEBOOK_FILE = sys.argv[1]

    if os.path.exists(NOTEBOOK_FILE):
        export_notebook_with_catppuccin(NOTEBOOK_FILE)
    else:
        print(f"Error: {NOTEBOOK_FILE} not found. Please update the filename in the script.")
