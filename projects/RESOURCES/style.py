## Styling
from cycler import cycler

# Manim-inspired color palette (the colors 3b1b actually uses)
manim_colors = {
    "blue":   "#58C4DD",
    "yellow": "#FFFF00",
    "green":  "#83C167",
    "red":    "#FC6255",
    "purple": "#9A72AC",
    "orange": "#FF862F",
    "teal":   "#5CD0B3",
    "pink":   "#EC92B0",
}

def init_plt(plt):
    plt.rcParams.update({
    # --- Backgrounds ---
    "figure.facecolor": "#000000",
    "axes.facecolor":   "#000000",
    "savefig.facecolor": "#000000",
    
    # --- Text & labels ---
    "text.color":       "#FFFFFF",
    "axes.labelcolor":  "#FFFFFF",
    "axes.titlecolor":  "#FFFFFF",
    "xtick.color":      "#FFFFFF",
    "ytick.color":      "#FFFFFF",
    "font.family":      "sans-serif",
    "font.sans-serif":  ["CMU Sans Serif", "DejaVu Sans", "Helvetica", "Arial"],
    "font.size":        10,
    
    # --- Axes / spines (Manim keeps things minimal) ---
    "axes.edgecolor":   "#FFFFFF",
    "axes.linewidth":   1,
    "axes.grid":        False,
    "axes.spines.top":  False,
    "axes.spines.right": False,
    
    # --- Grid, if you turn it on ---
    "grid.color":       "#333333",
    "grid.linestyle":   "--",
    "grid.linewidth":   0.5,
    
    # --- Lines ---
    "lines.linewidth":  1.5,
    "lines.solid_capstyle": "round",
    
    # --- Color cycle (order matters — matches typical 3b1b usage) ---
    "axes.prop_cycle": cycler(color=[
        manim_colors["blue"],
        manim_colors["yellow"],
        manim_colors["green"],
        manim_colors["red"],
        manim_colors["purple"],
        manim_colors["orange"],
        manim_colors["teal"],
        manim_colors["pink"],
    ]),
    
    # --- Legend ---
    "legend.frameon":   False,
    "legend.facecolor": "#000000",
    "legend.labelcolor": "#FFFFFF",
    
    # --- Figure ---
    "figure.figsize": (10, 10),
    "figure.dpi": 100,
    
    # --- Handle 3d ---
    "axes3d.xaxis.panecolor": "#000000",
    "axes3d.yaxis.panecolor": "#000000",
    "axes3d.zaxis.panecolor": "#000000",
    
    # --- Layout management (this fixes the cropping) ---
    "figure.constrained_layout.use": True,  # auto-adjusts spacing so labels/titles aren't cut off
    "figure.autolayout": False,             # don't use this alongside constrained_layout — they conflict
    
    # --- Saving (cropping often only shows up on savefig, not on-screen) ---
    "savefig.bbox": "tight",       # trims the saved figure to content, no clipped labels
    "savefig.pad_inches": 0.1,     # small margin so it's not *too* tight
    "savefig.dpi": 200,            # separate from display dpi — higher res on export only
    })

def style_3d_axes(ax, pane_color="#000000", grid_color="#333333", edge_color="#FFFFFF"):
    """Apply Manim/3b1b-style dark theme to a 3D axes object."""
    # Pane (wall) fills — set alpha=0 too if you want fully transparent walls
    ax.xaxis.set_pane_color((*_hex_to_rgb(pane_color), 1.0))
    ax.yaxis.set_pane_color((*_hex_to_rgb(pane_color), 1.0))
    ax.zaxis.set_pane_color((*_hex_to_rgb(pane_color), 1.0))

    # Pane edges (the lines outlining each wall)
    ax.xaxis.pane.set_edgecolor(edge_color)
    ax.yaxis.pane.set_edgecolor(edge_color)
    ax.zaxis.pane.set_edgecolor(edge_color)
    ax.xaxis.pane.set_alpha(1.0)
    ax.yaxis.pane.set_alpha(1.0)
    ax.zaxis.pane.set_alpha(1.0)

    # Grid lines on the panes
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis._axinfo["grid"]["color"] = grid_color
        axis._axinfo["grid"]["linewidth"] = 0.5

    # Tick and label colors (these usually DO pick up rcParams, but just in case)
    ax.tick_params(colors="#FFFFFF")
    ax.xaxis.label.set_color("#FFFFFF")
    ax.yaxis.label.set_color("#FFFFFF")
    ax.zaxis.label.set_color("#FFFFFF")

    # Axis lines themselves
    ax.xaxis.line.set_color(edge_color)
    ax.yaxis.line.set_color(edge_color)
    ax.zaxis.line.set_color(edge_color)


def _hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i+2], 16) / 255 for i in (0, 2, 4))