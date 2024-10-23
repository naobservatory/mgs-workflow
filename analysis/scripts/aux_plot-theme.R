#==============================================================================
# THEME INFORMATION FOR GGPLOT2 FIGURES
#==============================================================================

suppressMessages(library(ggplot2))

#------------------------------------------------------------------------------
# Auxiliary functions
#------------------------------------------------------------------------------

# Avoid always having to specify units
mmargin <- purrr::partial(ggplot2::margin, unit="mm")
lmargin <- purrr::partial(ggplot2::margin, unit="lines")
len <- purrr::partial(ggplot2::unit, units="mm")
lines <- purrr::partial(ggplot2::unit, units="lines")

# Avoid having to write "inherit.blank = TRUE" everywhere
element_text <- purrr::partial(ggplot2::element_text, inherit.blank = TRUE)
element_line <- purrr::partial(ggplot2::element_line, inherit.blank = TRUE)
element_rect <- purrr::partial(ggplot2::element_rect, inherit.blank = TRUE)

# Nice axis arrows (if needed)
axis_arrow <- purrr::partial(arrow, angle = 15, length = len(3.5),
                             type = "closed")

#------------------------------------------------------------------------------
# Specify font information
#------------------------------------------------------------------------------

# Family
font <- "sans" # Figure font family

# Sizes
fontsize_base <- 9 # Basic figure font size
fontsize_title <- 9 # Axis-title font size
fontsize_label <- 10 # Font size for subfigure legends

# Function for scaling geom sizes to match theme sizes
rescale_font <- function(size_pts) size_pts * 5/14
fontsize_base_rescaled <- rescale_font(fontsize_base)

#------------------------------------------------------------------------------
# Specify parameters
#------------------------------------------------------------------------------

# Core element parameters
line_size_base <- 0.3

# Legend parameters
legend_item_padding <- 0.5 # Extra spacing between adjacent legend items (in mm)
legend_box_fill <- "grey97"
legend_box_margin <- lmargin(rep(0.5, 4))
legend_series_spacing <- 0.2 # (in lines)
legend_key_size <- lines(1.2)
legend_box_spacing <- lines(1)

# Text parameters (in lines by default)
axis_title_spacing_x <- 0.55 # Space between x-axis text and title
axis_title_spacing_y <- 0.55 # Space between y-axis text and title
axis_text_spacing <- 0.25 # Space between axis text and ticks
axis_tick_length <- lines(0.15) # Trying doing this in lines; might go back to mm

# Panel parameters
panel_spacing <- lines(1) # Provisional, tinker later

# Strip parameters
strip_margin <- lmargin(rep(0.5, 4)) # Provisional, tinker later
strip_panel_spacing <- lines(0.2) # Ditto

# High-level plot parameters
plot_title_spacing <- 0.7 # Space between plot titles and main plot (in lines)
plot_tag_margin <- lmargin(b=0.6, r = 0.25)
plot_margin <- lmargin(0.4, 0.4, 0.4, 0.1)

# Additional parameters for internal-legend themes
legend_internal_plain_fill <- alpha("white", 0.5)
legend_internal_text_spacing <- 1.5 # (in lines)
legend_internal_title_margin <- lmargin(b=0.5)
legend_internal_box_margin_scale <- 0.86
legend_internal_spacing_y <- lines(0)
legend_internal_just_default <- c("right", "top")
legend_internal_position_default <- c(0.95, 0.95)

#------------------------------------------------------------------------------
# Define incomplete subthemes (sometimes you only want one)
#------------------------------------------------------------------------------

theme_core <- theme( # Define core elements for rest of theme
  text = element_text(family = font, face = "plain", colour = "black",
                      size = fontsize_base, hjust = 0.5, vjust = 0.5,
                      angle = 0, lineheight = 1, margin = lmargin(),
                      debug = FALSE),
  title = element_text(size = fontsize_title), # Inherits from text
  line = element_line(colour = "black", linewidth = line_size_base, linetype = 1,
                      lineend = "butt", arrow = FALSE),
  rect = element_rect(fill = "white", colour = NA, linewidth = 0.5,
                      linetype = 1),
)

theme_legend <- theme( # Legend formatting
  legend.text = element_text(margin = lmargin(r=legend_item_padding)), # Inherits from text
  legend.title = element_text(face = "bold", vjust = 0.5, hjust = 0,
                              margin = lmargin(l=legend_item_padding,
                                               r=legend_item_padding*1.5)), # Inherits from title
  legend.background = element_blank(), # No legend background (only box)
  legend.key = element_blank(), # No key background (only box)
  legend.box.background = element_rect(fill=legend_box_fill, colour = NA), # Inherits from rect
  legend.margin = lmargin(0),
  legend.box.margin = legend_box_margin,
  legend.spacing = lines(legend_series_spacing), # Spacing between different legends
  legend.spacing.x = lines(legend_series_spacing/2),
  legend.spacing.y = lines(legend_series_spacing/2),
  legend.key.size = legend_key_size,
  legend.key.height = NULL, # Inherits from legend.key.size
  legend.key.width = NULL, # Inherits from legend.key.size
  legend.box.spacing = legend_box_spacing, # Spacing between legend box and plot
  legend.text.align = NULL,
  legend.title.align = NULL,
  legend.position = "bottom",
  legend.direction = NULL, # Direction for arranging series within a legend
  legend.justification = "center",
  legend.box = NULL, # Direction for arranging multiple legends
  legend.box.just = NULL, # Justification of each legend within box
)

theme_axes <- theme( # Axis formatting
  axis.title = NULL, # Inherits from title
  axis.title.x = NULL, # Inherits from axis.title
  axis.title.y = NULL,
  axis.title.x.bottom = element_text(margin = lmargin(t=axis_title_spacing_x),
                                     vjust = 1),
  axis.title.x.top = element_text(margin = lmargin(b=axis_title_spacing_x),
                                  vjust = 0),
  axis.title.y.left = element_text(margin = lmargin(r=axis_title_spacing_y),
                                   vjust = 1, angle = 90),
  axis.title.y.right = element_text(margin = lmargin(l=axis_title_spacing_y),
                                    vjust = 0, angle = -90),
  axis.text = NULL, # Inherits from text
  axis.text.x = NULL, # Inherits from axis.text
  axis.text.x.bottom = element_text(margin = lmargin(t=axis_text_spacing),
                                    vjust = 1),
  axis.text.x.top = element_text(margin = lmargin(b=axis_text_spacing),
                                 vjust = 0),
  axis.text.y = NULL, # Inherits from axis.text
  axis.text.y.left = element_text(margin = lmargin(r=axis_text_spacing),
                                  hjust = 1),
  axis.text.y.right = element_text(margin = lmargin(l=axis_text_spacing),
                                   hjust = 0),
  axis.ticks = NULL, # Inherits from line
  axis.ticks.x = NULL,
  axis.ticks.x.top = NULL,
  axis.ticks.x.bottom = NULL,
  axis.ticks.y = NULL,
  axis.ticks.y.left = NULL,
  axis.ticks.y.right = NULL,
  axis.line = NULL, # Inherits from line
  axis.line.x = NULL,
  axis.line.x.top = NULL,
  axis.line.x.bottom = NULL,
  axis.line.y = NULL,
  axis.line.y.left = NULL,
  axis.line.y.right = NULL,
  axis.ticks.length = axis_tick_length,
  axis.ticks.length.x = NULL,
  axis.ticks.length.x.top = NULL,
  axis.ticks.length.x.bottom = NULL,
  axis.ticks.length.y = NULL,
  axis.ticks.length.y.left = NULL,
  axis.ticks.length.y.right = NULL,
)

theme_panel <- theme( # Panel formatting (except for axes)
  panel.grid = element_line(colour = "grey92"), # Inherits from line
  panel.grid.major = NULL,
  panel.grid.major.x = NULL,
  panel.grid.major.y = NULL,
  panel.grid.minor = element_line(linewidth = rel(0.5)),
  panel.grid.minor.x = NULL,
  panel.grid.minor.y = NULL,
  panel.background = NULL, # Inherits from rect
  panel.border = element_blank(), # No panel border
  panel.spacing = panel_spacing, # Spacing between plot panels
  panel.spacing.x = NULL,
  panel.spacing.y = NULL,
  panel.ontop = FALSE, # Put panel (gridlines etc.) on top of data?
)

theme_strips <- theme( # Facet strip formatting
  strip.text = element_text(size = fontsize_title, face = "bold",
                            margin = strip_margin), # Inherits from text
  strip.text.x = NULL,
  strip.text.y = element_text(angle = -90),
  strip.text.y.left = element_text(angle = 90),
  strip.background = element_blank(), # No strip background
  strip.background.x = NULL,
  strip.background.y = NULL,
  strip.switch.pad.grid = strip_panel_spacing,
  strip.switch.pad.wrap = strip_panel_spacing,
  strip.placement = "outside", # Placement of strips when axis on same side
)

theme_outer <- theme( # High-level plot formatting
  plot.title = element_text(hjust=0, face="plain", margin = lmargin(b=plot_title_spacing),
                            size = rel(1.2)), # Inherits from title
  plot.subtitle = element_text(margin = lmargin(b=plot_title_spacing)), # Inherits from title
  plot.caption = element_text(margin = lmargin(t=plot_title_spacing)), # Inherits from title
  plot.tag = element_text(face = "bold", size = fontsize_label,
                          margin = plot_tag_margin), # Inherits from title
  plot.background = element_blank(), # No plot background (pointless white square)
  plot.margin = plot_margin,
  plot.title.position = "panel", # Align title to panel or entire plot?
  plot.caption.position = "panel", # Align caption to panel or entire plot?
  plot.tag.position = "topleft", # Alignment of tag (e.g. subfigure letter)
  aspect.ratio = NULL # Aspect ratio of plot (height/width)
)

#------------------------------------------------------------------------------
# Complete themes
#------------------------------------------------------------------------------

# Main theme (for most single plots)
theme_base <- theme(complete = TRUE) + theme_core + theme_legend +
  theme_axes + theme_panel + theme_strips + theme_outer

# Alternative theme for tilted x-axis labels
theme_tilt <- theme_base + theme(
  axis.text.x = element_text(hjust = 1, angle = 45),
  axis.title.x = element_blank(),
)

# Alternative theme for internal legends
theme_internal_plain <- theme_base + theme(
  legend.box.background = element_rect(fill = legend_internal_plain_fill),
  legend.text = element_text(margin = lmargin(t=legend_internal_text_spacing,
                                              b=legend_internal_text_spacing,
                                              r=legend_item_padding)),
  legend.title = element_text(margin = legend_internal_title_margin),
  legend.justification = legend_internal_just_default,
  legend.position = legend_internal_position_default,
  legend.spacing.y = legend_internal_spacing_y,
  legend.box.margin = legend_box_margin * legend_internal_box_margin_scale,
)

# Another internal-legend theme, this time with an outlined box
theme_internal_strong <- theme_internal_plain + theme(
  legend.box.background = element_rect(fill = "white", colour = "black")
)

# Blank theme (e.g. for graphics from other libraries)
theme_blank <- theme_base + theme(
  axis.title = element_blank(),
  axis.title.x.bottom = NULL,
  axis.title.x.top = NULL,
  axis.title.y.left = NULL,
  axis.title.y.right = NULL,
  axis.text = element_blank(),
  axis.text.x.bottom = NULL,
  axis.text.x.top = NULL,
  axis.text.y.left = NULL,
  axis.text.y.right = NULL,
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  panel.grid = element_blank(),
  panel.grid.minor = NULL,
  legend.box.background = element_blank(),
  panel.background = element_blank(),
  plot.margin = lmargin(0),
  legend.margin = lmargin(0),
  legend.box.margin = lmargin(0),
  legend.position = "none",
  aspect.ratio = NULL)