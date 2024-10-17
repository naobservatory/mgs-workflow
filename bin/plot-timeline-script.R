library(pacman)
p_load(tidyverse, lubridate, plotly, htmlwidgets, rmarkdown, optparse)

convert_date_to_seconds <- function(date_str) {
  if (grepl("ms", date_str, ignore.case = TRUE)) {
    numeric_value <- as.numeric(gsub("[^0-9.]", "", date_str))
    return(numeric_value / 1000)

  } else if (grepl("^(\\d+h\\s*)?(\\d+m\\s*)?(\\d+(\\.\\d+)?s)?$", date_str)) {
    hours <- suppressWarnings(as.numeric(sub("^(\\d+)h.*$", "\\1", date_str, perl = TRUE)))
    minutes <- suppressWarnings(as.numeric(sub("^(?:\\d+h\\s*)?(\\d+)m.*$", "\\1", date_str, perl = TRUE)))
    seconds <- suppressWarnings(as.numeric(sub("^(?:\\d+h\\s*)?(?:\\d+m\\s*)?(\\d+(?:\\.\\d+)?)s$", "\\1", date_str, perl = TRUE)))

    hours <- ifelse(is.na(hours), 0, hours)
    minutes <- ifelse(is.na(minutes), 0, minutes)
    seconds <- ifelse(is.na(seconds), 0, seconds)
    return(hours * 3600 + minutes * 60 + seconds)
  } else {
    return(NULL)
  }
}

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input file path to trace file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output file path to html file", metavar="character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output)) {
  stop("Input and output file paths must be provided. Use -h for help.")
}

data <- read_tsv(opt$input)

selective_data <- data %>%
  select(task_id, name, realtime, submit, start, complete, duration, realtime, rss) %>%
  mutate(
    ending = as_datetime(start) + seconds(lapply(realtime, convert_date_to_seconds)),
    midpoint = start + (ending - start) / 2
  )

selective_data <- selective_data %>%
  arrange(desc(task_id)) %>%
  mutate(name = factor(name, levels = unique(name)))


g <- ggplot(selective_data) +
  geom_segment(aes(y = name, yend = name, x = submit, xend = start, color = "Time from submission to task starting"),
               size = 8) +
  geom_segment(aes(y = name, yend = name, x = start, xend = ending, color = "Time from task starting to task finishing"),
               size = 8) +
  geom_segment(aes(y = name, yend = name, x = ending, xend = complete, color = "Time spent staging and cleaning/unstaging input and output data"),
               size = 8) +
  geom_text(aes(x = midpoint, y = name, label = sprintf("%s/%s", duration, rss)), size = 3) +
  scale_color_manual(values = c("Time from submission to task starting" = "#fc8d59",
                                "Time from task starting to task finishing" = "#ffffbf",
                                "Time spent staging and cleaning/unstaging input and output data" = "#91bfdb")) +
  theme_light() +
  labs(x = "", y = "", color = "Task Status") +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black", face = "bold"),
    axis.title.x = element_text(size = 10, color = "black", face = "bold"),
    axis.title.y = element_text(size = 10, color = "black", face = "bold"),
    panel.grid.major.x = element_line(size = 1, color = "black", linetype = "dotted"),
    panel.grid.minor.x = element_line(color = "black", linetype = "dotted"),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

g

interactive_plot <- ggplotly(g, height = nrow(selective_data) * 25)
saveWidget(interactive_plot, opt$output)

cat("Plot saved to:", opt$output, "\n")
