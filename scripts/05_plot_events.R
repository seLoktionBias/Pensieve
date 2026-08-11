#!/usr/bin/env Rscript
# Pensieve event figures.
#
# Two figures, because one cannot do both jobs on a 100-tip tree:
#
#   <GENE>.pseudogenization_tree      the headline figure. Branches coloured by
#                                     ORF status: saturated red where the gene
#                                     broke, pale red for everything that
#                                     inherited it, grey for intact. Events
#                                     labelled with alignment start above the
#                                     branch and length below.
#
#   <GENE>.event_map                  every event as a bar at its real ALIGNMENT
#                                     POSITION, on rows shared with the tree.
#                                     Events on one branch separate along the
#                                     alignment axis instead of stacking, which
#                                     is what makes a 14-event branch readable.
#
# Only internodes that carry an event are labelled. Labelling all ~100 internal
# nodes, and stacking two text labels per event on a short branch, is what made
# earlier figures unreadable.
#
# Figure size, tip font and the tip-label strip are derived from the tip count
# and the longest name. Requires only ape + ggplot2.

suppressPackageStartupMessages({ library(ape); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  hit <- hit[1]
  if (hit >= length(args)) return(default)
  args[hit + 1]
}
is_yes <- function(x) tolower(as.character(x)) %in% c("yes", "true", "t", "1")
to_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes")

gene       <- get_arg("--gene")
tree_file  <- get_arg("--tree")
events_file<- get_arg("--events")
orf_file   <- get_arg("--orf-transitions", NA_character_)
outdir     <- get_arg("--outdir", ".")
dated      <- is_yes(get_arg("--dated", "yes"))
aln_len_arg<- get_arg("--alignment-length", NA_character_)
min_sites  <- as.numeric(get_arg("--min-sites", "1"))
show_tips  <- is_yes(get_arg("--show-tips", "yes"))
width_arg  <- get_arg("--width", "auto")
height_arg <- get_arg("--height", "auto")
tipfont_arg<- get_arg("--tip-font", "auto")

# A typo'd or renamed flag (e.g. "--transitions" instead of "--orf-transitions")
# does not match any get_arg() lookup above and would otherwise be silently
# ignored, leaving that setting at its default -- for --orf-transitions
# specifically, that default is "no ORF history at all", which renders every
# branch as plain grey "intact" with no error or warning anywhere. Reject any
# "--xxx"-shaped token that isn't one of the flags this script recognises so
# a bad flag name fails loudly here instead of silently producing a wrong,
# uniformly grey plot.
known_flags <- c("--gene", "--tree", "--events", "--orf-transitions", "--outdir",
                  "--dated", "--alignment-length", "--min-sites", "--show-tips",
                  "--width", "--height", "--tip-font")
flag_positions <- grep("^--", args)
unknown_flags <- setdiff(args[flag_positions], known_flags)
if (length(unknown_flags) > 0)
  stop("Unrecognised flag(s): ", paste(unknown_flags, collapse = ", "),
       ". Known flags: ", paste(known_flags, collapse = ", "))

if (is.na(gene))       stop("--gene is required (figures must be identifiable)")
if (is.na(tree_file) || is.na(events_file)) stop("--tree and --events are required")
for (f in c(tree_file, events_file)) if (!file.exists(f)) stop("File not found: ", f)
if (!is.na(orf_file) && !file.exists(orf_file))
  stop("--orf-transitions file not found: ", orf_file)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_tsv_base <- function(path) {
  if (is.na(path) || !file.exists(path) || file.info(path)$size == 0) return(data.frame())
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
             check.names = FALSE, quote = "", comment.char = "")
}

# ---------------------------------------------------------------- tree layout
tr <- read.tree(tree_file)
if (is.null(tr$node.label)) tr$node.label <- rep(NA_character_, tr$Nnode)
blank <- is.na(tr$node.label) | tr$node.label == ""
if (any(blank)) tr$node.label[blank] <- paste0("Node", which(blank))
if (!dated || is.null(tr$edge.length)) tr$edge.length <- rep(1, nrow(tr$edge))
ntip <- length(tr$tip.label)

tmpf <- tempfile(fileext = ".pdf"); pdf(tmpf, width = 8, height = 8)
plot(tr, show.tip.label = FALSE, no.margin = TRUE, direction = "rightwards")
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
dev.off(); unlink(tmpf)

coord <- data.frame(node = seq_along(pp$xx), x = pp$xx, y = pp$yy,
                    label = c(tr$tip.label, tr$node.label),
                    isTip = c(rep(TRUE, ntip), rep(FALSE, tr$Nnode)),
                    stringsAsFactors = FALSE)
edge_df <- data.frame(parent = tr$edge[, 1], child = tr$edge[, 2])
edge_df <- merge(edge_df, coord[, c("node", "x", "y", "label")], by.x = "parent", by.y = "node", all.x = TRUE)
names(edge_df)[names(edge_df) %in% c("x", "y", "label")] <- c("parent_x", "parent_y", "parent_label")
edge_df <- merge(edge_df, coord[, c("node", "x", "y", "label")], by.x = "child", by.y = "node", all.x = TRUE)
names(edge_df)[names(edge_df) %in% c("x", "y", "label")] <- c("child_x", "child_y", "child_label")
edge_df$branch <- paste0(edge_df$parent_label, "->", edge_df$child_label)
tree_w <- max(coord$x, na.rm = TRUE); if (!is.finite(tree_w) || tree_w <= 0) tree_w <- 1

# ---------------------------------------------------------------- events
ev_all <- read_tsv_base(events_file)
for (col in c("shared_event", "origin_is_tip", "reversal_below_origin", "ambiguous_origin"))
  if (col %in% names(ev_all)) ev_all[[col]] <- to_bool(ev_all[[col]])
aln_len <- if (!is.na(aln_len_arg)) as.numeric(aln_len_arg) else if (nrow(ev_all) > 0) max(ev_all$alignment_end) else 1
aln_len <- max(aln_len, 1)

ev <- ev_all[ev_all$event_length >= min_sites, , drop = FALSE]
ev <- merge(ev, edge_df[, c("branch", "parent_x", "child_x", "child_y")], by = "branch", all.x = TRUE)
lost <- ev[is.na(ev$child_y), , drop = FALSE]
if (nrow(lost) > 0) {
  warning(nrow(lost), " event(s) reference a branch absent from the tree: ",
          paste(head(unique(lost$branch), 5), collapse = ", "))
  ev <- ev[!is.na(ev$child_y), , drop = FALSE]
}

# ---------------------------------------------------------------- ORF status
orf <- read_tsv_base(orf_file)
edge_df$orf_class_detail <- "intact"
if (nrow(orf) > 0 && "orf_transition" %in% names(orf)) {
  m <- match(edge_df$branch, orf$branch)
  tr_class <- orf$orf_transition[m]
  edge_df$orf_class_detail <- ifelse(is.na(tr_class), "unknown", tr_class)
}
# Collapse the fine-grained orf_transition categories written to the TSVs
# (which keep full nuance, e.g. "confirmed disabling event but this may not
# be the very first loss because the root itself is never assigned a
# fabricated sequence") down to the four states a viewer actually needs to
# tell apart at a glance. Six different underlying categories used to all
# render in the same amber, including ones that are in fact CONFIRMED
# disabling events (just with an unresolved first-loss footnote) -- that
# made a confirmed lesion visually indistinguishable from genuine
# uncertainty. Only categories with no confident disabling signal either way
# stay amber now.
orf_display_map <- c(
  pseudogenization = "pseudogenization",
  confirmed_disabling_event_first_loss_unresolved = "pseudogenization",
  root_adjacent_disruption_first_loss_unresolved = "pseudogenization",
  sequence_disruption_first_loss_unresolved = "pseudogenization",
  already_pseudogenic = "already_pseudogenic",
  apparent_orf_restoration = "already_pseudogenic",
  ambiguous_disabling_event = "uncertain",
  sequence_state_uncertain = "uncertain",
  history_uncertain_no_confirmed_event = "uncertain",
  sequence_state_unavailable = "unavailable",
  intact = "intact",
  unknown = "unknown"
)
orf_colours <- c(
  pseudogenization = "#c0392b",
  already_pseudogenic = "#e8a6a0",
  uncertain = "#d68910",
  unavailable = "grey80",
  intact = "grey55",
  unknown = "grey80"
)
edge_df$orf_class <- ifelse(edge_df$orf_class_detail %in% names(orf_display_map),
                            orf_display_map[edge_df$orf_class_detail], "unknown")

# ---------------------------------------------------------------- sizing
max_chars <- if (show_tips) max(nchar(gsub("_", " ", tr$tip.label))) else 0
tip_font <- if (tipfont_arg == "auto") max(1.2, min(3.2, 320 / max(ntip, 1))) else as.numeric(tipfont_arg)
fig_h <- if (height_arg == "auto") max(7, ntip * 0.155 + 2.5) else as.numeric(height_arg)
fig_w <- if (width_arg == "auto") max(15, 11 + max_chars * 0.075) else as.numeric(width_arg)
gap <- tree_w * 0.03
usable_in <- max(fig_w - 2.2, 2)
# ggplot2 text `size` is in MILLIMETRES; .pt = 72.27/25.4. Treating it as points
# under-reserves the tip-label strip ~2.8x and the labels run into the map.
label_in <- if (show_tips) max_chars * (tip_font * 2.845) * 0.60 / 72 else 0
denom <- usable_in - label_in
label_w <- if (denom > 0.2) label_in * (2 * tree_w + 2 * gap) / denom else tree_w * 0.8

# Event markers are classified and coloured ONCE here, shared by both
# figures (Figure 1's per-branch tick marks and Figure 2's alignment-position
# bars) so the exact same event is never shown in two different colours
# depending on which figure it happens to be looked at in. Classification
# uses the same character_class/biological_interpretation/length_mod_3
# fields 04_ancestral_orf_walk.py itself uses to decide a branch's own
# pseudogenization call, plus origin_is_tip for "happened on an ancestral
# node" (shared, non-disrupting indels only -- an ambiguous or disabling
# event keeps its own category regardless of where it originates).
mc_len3 <- suppressWarnings(as.numeric(ev$length_mod_3))
is_stop_gain          <- ev$character_class == "stop_mask" & ev$biological_interpretation == "premature_stop_gained"
is_stop_lost           <- ev$character_class == "stop_mask" & ev$biological_interpretation == "premature_stop_lost"
is_frameshift_indel   <- ev$character_class == "indel" & !is.na(mc_len3) & mc_len3 != 0
is_inframe_indel      <- ev$character_class == "indel" & !is.na(mc_len3) & mc_len3 == 0
is_insertion          <- ev$biological_interpretation == "insertion_or_restoration"
is_deletion            <- ev$biological_interpretation == "deletion"
is_ancestral_origin    <- !ev$origin_is_tip
# "Ambiguous" here means the origin/direction of a potentially disabling
# change (a STOP allele, or a frameshifting-length indel) is unresolved --
# matching exactly which kinds of ambiguity 04_ancestral_orf_walk.py itself
# treats as ORF-relevant. An ambiguous IN-FRAME indel's polarity tie cannot
# flip intact/disrupted either way, so it is left in the ordinary in-frame
# categories below rather than flagged as if it might be disabling.
is_ambiguous_disabling <- (ev$character_class == "stop_mask" & ev$biological_interpretation == "ambiguous_stop_change") |
                          (is_frameshift_indel & ev$biological_interpretation == "ambiguous_indel_change")

ev$marker_class <- "other"
ev$marker_class[is_inframe_indel & is_insertion & is_ancestral_origin] <- "shared_inframe_insertion"
ev$marker_class[is_inframe_indel & is_deletion & is_ancestral_origin]  <- "shared_inframe_deletion"
ev$marker_class[is_ambiguous_disabling] <- "ambiguous_disabling_candidate"
ev$marker_class[is_stop_lost] <- "stop_lost_reversion"
ev$marker_class[is_frameshift_indel & !is_ambiguous_disabling] <- "pseudogenizing_indel"
ev$marker_class[is_stop_gain] <- "pseudogenizing_substitution"

# Deliberately NOT red: pseudogenizing markers can sit directly on top of
# branches that are already red (bold for "confirmed disabling", pale for
# "inherited"), so a same-family red marker on a red branch is nearly
# invisible. Both pseudogenizing marker colours are chosen to read clearly
# against both red branch shades. "Other" (ordinary, non-notable in-frame
# indels -- most events on most genes) is a pale, receding grey so it does
# not read as important; it deliberately does not try to stand out.
marker_colours <- c(
  pseudogenizing_substitution = "#17202a",
  pseudogenizing_indel = "#f1c40f",
  ambiguous_disabling_candidate = "#d68910",
  stop_lost_reversion = "#16a085",
  shared_inframe_insertion = "#2471a3",
  shared_inframe_deletion = "#bdb76b",
  other = "grey75"
)
marker_labels <- c(
  pseudogenizing_substitution = "Pseudogenizing substitution (premature stop)",
  pseudogenizing_indel = "Pseudogenizing indel (frameshift)",
  ambiguous_disabling_candidate = "Ambiguous disabling candidate (uncertain origin/direction)",
  stop_lost_reversion = "Premature stop lost (reversion)",
  shared_inframe_insertion = "Shared in-frame insertion (ancestral, non-disrupting)",
  shared_inframe_deletion = "Shared in-frame deletion (ancestral, non-disrupting)",
  other = "Other event (e.g. lineage-specific in-frame indel)"
)

subtitle <- paste0(nrow(ev), " events on ", length(unique(ev$branch)), " branches; ",
                   sum(ev$shared_event), " shared; ", ntip, " tips; ",
                   aln_len, " alignment columns",
                   if (any(ev$ambiguous_origin)) paste0("; ", sum(ev$ambiguous_origin),
                                                        " with an exact parsimony tie") else "")

base_tree <- function(p) {
  p + geom_segment(data = edge_df, aes(x = parent_x, xend = parent_x, y = parent_y, yend = child_y),
                   linewidth = 0.35, colour = "grey55") +
    # Only the horizontal segment is the evolutionary branch; the vertical one is
    # a drawing connector shared between siblings and must never take a child's colour.
    geom_segment(data = edge_df, aes(x = parent_x, xend = child_x, y = child_y, yend = child_y,
                                     colour = orf_class), linewidth = 1.0)
}

tip_df <- coord[coord$isTip, , drop = FALSE]
tip_df$plot_label <- gsub("_", " ", tip_df$label)
broken_tips <- if (nrow(orf) > 0) orf$child_node[
  orf$known_pseudogenic_history %in% c(TRUE, "True", "TRUE", "true", "1")
] else character(0)
tip_df$broken <- tip_df$label %in% broken_tips

time_breaks <- pretty(c(0, tree_w), n = 6); time_breaks <- time_breaks[time_breaks >= 0 & time_breaks <= tree_w]

# ================================================================ FIGURE 1
ev1 <- ev[order(ev$branch, ev$alignment_start), ]
ev1$rank <- ave(seq_len(nrow(ev1)), ev1$branch, FUN = seq_along)
ev1$non  <- ave(seq_len(nrow(ev1)), ev1$branch, FUN = length)
ev1$ex   <- ev1$parent_x + (ev1$child_x - ev1$parent_x) * (ev1$rank / (ev1$non + 1))
ev1$pos_label <- as.character(ev1$alignment_start)
ev1$len_label <- as.character(ev1$event_length)

# marker_class/marker_colours/marker_labels were already computed once above
# (shared with Figure 2) so the exact same event is never a different colour
# depending on which figure it is looked at in.
bar_half_w <- tree_w * 0.0008
bar_half_h <- 0.16
ev1$mxmin <- ev1$ex - bar_half_w; ev1$mxmax <- ev1$ex + bar_half_w
ev1$mymin <- ev1$child_y - bar_half_h; ev1$mymax <- ev1$child_y + bar_half_h

event_nodes <- edge_df[edge_df$child_label %in% unique(ev$origin_node[!ev$origin_is_tip]), , drop = FALSE]
if (nrow(event_nodes) > 0) {
  event_nodes$lx <- (event_nodes$parent_x + event_nodes$child_x) / 2
  event_nodes$ly <- event_nodes$child_y + 0.95
  event_nodes$label <- event_nodes$child_label
}

p1 <- base_tree(ggplot()) +
  geom_rect(data = ev1, aes(xmin = mxmin, xmax = mxmax, ymin = mymin, ymax = mymax, fill = marker_class),
            colour = NA) +
  # Fixed offset from the branch's own row for every event, regardless of how
  # many events share that branch (no longer scaled by `rank`): events on one
  # branch are already kept apart along the x-axis via `ex` above, so a
  # constant, small y-offset keeps every label on its own branch's axis
  # instead of fanning upward/downward into neighbouring branches' rows.
  # vjust anchors each label at the edge closest to the bar and grows it
  # AWAY from the bar (up for the position label, down for the length
  # label), on top of a gap already bigger than the bar's own half-height,
  # so the text can never sit on top of/inside the coloured marker.
  geom_text(data = ev1, aes(x = ex, y = child_y + 0.22, label = pos_label),
            colour = "grey10", size = 1.5, fontface = "bold", vjust = 0) +
  geom_text(data = ev1, aes(x = ex, y = child_y - 0.22, label = len_label),
            colour = "grey35", size = 1.3, vjust = 1)
if (nrow(event_nodes) > 0)
  p1 <- p1 + geom_label(data = event_nodes, aes(x = lx, y = ly, label = label),
                        size = max(1.4, tip_font * 0.8), colour = "#7b241c", fill = "white",
                        label.size = 0.15, label.padding = grid::unit(0.08, "lines"))
if (show_tips)
  p1 <- p1 +
    geom_text(data = tip_df[!tip_df$broken, , drop = FALSE], aes(x = tree_w + gap, y = y, label = plot_label),
              hjust = 0, size = tip_font, colour = "grey25", fontface = "italic") +
    geom_text(data = tip_df[tip_df$broken, , drop = FALSE], aes(x = tree_w + gap, y = y, label = plot_label),
              hjust = 0, size = tip_font, colour = "#c0392b", fontface = "bold.italic")

p1 <- p1 +
  scale_colour_manual(name = "Pseudogenization history", values = orf_colours,
                      breaks = names(orf_colours)[names(orf_colours) %in% edge_df$orf_class],
                      labels = c(
                        pseudogenization = "Confirmed disabling event",
                        already_pseudogenic = "Inherited pseudogenic history",
                        uncertain = "Uncertain; no confirmed event",
                        unavailable = "Sequence unavailable",
                        intact = "No inferred disabling history",
                        unknown = "No ancestral sequence")[
                        names(orf_colours)[names(orf_colours) %in% edge_df$orf_class]]) +
  scale_fill_manual(name = "Event marker", values = marker_colours,
                    breaks = names(marker_colours)[names(marker_colours) %in% ev1$marker_class],
                    labels = marker_labels[names(marker_colours)[names(marker_colours) %in% ev1$marker_class]]) +
  guides(fill = guide_legend(order = 1, nrow = 2, byrow = TRUE, title.position = "top"),
        colour = guide_legend(order = 2, nrow = 2, byrow = TRUE, title.position = "top")) +
  scale_x_continuous(name = if (dated) "Time before present" else "Node depth",
                     breaks = time_breaks,
                     labels = if (dated) round(tree_w - time_breaks, 1) else round(time_breaks, 1),
                     limits = c(0, tree_w + gap + label_w * 1.05),
                     expand = expansion(mult = c(0.005, 0.005))) +
  scale_y_continuous(name = NULL, expand = expansion(mult = c(0.01, 0.02))) +
  coord_cartesian(clip = "off") +
  ggtitle(paste0(gene, ": reconstructed indel and pseudogenizing events"),
          subtitle = paste0(subtitle,
                            "\nEvent labels: alignment start position above the branch, event length below.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 7.5, colour = "grey30"),
        legend.position = "bottom", legend.box = "vertical",
        legend.key.size = unit(0.32, "cm"), legend.text = element_text(size = 6.8),
        legend.title = element_text(size = 7.5, face = "bold"),
        legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(0.05, "cm"),
        plot.margin = margin(8, 12, 8, 8))

f1 <- file.path(outdir, paste0(gene, ".pseudogenization_tree"))
ggsave(paste0(f1, ".pdf"), p1, width = fig_w * 0.8, height = fig_h, units = "in", limitsize = FALSE)
ggsave(paste0(f1, ".png"), p1, width = fig_w * 0.8, height = fig_h, units = "in", dpi = 200, limitsize = FALSE)
message("Wrote ", f1, ".pdf / .png")

# ================================================================ FIGURE 2
map_x0 <- tree_w + gap + label_w + gap; map_w <- tree_w; map_x1 <- map_x0 + map_w
site_to_x <- function(s) map_x0 + (s - 1) / max(aln_len - 1, 1) * map_w
ev2 <- ev
ev2$xmin <- site_to_x(ev2$alignment_start); ev2$xmax <- site_to_x(ev2$alignment_end + 1)
min_bar <- map_w * 0.0035
short <- (ev2$xmax - ev2$xmin) < min_bar; ev2$xmax[short] <- ev2$xmin[short] + min_bar

# Greedy sub-lane packing: two events on one branch collide either genuinely
# (a 1 bp mask inside an insertion) or because a 1 bp bar was widened to stay
# visible. Overlapping bars go to separate lanes instead of on top of each other.
pad <- map_w * 0.002; ev2$lane <- 0L; ev2$n_lanes <- 1L
for (yy in unique(ev2$child_y)) {
  idx <- which(ev2$child_y == yy); idx <- idx[order(ev2$xmin[idx])]; lane_end <- numeric(0)
  for (i in idx) {
    slot <- which(lane_end + pad <= ev2$xmin[i])
    if (length(slot) > 0) { k <- slot[1]; ev2$lane[i] <- k - 1L; lane_end[k] <- ev2$xmax[i] }
    else { lane_end <- c(lane_end, ev2$xmax[i]); ev2$lane[i] <- length(lane_end) - 1L }
  }
  ev2$n_lanes[idx] <- length(lane_end)
}
row_budget <- 0.86
ev2$lane_h <- row_budget / ev2$n_lanes
ev2$ymin <- ev2$child_y - row_budget / 2 + ev2$lane * ev2$lane_h
ev2$ymax <- ev2$ymin + ev2$lane_h * 0.9

site_breaks <- pretty(c(1, aln_len), n = 8); site_breaks <- site_breaks[site_breaks >= 1 & site_breaks <= aln_len]
p2 <- ggplot() +
  annotate("rect", xmin = map_x0, xmax = map_x1, ymin = 0.3, ymax = ntip + 0.7, fill = "grey96", colour = NA) +
  annotate("segment", x = site_to_x(site_breaks), xend = site_to_x(site_breaks),
           y = 0.3, yend = ntip + 0.7, colour = "white", linewidth = 0.4) +
  geom_segment(data = data.frame(y = unique(ev2$child_y)),
               aes(x = tree_w + gap + label_w * 0.99, xend = map_x1, y = y, yend = y),
               colour = "grey85", linewidth = 0.25)
p2 <- base_tree(p2) +
  geom_rect(data = ev2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                            fill = marker_class), colour = NA) +
  geom_rect(data = ev2[ev2$shared_event, , drop = FALSE],
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, colour = "black", linewidth = 0.3) +
  scale_fill_manual(name = "Event marker", values = marker_colours,
                    breaks = names(marker_colours)[names(marker_colours) %in% ev2$marker_class],
                    labels = marker_labels[names(marker_colours)[names(marker_colours) %in% ev2$marker_class]]) +
  scale_colour_manual(name = "Pseudogenization history", values = orf_colours,
                      breaks = names(orf_colours)[names(orf_colours) %in% edge_df$orf_class], guide = "none") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, title.position = "top"))
# Every internal node, not just event-carrying ones: on a 100+ tip tree the
# unlabelled backbone made it impossible to tell which branch in the event
# map corresponds to which node in alignment_events.tsv/
# orf_transitions_by_branch.tsv, especially for deep nodes far from any tip.
# Figure 1 deliberately labels only event-carrying nodes (it already has
# event markers/labels competing for the same vertical space per branch);
# Figure 2 has no such competition in the tree panel, so every node gets one.
internal_df <- coord[!coord$isTip, , drop = FALSE]
p2 <- p2 +
  geom_label(data = internal_df, aes(x = x, y = y, label = label),
            size = max(1.4, tip_font * 0.7), colour = "grey20", fill = "white",
            label.size = 0.1, label.padding = grid::unit(0.08, "lines"), alpha = 0.9)
if (show_tips)
  p2 <- p2 +
    geom_text(data = tip_df[!tip_df$broken, , drop = FALSE], aes(x = tree_w + gap, y = y, label = plot_label),
              hjust = 0, size = tip_font, colour = "grey25", fontface = "italic") +
    geom_text(data = tip_df[tip_df$broken, , drop = FALSE], aes(x = tree_w + gap, y = y, label = plot_label),
              hjust = 0, size = tip_font, colour = "#c0392b", fontface = "bold.italic")

p2 <- p2 +
  annotate("text", x = tree_w / 2, y = ntip + 2.4, label = if (dated) "Dated tree" else "Cladogram",
           fontface = "bold", size = 3.2, colour = "grey20") +
  annotate("text", x = (map_x0 + map_x1) / 2, y = ntip + 2.4, label = "Event position in alignment",
           fontface = "bold", size = 3.2, colour = "grey20") +
  scale_x_continuous(name = paste(if (dated) "Time before present" else "Node depth",
                                  "           |           Alignment position (bp)"),
                     breaks = c(time_breaks, site_to_x(site_breaks)),
                     labels = c(if (dated) round(tree_w - time_breaks, 1) else round(time_breaks, 1), site_breaks),
                     limits = c(0, map_x1 + map_w * 0.01), expand = expansion(mult = c(0.005, 0.005))) +
  scale_y_continuous(name = NULL, limits = c(0, ntip + 3.4), expand = expansion(mult = c(0.004, 0.004))) +
  coord_cartesian(clip = "off") +
  ggtitle(paste0(gene, ": event map"), subtitle = subtitle) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 8, colour = "grey30"),
        legend.position = "bottom",
        legend.key.size = unit(0.32, "cm"), legend.text = element_text(size = 6.8),
        legend.title = element_text(size = 7.5, face = "bold"),
        legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(0.05, "cm"),
        plot.margin = margin(8, 12, 8, 8))

f2 <- file.path(outdir, paste0(gene, ".event_map"))
ggsave(paste0(f2, ".pdf"), p2, width = fig_w, height = fig_h, units = "in", limitsize = FALSE)
ggsave(paste0(f2, ".png"), p2, width = fig_w, height = fig_h, units = "in", dpi = 200, limitsize = FALSE)
message("Wrote ", f2, ".pdf / .png")

write.table(ev[order(ev$alignment_start), ],
            file.path(outdir, paste0(gene, ".plotted_events.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)
message(gene, ": ", nrow(ev), " events plotted, ", sum(ev$shared_event), " shared.")
