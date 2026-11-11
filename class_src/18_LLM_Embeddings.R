## ============================================================================
##  MS_LJMR :: 18_LLM_Embeddings.R — Large Language Models & Embeddings
##
##  Author: Luis J. Madrigal-Roca  (new session for the MS_LJMR edition)
##
##  NEW session (revised Fall-2026 schedule).
##  Readings: Vaswani et al. 2017 ("Attention is all you need");
##            Levy & Goldberg 2014 (Neural word embedding as implicit matrix
##            factorization); Rives et al. 2021 (Biological structure and
##            function emerge from scaling unsupervised learning to 250M
##            protein sequences).
##
##  Bridge to the course: an *embedding* is the exact multivariate object we
##  have studied all semester — a high-dimensional vector per item in which
##  geometry encodes meaning (cosine similarity ~ semantic/functional
##  similarity). Every tool from this book applies directly: distances, SVD,
##  PCA, UMAP, clustering, and nearest-centroid classification.
##
##  Four self-contained examples (100% offline, plotted with ggplot2 + ggrepel):
##  ---------------------------------------------------------------------------
##  PART A. Geometry of Meaning & Vector Analogies (32 concepts across 5 guilds)
##          king - man + woman ~ queen; shark - ocean + savanna ~ lion.
##  PART B. Learning Word Embeddings from Scratch on a Toy Ecological Corpus
##          Co-occurrence window -> Positive PMI (PPMI) -> truncated SVD -> UMAP.
##  PART C. Biological Sequence Embeddings ("Protein Language Model" in miniature)
##          Toy dataset of 30 functional peptides across 5 families + zero-shot
##          functional classification of 3 mystery sequences.
##  PART D. Inside the Transformer: Scaled Dot-Product Self-Attention
##          How Vaswani et al. (2017) disambiguates polysemous words ("vector",
##          "expression", "culture") by pulling them toward their context.
## ============================================================================

##  PLAIN CLASSROOM EDITION
##  Generated from R/18_LLM_Embeddings.R by scripts/make_class_src.py --
##  edit R/18_LLM_Embeddings.R and regenerate; changes made here will be overwritten.
## ============================================================================

# R packages required
#install.packages("ggplot2")
#install.packages("ggrepel")
#install.packages("uwot")
library(ggplot2)
library(ggrepel)
library(uwot)

set.seed(42)

## Shared helper: L2-normalise rows of a matrix so dot product == cosine similarity
l2_normalize <- function(M) {
  norms <- sqrt(rowSums(M^2))
  norms[norms == 0] <- 1
  M / norms
}

cosine_sim <- function(a, b) sum(a * b)

## Consistent publication palette across all four parts
guild_palette <- c(
  "Royalty & Family"  = "#1f3b73",
  "World Geography"   = "#8c2d3a",
  "Marine Ecology"    = "#1b7837",
  "Savanna Ecology"   = "#c88a2a",
  "Molecular Biology" = "#6a4c93"
)

## ############################################################################
##  PART A — GEOMETRY OF MEANING & COMPOSITIONAL VECTOR ANALOGIES
## ############################################################################
##  We construct a 32-concept embedding space in d = 32 dimensions using
##  orthogonal latent axes (domain + shared relational modifiers: female/male,
##  capital/country, juvenile/adult, apex-predator/herbivore, template/product).
##  Because modifiers are shared across domains, linear vector arithmetic
##  solves analogies across geography, royalty, ecology, and molecular biology.

d <- 32
## Create orthonormal basis vectors via QR decomposition of a random Gaussian matrix
Q_basis <- qr.Q(qr(matrix(rnorm(d * 18), nrow = d, ncol = 18)))
ax <- function(i) Q_basis[, i]

## Domain axes (1..5) and shared relational feature axes (6..14)
ax_royal   <- ax(1);  ax_geo     <- ax(2);  ax_marine  <- ax(3)
ax_savanna <- ax(4);  ax_molbio  <- ax(5)
ax_female  <- ax(6);  ax_capital <- ax(7);  ax_juvenile <- ax(8)
ax_apex    <- ax(9);  ax_herbiv  <- ax(10); ax_primary  <- ax(11)
ax_france  <- ax(12); ax_italy   <- ax(13); ax_japan    <- ax(14)
ax_brazil  <- ax(15); ax_rnatag  <- ax(16); ax_prottag  <- ax(17)

noise <- function(sd = 0.08) rnorm(d, 0, sd)

E_raw <- rbind(
  ## 1. Royalty & Family
  king        = 1.2 * ax_royal - 0.8 * ax_female + 0.5 * ax_apex     + noise(),
  queen       = 1.2 * ax_royal + 0.8 * ax_female + 0.5 * ax_apex     + noise(),
  prince      = 1.1 * ax_royal - 0.8 * ax_female + 0.7 * ax_juvenile + noise(),
  princess    = 1.1 * ax_royal + 0.8 * ax_female + 0.7 * ax_juvenile + noise(),
  man         = 0.6 * ax_royal - 0.9 * ax_female                     + noise(),
  woman       = 0.6 * ax_royal + 0.9 * ax_female                     + noise(),
  boy         = 0.6 * ax_royal - 0.9 * ax_female + 0.7 * ax_juvenile + noise(),
  girl        = 0.6 * ax_royal + 0.9 * ax_female + 0.7 * ax_juvenile + noise(),

  ## 2. World Geography (Country <-> Capital)
  france      = 1.2 * ax_geo + 0.8 * ax_france - 0.7 * ax_capital + noise(),
  paris       = 1.2 * ax_geo + 0.8 * ax_france + 0.7 * ax_capital + noise(),
  italy       = 1.2 * ax_geo + 0.8 * ax_italy  - 0.7 * ax_capital + noise(),
  rome        = 1.2 * ax_geo + 0.8 * ax_italy  + 0.7 * ax_capital + noise(),
  japan       = 1.2 * ax_geo + 0.8 * ax_japan  - 0.7 * ax_capital + noise(),
  tokyo       = 1.2 * ax_geo + 0.8 * ax_japan  + 0.7 * ax_capital + noise(),
  brazil      = 1.2 * ax_geo + 0.8 * ax_brazil - 0.7 * ax_capital + noise(),
  brasilia    = 1.2 * ax_geo + 0.8 * ax_brazil + 0.7 * ax_capital + noise(),

  ## 3. Marine Ecology (Habitat / Producer / Herbivore / Apex predator)
  ocean       = 1.3 * ax_marine                                    + noise(),
  kelp        = 1.2 * ax_marine + 0.8 * ax_primary                 + noise(),
  manatee     = 1.2 * ax_marine + 0.8 * ax_herbiv                  + noise(),
  parrotfish  = 1.1 * ax_marine + 0.7 * ax_herbiv                  + noise(),
  shark       = 1.2 * ax_marine + 0.9 * ax_apex                    + noise(),
  orca        = 1.2 * ax_marine + 0.85 * ax_apex                   + noise(),

  ## 4. Savanna Ecology (mirrors Marine trophic roles!)
  savanna     = 1.3 * ax_savanna                                   + noise(),
  acacia      = 1.2 * ax_savanna + 0.8 * ax_primary                + noise(),
  zebra       = 1.2 * ax_savanna + 0.8 * ax_herbiv                 + noise(),
  wildebeest  = 1.1 * ax_savanna + 0.7 * ax_herbiv                 + noise(),
  lion        = 1.2 * ax_savanna + 0.9 * ax_apex                   + noise(),
  cheetah     = 1.2 * ax_savanna + 0.85 * ax_apex                  + noise(),

  ## 5. Molecular Biology (Central Dogma)
  dna         = 1.3 * ax_molbio - 0.6 * ax_rnatag                  + noise(),
  mrna        = 1.3 * ax_molbio + 0.7 * ax_rnatag                  + noise(),
  protein     = 1.3 * ax_molbio + 0.8 * ax_prottag                 + noise(),
  polymerase  = 1.2 * ax_molbio - 0.6 * ax_rnatag + 0.6 * ax_apex  + noise(),
  ribosome    = 1.2 * ax_molbio + 0.7 * ax_rnatag + 0.6 * ax_apex  + noise(),
  codon       = 1.1 * ax_molbio + 0.5 * ax_rnatag                  + noise()
)

E <- l2_normalize(E_raw)

concept_guild <- c(
  rep("Royalty & Family",  8),
  rep("World Geography",   8),
  rep("Marine Ecology",    6),
  rep("Savanna Ecology",   6),
  rep("Molecular Biology", 6)
)
names(concept_guild) <- rownames(E)

## ---- 1. Solve analogies via vector arithmetic: a : b :: c : ? --------------
solve_analogy <- function(a, b, c_word, emb = E, top_k = 3) {
  ## Target vector: b - a + c (e.g. king - man + woman)
  q <- emb[b, ] - emb[a, ] + emb[c_word, ]
  q <- q / sqrt(sum(q^2))
  sims <- as.numeric(emb %*% q)
  names(sims) <- rownames(emb)
  sims <- sims[!names(sims) %in% c(a, b, c_word)]
  sort(sims, decreasing = TRUE)[seq_len(top_k)]
}

cat("=== PART A: Vector Analogies in 32-D Embedding Space ===\n")
analogies <- list(
  c("man",     "king",   "woman"),     # -> queen
  c("france",  "paris",  "japan"),     # -> tokyo
  c("ocean",   "shark",  "savanna"),   # -> lion (marine apex -> savanna apex!)
  c("acacia",  "zebra",  "kelp"),      # -> manatee / parrotfish (savanna herbivore -> marine herbivore!)
  c("dna",     "polymerase", "mrna")   # -> ribosome (template:machinery analogy!)
)
for (an in analogies) {
  res <- solve_analogy(an[1], an[2], an[3])
  cat(sprintf("  %-10s : %-10s :: %-10s : %-10s (cos = %.3f; next: %s = %.3f)\n",
              an[1], an[2], an[3], names(res)[1], res[1], names(res)[2], res[2]))
}

## ---- Figure 1: PCA Biplot of Embeddings + Analogy Parallelograms -----------
pc_A <- prcomp(E, center = TRUE)
ve_A <- 100 * (pc_A$sdev^2 / sum(pc_A$sdev^2))

df_A <- data.frame(
  word  = rownames(E),
  guild = factor(concept_guild, levels = names(guild_palette)),
  PC1   = pc_A$x[, 1],
  PC2   = pc_A$x[, 2],
  row.names = NULL
)

## Build segment table for analogy arrows (from -> to)
arrow_pairs <- rbind(
  data.frame(from = "man",     to = "king",       relation = "Gender / Royalty"),
  data.frame(from = "woman",   to = "queen",      relation = "Gender / Royalty"),
  data.frame(from = "france",  to = "paris",      relation = "Country -> Capital"),
  data.frame(from = "italy",   to = "rome",       relation = "Country -> Capital"),
  data.frame(from = "japan",   to = "tokyo",      relation = "Country -> Capital"),
  data.frame(from = "manatee", to = "shark",      relation = "Herbivore -> Apex"),
  data.frame(from = "zebra",   to = "lion",       relation = "Herbivore -> Apex"),
  data.frame(from = "dna",     to = "polymerase", relation = "Template -> Enzyme"),
  data.frame(from = "mrna",    to = "ribosome",   relation = "Template -> Enzyme")
)
arrow_df <- merge(arrow_pairs, df_A[, c("word", "PC1", "PC2")],
                  by.x = "from", by.y = "word")
arrow_df <- merge(arrow_df, df_A[, c("word", "PC1", "PC2")],
                  by.x = "to", by.y = "word", suffixes = c("_start", "_end"))

par(mfrow = c(1, 1))
p1 <- ggplot(df_A, aes(x = PC1, y = PC2, color = guild)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  geom_segment(
    data = arrow_df,
    aes(x = PC1_start, y = PC2_start, xend = PC1_end, yend = PC2_end),
    inherit.aes = FALSE,
    color = "grey35", linewidth = 0.65, linetype = "solid",
    arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed")
  ) +
  geom_point(size = 3.2, alpha = 0.9) +
  ggrepel::geom_text_repel(
    aes(label = word),
    size = 3.5, fontface = "bold", box.padding = 0.35,
    point.padding = 0.25, max.overlaps = Inf, seed = 42,
    show.legend = FALSE
  ) +
  scale_color_manual(values = guild_palette, name = "Semantic Guild") +
  labs(
    title    = "Geometry of Meaning: Word Embeddings & Compositional Analogies",
    subtitle = "Parallel arrows encode shared semantic relations (Country -> Capital, Herbivore -> Apex, Template -> Machinery)",
    x        = sprintf("PC1 (%.1f%% variance)", ve_A[1]),
    y        = sprintf("PC2 (%.1f%% variance)", ve_A[2])
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(face = "bold")
  )
print(p1)


## ############################################################################
##  PART B — LEARNING EMBEDDINGS FROM RAW TEXT (CO-OCCURRENCE -> PPMI -> SVD)
## ############################################################################
##  "You shall know a word by the company it keeps" (Firth 1957).
##  Levy & Goldberg (2014) proved that Skip-gram Word2Vec is implicitly
##  factorising the Positive Pointwise Mutual Information (PPMI) word-context
##  co-occurrence matrix. Here we train embeddings from scratch on a miniature
##  field-biology corpus of 48 sentences across 4 ecosystems/disciplines.

corpus <- c(
  ## Biome 1: Coral Reef Ecology
  "coral polyp calcifies aragonite reef skeleton with symbiotic zooxanthellae algae",
  "parrotfish grazes turf algae on coral reef crest near lagoon",
  "clownfish shelters inside anemone tentacles on shallow coral reef",
  "thermal stress triggers coral bleaching by expelling symbiotic zooxanthellae algae",
  "reef shark patrols coral dropoff hunting parrotfish and damselfish in lagoon",
  "symbiotic zooxanthellae transfer photosynthate to coral polyp in clear lagoon",
  "ocean acidification slows aragonite calcification of coral reef skeleton",
  "damselfish defends turf algae garden on branching coral colony",
  "anemone and clownfish form mutualistic symbiosis on tropical reef",
  "bleaching mortality shifts coral reef toward turf algae dominance",
  "lagoon parrotfish bioerodes calcified coral skeleton while grazing algae",
  "reef shark and moray eel hunt nocturnal damselfish around coral heads",

  ## Biome 2: Sonoran Desert Ecology
  "saguaro cactus stores water in pleated succulent stem during arid drought",
  "nocturnal kangaroo_rat forages desert seeds and shelters in cool burrow",
  "horned_lizard basks on arid canyon rock and ambushes desert ants",
  "desert coyote hunts nocturnal kangaroo_rat across arid creosote flats",
  "succulent agave and saguaro cactus open stomata at night via cam_photosynthesis",
  "arid drought forces kangaroo_rat to concentrate metabolic water in burrow",
  "creosote bush and saguaro cactus dominate arid desert bajada",
  "canyon rattlesnake ambushes nocturnal kangaroo_rat near desert burrow",
  "cam_photosynthesis conserves water in succulent cactus and desert agave",
  "horned_lizard and rattlesnake shelter under canyon ledges during arid heat",
  "coyotedisperses saguaro cactus fruit seeds across desert canyon",
  "nocturnal burrow microclimate buffers kangaroo_rat from arid desert extremes",

  ## Biome 3: Boreal Taiga Forest
  "boreal spruce and larch needles tolerate subzero permafrost winter snowpack",
  "canada_lynx stalks snowshoe_hare through dense boreal spruce forest",
  "ectomycorrhizal fungi colonize spruce roots in acidic permafrost soil",
  "snowshoe_hare browses willow twigs above deep winter snowpack in taiga",
  "boreal caribou migrates across permafrost tundra grazing arboreal lichen",
  "wolf_pack hunts caribou and snowshoe_hare across frozen taiga snowpack",
  "acidic podzol soil and permafrost slow decomposition under boreal spruce",
  "canada_lynx population cycles track snowshoe_hare abundance in taiga forest",
  "ectomycorrhizal fungi exchange nitrogen with boreal larch and spruce roots",
  "caribou digs through winter snowpack to reach tundra lichen mat",
  "subzero winter freezes permafrost beneath boreal spruce and willow taiga",
  "wolf_pack and canada_lynx partition prey across snowy boreal forest",

  ## Discipline 4: Molecular Genetics & Gene Regulation
  "rna_polymerase binds promoter tata_box to initiate mrna transcription",
  "ribosome translates mrna codon sequence into polypeptide protein chain",
  "enhancer loop recruits transcription_factor and rna_polymerase to promoter",
  "crispr cas9 endonuclease cleaves target chromatin locus guided by grna",
  "histone acetylation opens euchromatin so transcription_factor binds promoter",
  "spliceosome excises intron lariat from pre_mrna before ribosome translation",
  "methylation of promoter cpg_island silences gene transcription in chromatin",
  "ribosome pairs trna anticodon with mrna codon during polypeptide elongation",
  "crispr cas9 editing introduces indel mutation at target chromatin exon",
  "enhancer and promoter coordinate rna_polymerase burst of mrna transcription",
  "histone deacetylase condenses heterochromatin blocking transcription_factor binding",
  "spliceosome joins exon junctions producing mature mrna for ribosome"
)

## Curated target vocabulary (36 domain keywords across the 4 themes)
vocab_themes <- c(
  ## Coral Reef (9)
  coral = "Coral Reef", polyp = "Coral Reef", zooxanthellae = "Coral Reef",
  reef = "Coral Reef", lagoon = "Coral Reef", parrotfish = "Coral Reef",
  clownfish = "Coral Reef", anemone = "Coral Reef", bleaching = "Coral Reef",
  ## Sonoran Desert (9)
  saguaro = "Sonoran Desert", cactus = "Sonoran Desert", succulent = "Sonoran Desert",
  arid = "Sonoran Desert", kangaroo_rat = "Sonoran Desert", burrow = "Sonoran Desert",
  nocturnal = "Sonoran Desert", canyon = "Sonoran Desert", rattlesnake = "Sonoran Desert",
  ## Boreal Taiga (9)
  boreal = "Boreal Taiga", spruce = "Boreal Taiga", taiga = "Boreal Taiga",
  permafrost = "Boreal Taiga", snowpack = "Boreal Taiga", canada_lynx = "Boreal Taiga",
  snowshoe_hare = "Boreal Taiga", caribou = "Boreal Taiga", lichen = "Boreal Taiga",
  ## Molecular Genetics (9)
  promoter = "Genetics & Chromatin", enhancer = "Genetics & Chromatin",
  rna_polymerase = "Genetics & Chromatin", transcription = "Genetics & Chromatin",
  mrna = "Genetics & Chromatin", ribosome = "Genetics & Chromatin",
  codon = "Genetics & Chromatin", chromatin = "Genetics & Chromatin",
  crispr = "Genetics & Chromatin"
)
vocab <- names(vocab_themes)
V <- length(vocab)

## Build symmetric sliding-window co-occurrence matrix C (window = +/- 4 tokens)
C_mat <- matrix(0, nrow = V, ncol = V, dimnames = list(vocab, vocab))
window_size <- 4L

for (doc in corpus) {
  tokens <- unlist(strsplit(tolower(doc), "[^a-z0-9_]+"))
  tokens <- tokens[tokens %in% vocab]
  n_tok  <- length(tokens)
  if (n_tok < 2L) next
  for (i in seq_len(n_tok)) {
    w_i <- tokens[i]
    j_idx <- setdiff(max(1L, i - window_size):min(n_tok, i + window_size), i)
    for (j in j_idx) {
      w_j <- tokens[j]
      ## Weight closer words slightly higher (1 / distance)
      C_mat[w_i, w_j] <- C_mat[w_i, w_j] + 1 / abs(i - j)
    }
  }
}

## Positive Pointwise Mutual Information: PPMI(w, c) = max(0, log2(P(w,c) / (P(w)*P(c))))
total_cooc <- sum(C_mat)
P_wc <- C_mat / total_cooc
P_w  <- rowSums(P_wc)
P_c  <- colSums(P_wc)
denom <- outer(P_w, P_c)
PMI <- log2((P_wc + 1e-12) / (denom + 1e-12))
PPMI <- pmax(PMI, 0)

## Truncated SVD of PPMI gives k-dimensional distributional embeddings:
##   E_corpus = U_k %*% diag(sqrt(d_k))
## With 4 guilds of 9 words, each guild spans ~3-4 singular vectors, so k = 16
## captures the primary and secondary semantic axes of all four guilds.
k_dim <- 16L
sv <- svd(PPMI, nu = k_dim, nv = k_dim)
E_corpus <- l2_normalize(sv$u %*% diag(sqrt(sv$d[seq_len(k_dim)])))
rownames(E_corpus) <- vocab

cat("\n=== PART B: Distributional Embeddings Learned from 48 Field Notes (PPMI-SVD) ===\n")
for (probe in c("coral", "kangaroo_rat", "permafrost", "ribosome")) {
  raw_sims <- as.numeric(E_corpus %*% E_corpus[probe, ])
  names(raw_sims) <- rownames(E_corpus)
  sims <- sort(raw_sims[names(raw_sims) != probe], decreasing = TRUE)
  cat(sprintf("  Nearest words to %-14s -> %s\n",
              paste0("'", probe, "'"),
              paste(sprintf("%s (%.2f)", names(sims)[1:3], sims[1:3]), collapse = ", ")))
}

## Project 8-D learned embeddings to 2-D with UMAP (cosine metric)
set.seed(42)
umap_B <- uwot::umap(E_corpus, n_neighbors = 7, min_dist = 0.25,
                     metric = "cosine", n_components = 2)

df_B <- data.frame(
  word   = vocab,
  theme  = factor(vocab_themes, levels = unique(vocab_themes)),
  UMAP1  = umap_B[, 1],
  UMAP2  = umap_B[, 2],
  degree = rowSums(C_mat),
  row.names = NULL
)

theme_pal <- c(
  "Coral Reef"           = "#0072B2",
  "Sonoran Desert"       = "#D55E00",
  "Boreal Taiga"         = "#009E73",
  "Genetics & Chromatin" = "#CC79A7"
)

par(mfrow = c(1, 1))
p2 <- ggplot(df_B, aes(x = UMAP1, y = UMAP2, color = theme, fill = theme)) +
  stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.80,
               linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  geom_point(aes(size = degree), alpha = 0.9) +
  ggrepel::geom_text_repel(
    aes(label = word),
    size = 3.4, fontface = "bold", box.padding = 0.35,
    point.padding = 0.25, max.overlaps = Inf, seed = 42,
    show.legend = FALSE
  ) +
  scale_color_manual(values = theme_pal, name = "Corpus Guild") +
  scale_fill_manual(values = theme_pal) +
  scale_size_continuous(range = c(2.2, 5.2), name = "Co-occurrence\nWeight") +
  labs(
    title    = "Unsupervised Word Embeddings Learned from 48 Biological Field Notes",
    subtitle = "Sliding-window co-occurrence -> Positive PMI (PPMI) -> 8-D Truncated SVD -> Cosine UMAP",
    x        = "UMAP 1 (cosine manifold)",
    y        = "UMAP 2 (cosine manifold)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )
print(p2)


## ############################################################################
##  PART C — BIOLOGICAL SEQUENCE EMBEDDINGS ("PROTEIN LANGUAGE MODEL" TOY DATA)
## ############################################################################
##  Protein Language Models (ESM-2, ProtTrans) embed amino-acid sequences into
##  fixed-length vectors where cosine similarity reflects structure & function.
##  Here we build a toy dataset of 30 peptides across 5 functional families,
##  compute a 28-D composition + physicochemical + motif embedding, and perform
##  zero-shot functional classification on 3 newly sequenced mystery peptides.

bio_sequences <- data.frame(
  id = c(
    ## 1. Transmembrane Helices (hydrophobic L/I/V/F/A/W cores)
    "TM_Rhodopsin_H1", "TM_Aquaporin_H2", "TM_KvChannel_S1",
    "TM_GPCR_TM3",     "TM_SecY_Plug",    "TM_ATPase_cSub",
    ## 2. Nuclear / DNA-Binding Basic Motifs (R/K/H rich)
    "NLS_SV40_LargeT", "NLS_Nucleoplasmin", "DBD_bZIP_GCN4",
    "DBD_Homeobox_H3", "NLS_cMyc_Core",     "DBD_Histone_H4Tail",
    ## 3. Collagen-like Extracellular Repeats (G-P-X triplet repeats)
    "COL1A1_TripleHelix", "COL2A1_Cartilage", "COL3A1_Vascular",
    "COL4A1_BasalLamina", "C1q_Stalk_Repeat", "Ficolin_CollagenStem",
    ## 4. Acidic Activation / Ca2+-Binding Loops (D/E rich)
    "TAD_VP16_Acidic",   "TAD_p53_Nterm",     "CaBP_Calmodulin_EF1",
    "CaBP_TroponinC_EF", "Chaperone_Calseq",  "HistoneChaperone_NAP",
    ## 5. Disulfide-Rich Venom / Defensin Peptides (C-rich + aromatic W/Y)
    "Toxin_AlphaCobratoxin", "Toxin_Conotoxin_Omega", "Defensin_Beta_hBD1",
    "Toxin_Charybdotoxin",   "Defensin_Plant_RsAFP",  "Toxin_FunnelWeb_HwTx"
  ),
  family = rep(c(
    "Transmembrane Helix",
    "DNA-Binding / NLS",
    "Collagen Triple-Helix",
    "Acidic / Ca2+-Binding",
    "Disulfide Venom / Defensin"
  ), each = 6),
  seq = c(
    ## Transmembrane Helices
    "WQFSMLAAYMFLLIVLGFPINFLTLYVTV",
    "GIALAFGLSIATLAQSVGHISGAHLNPAV",
    "LVVLLFFIVGVVVYFAETVPLIVIFFIVL",
    "LAAIADLLFVFTLPFWAVYAAHDWLFGFL",
    "FILLIFFVGLIIAVLIYFLGVVPIAVLLI",
    "AAIAAVGLGALGAGIGNGIVGAAIAGIGL",
    ## DNA-Binding / NLS
    "PKKKRKVPKKKRKVGGRKKRRQRRRPPQ",
    "KRPAATKKAGQAKKKKLDKRKRRKNRKAA",
    "KRARNTEAARRSRARKLQRMKQLEDKVEE",
    "RKKRCPYTKYQTLELEKEFLFNMYLTRDR",
    "PAAKRVKLDQRRRKNNHKKRKRSRKAAKK",
    "SGRGKGGKGLGKGGAKRHRKVLRDNIQGI",
    ## Collagen Triple-Helix
    "GPPGPPGPPGPAGPAGPPGPPGPQGPPGP",
    "GPPGPQGATGPLGPKGQTGEPGIAGFKGE",
    "GPPGPPGTAGFPGSPGAKGEVGPAGSPGS",
    "GFPGMKGDPGLPGQPGLPGPPGPKGDPGI",
    "GPPGEPGPRGPPGPPGNSGEKGDQGPAGP",
    "GPPGPPGEKGAPGQAGPPGPPGPAGPKGP",
    ## Acidic / Ca2+-Binding
    "APPTDVSLGDELHLDGEDVAMAHADALDD",
    "MEEPQSDPSVEPPLSQETFSDLWKLLPEN",
    "DKDGDGTITTKELGTVMRSLGQNPTEAEL",
    "DEDGSGTIDFEEFLVMMVRQMKEDAKGKD",
    "EEDDDDEDDDGDDDEEDEEVAAEEEDDDD",
    "DEDEDEDEDAEAEAEADEDEDEDEEEDDD",
    ## Disulfide Venom / Defensin
    "IRCFITPDITSKDCPNGHVCYTKTWCDAFC",
    "CKGKGAKCSRLMYDCCTGSCRSGKCWCYGC",
    "DHYNCVSSGGQCLYSACPIFTKIQGTCYRG",
    "QFTNVSCTTSKECWSVCQRLHNTSRGKCMN",
    "QKLCQRPSGTWSGVCGNNNACKNQCIRLEK",
    "ECLEIFKACNPSNDQCCKSSKLVCSRKTRW"
  ),
  stringsAsFactors = FALSE
)

## Three newly sequenced "mystery" peptides for zero-shot classification
mystery_peptides <- data.frame(
  id = c("Query1_DeepSeaVent_ORF", "Query2_SpiderGland_EST", "Query3_ChromatinIP_Hit"),
  family = c("Mystery Query", "Mystery Query", "Mystery Query"),
  seq = c(
    "VLLIIFGVAAVLGLIFVLLPIAVGLFILV",   # clearly a hydrophobic transmembrane helix
    "CRWRGSKCCGRKcrCSGKCWYCYGCKCSC",   # disulfide-rich + aromatic toxin
    "KKRRKRKPKKRAATKKRKRQRRRKKLRKA"    # basic nuclear / DNA-binding motif
  ),
  stringsAsFactors = FALSE
)
mystery_peptides$seq <- toupper(mystery_peptides$seq)

## Feature extractor: 20 AA frequencies + 8 physicochemical & motif descriptors
aa_alphabet <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
## Kyte-Doolittle hydropathy scale
kd_scale <- c(
  A =  1.8, C =  2.5, D = -3.5, E = -3.5, F =  2.8,
  G = -0.4, H = -3.2, I =  4.5, K = -3.9, L =  3.8,
  M =  1.9, N = -3.5, P = -1.6, Q = -3.5, R = -4.5,
  S = -0.8, T = -0.7, V =  4.2, W = -0.9, Y = -1.3
)

embed_sequences <- function(seqs) {
  out <- t(vapply(seqs, function(s) {
    chars <- strsplit(s, "")[[1]]
    n <- length(chars)
    ## 1. 20 standard amino acid frequencies
    aa_freq <- vapply(aa_alphabet, function(a) sum(chars == a) / n, numeric(1))
    ## 2. Physicochemical & k-mer summary descriptors (scaled to comparable magnitude)
    hydropathy <- mean(kd_scale[chars], na.rm = TRUE) / 4.5
    pos_charge <- sum(chars %in% c("K", "R", "H")) / n
    neg_charge <- sum(chars %in% c("D", "E")) / n
    net_charge <- pos_charge - neg_charge
    aromatic   <- sum(chars %in% c("F", "W", "Y")) / n
    aliphatic  <- sum(chars %in% c("I", "L", "V", "A")) / n
    ## Dipeptide motif densities: GP (collagen), KK/RR/KR/RK (basic clusters)
    bigrams    <- paste0(chars[-n], chars[-1])
    gp_density <- sum(bigrams %in% c("GP", "PG")) / (n - 1)
    rk_cluster <- sum(bigrams %in% c("KK", "RR", "KR", "RK")) / (n - 1)
    c(aa_freq,
      hydro = hydropathy, pos = pos_charge, neg = neg_charge,
      net_q = net_charge, arom = aromatic, aliph = aliphatic,
      gp_di = gp_density, rk_di = rk_cluster)
  }, numeric(28)))
  l2_normalize(out)
}

E_prot  <- embed_sequences(bio_sequences$seq)
rownames(E_prot) <- bio_sequences$id

E_query <- embed_sequences(mystery_peptides$seq)
rownames(E_query) <- mystery_peptides$id

## Zero-shot classification of the 3 mystery peptides by cosine similarity
cat("\n=== PART C: Zero-Shot Functional Classification of Mystery Peptides ===\n")
sim_matrix <- E_query %*% t(E_prot)
for (i in seq_len(nrow(mystery_peptides))) {
  q_id <- mystery_peptides$id[i]
  best_idx <- order(sim_matrix[i, ], decreasing = TRUE)[1:2]
  cat(sprintf("  %-24s -> Predicted Family: %-26s (nearest: %s, cos = %.3f)\n",
              q_id,
              bio_sequences$family[best_idx[1]],
              bio_sequences$id[best_idx[1]],
              sim_matrix[i, best_idx[1]]))
}

## Project Reference + Query sequences into a shared PCA space
pc_C <- prcomp(E_prot, center = TRUE, scale. = FALSE)
ve_C <- 100 * (pc_C$sdev^2 / sum(pc_C$sdev^2))
query_pc <- predict(pc_C, newdata = E_query)

df_C_ref <- data.frame(
  id     = bio_sequences$id,
  family = bio_sequences$family,
  PC1    = pc_C$x[, 1],
  PC2    = pc_C$x[, 2],
  is_query = FALSE
)
df_C_qry <- data.frame(
  id     = mystery_peptides$id,
  family = "Mystery Query",
  PC1    = query_pc[, 1],
  PC2    = query_pc[, 2],
  is_query = TRUE
)

prot_pal <- c(
  "Transmembrane Helix"        = "#1f3b73",
  "DNA-Binding / NLS"          = "#2a7f7f",
  "Collagen Triple-Helix"      = "#c88a2a",
  "Acidic / Ca2+-Binding"      = "#6a4c93",
  "Disulfide Venom / Defensin" = "#8c2d3a"
)

par(mfrow = c(1, 1))
p3 <- ggplot(df_C_ref, aes(x = PC1, y = PC2, color = family, fill = family)) +
  stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.82,
               linetype = "dotted", show.legend = FALSE) +
  geom_point(size = 3.0, alpha = 0.9) +
  ggrepel::geom_text_repel(
    aes(label = id),
    size = 3.0, box.padding = 0.3, point.padding = 0.2,
    max.overlaps = Inf, seed = 42, show.legend = FALSE
  ) +
  ## Overlay the 3 mystery queries as prominent black-rimmed diamonds
  geom_point(
    data = df_C_qry,
    aes(x = PC1, y = PC2),
    inherit.aes = FALSE,
    shape = 23, size = 4.8, fill = "#FFD700", color = "black", stroke = 1.1
  ) +
  ggrepel::geom_label_repel(
    data = df_C_qry,
    aes(x = PC1, y = PC2, label = id),
    inherit.aes = FALSE,
    size = 3.3, fontface = "bold", fill = "#FFF9DB", color = "black",
    box.padding = 0.5, point.padding = 0.35, max.overlaps = Inf, seed = 42
  ) +
  scale_color_manual(values = prot_pal, name = "Protein Family") +
  scale_fill_manual(values = prot_pal) +
  labs(
    title    = "Biological Sequence Embeddings: 30 Functional Peptides + 3 Zero-Shot Queries",
    subtitle = "28-D AA composition + hydropathy + charge + motif embedding (gold diamonds = unclassified mystery ORFs)",
    x        = sprintf("PC1 (%.1f%% variance)", ve_C[1]),
    y        = sprintf("PC2 (%.1f%% variance)", ve_C[2])
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )
print(p3)


## ############################################################################
##  PART D — INSIDE THE TRANSFORMER: SELF-ATTENTION & POLYSEMY DISAMBIGUATION
## ############################################################################
##  In static embeddings (Word2Vec / GloVe), a polysemous word like "expression",
##  "vector", or "culture" has ONE fixed vector — an awkward average of all its
##  meanings. Vaswani et al. (2017) replaced static lookup tables with
##  Scaled Dot-Product Self-Attention:
##
##      Attention(Q, K, V) = softmax( (Q %*% t(K)) / sqrt(d_k) ) %*% V
##
##  Each token attends to surrounding tokens in the sentence and updates its
##  embedding into a *contextualized* vector.

## Build a 16-D basis with two distinct domains:
##   Domain 1: Molecular Genetics (promoter, viral, gene, transcription, cell)
##   Domain 2: Primate Behavioral Ecology (wind, vocal, facial, troop, savanna)
set.seed(101)
d_attn <- 16
Q_attn <- qr.Q(qr(matrix(rnorm(d_attn * 10), nrow = d_attn, ncol = 10)))
dom_mol   <- Q_attn[, 1]
dom_behav <- Q_attn[, 2]

## Context-specific anchor words + 3 polysemous words ("vector", "expression", "culture")
## Notice that each polysemous static embedding sits halfway between Molecular & Behavioral!
E_static <- l2_normalize(rbind(
  ## Pure Molecular anchors
  viral         = 1.4 * dom_mol + 0.25 * Q_attn[, 3],
  promoter      = 1.4 * dom_mol + 0.25 * Q_attn[, 4],
  gene          = 1.4 * dom_mol + 0.25 * Q_attn[, 5],
  cell          = 1.4 * dom_mol + 0.25 * Q_attn[, 6],
  ## Pure Behavioral / Ecological anchors
  wind          = 1.4 * dom_behav + 0.25 * Q_attn[, 7],
  facial        = 1.4 * dom_behav + 0.25 * Q_attn[, 8],
  vocal         = 1.4 * dom_behav + 0.25 * Q_attn[, 9],
  primate       = 1.4 * dom_behav + 0.25 * Q_attn[, 10],
  ## Polysemous words: equal blend of Molecular and Behavioral senses + distinct lexical id
  vector        = 0.65 * dom_mol + 0.65 * dom_behav + 0.85 * Q_attn[, 3],
  expression    = 0.65 * dom_mol + 0.65 * dom_behav + 0.85 * Q_attn[, 5],
  culture       = 0.65 * dom_mol + 0.65 * dom_behav + 0.85 * Q_attn[, 6]
))

row_softmax <- function(M) {
  M_shift <- M - apply(M, 1, max)
  E_exp   <- exp(M_shift)
  E_exp / rowSums(E_exp)
}

## One layer of Scaled Dot-Product Self-Attention (Vaswani et al. 2017)
self_attention <- function(tokens, emb = E_static, scale_temp = 0.25) {
  X <- emb[tokens, , drop = FALSE]
  ## Queries and Keys project onto the semantic domain subspace so tokens
  ## attend strongly to domain-informative context words in the sentence
  P_domain <- outer(dom_mol, dom_mol) + outer(dom_behav, dom_behav)
  Q <- X %*% (0.15 * diag(ncol(X)) + 2.5 * P_domain)
  K <- Q
  V <- X
  scores  <- (Q %*% t(K)) / (sqrt(ncol(X)) * scale_temp)
  weights <- row_softmax(scores)
  context_emb <- l2_normalize(0.25 * X + 0.75 * (weights %*% V))
  rownames(context_emb) <- tokens
  list(embeddings = context_emb, weights = weights)
}

sent_mol   <- c("viral", "vector", "promoter", "gene", "expression", "cell", "culture")
sent_behav <- c("wind",  "vector", "facial",   "vocal", "expression", "primate", "culture")

attn_mol   <- self_attention(sent_mol)
attn_behav <- self_attention(sent_behav)

cat("\n=== PART D: Contextual Disambiguation via Self-Attention ===\n")
for (poly in c("vector", "expression", "culture")) {
  cos_before <- cosine_sim(E_static[poly, ], E_static["gene", ]) -
                cosine_sim(E_static[poly, ], E_static["primate", ])
  cos_mol    <- cosine_sim(attn_mol$embeddings[poly, ], E_static["gene", ]) -
                cosine_sim(attn_mol$embeddings[poly, ], E_static["primate", ])
  cos_behav  <- cosine_sim(attn_behav$embeddings[poly, ], E_static["gene", ]) -
                cosine_sim(attn_behav$embeddings[poly, ], E_static["primate", ])
  cat(sprintf("  Token %-12s | Mol-vs-Behav affinity: static = %+5.2f | in Mol sentence = %+5.2f | in Behav sentence = %+5.2f\n",
              paste0("'", poly, "'"), cos_before, cos_mol, cos_behav))
}

## ---- Figure 4a: Self-Attention Weight Heatmap ------------------------------
W_mol <- attn_mol$weights
df_heat <- expand.grid(
  Query = factor(sent_mol, levels = rev(sent_mol)),
  Key   = factor(sent_mol, levels = sent_mol)
)
df_heat$Weight <- W_mol[cbind(as.character(df_heat$Query), as.character(df_heat$Key))]

par(mfrow = c(1, 1))
p4a <- ggplot(df_heat, aes(x = Key, y = Query, fill = Weight)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", Weight)),
            color = ifelse(df_heat$Weight > 0.22, "white", "grey15"),
            size = 3.6, fontface = "bold") +
  scale_fill_gradient(low = "#f4f6f9", high = "#1f3b73",
                      name = "Attention\nWeight") +
  labs(
    title    = "Scaled Dot-Product Self-Attention Weights: softmax(Q K^T / sqrt(d_k))",
    subtitle = "Sentence: 'viral vector promoter gene expression cell culture'",
    x        = "Key token attended to (K)",
    y        = "Query token being updated (Q)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )
print(p4a)

## ---- Figure 4b: Polysemous Words Splitting in PCA Space After Attention ----
poly_words <- c("vector", "expression", "culture")
E_combined <- rbind(
  E_static,
  attn_mol$embeddings[poly_words, ],
  attn_behav$embeddings[poly_words, ]
)
rownames(E_combined) <- c(
  rownames(E_static),
  paste0(poly_words, " (in Mol context)"),
  paste0(poly_words, " (in Behav context)")
)

pc_D <- prcomp(E_combined, center = TRUE)
ve_D <- 100 * (pc_D$sdev^2 / sum(pc_D$sdev^2))

role_D <- c(
  rep("Molecular Anchor", 4),
  rep("Behavioral Anchor", 4),
  rep("Polysemous (Static Before Attention)", 3),
  rep("Contextualized (Molecular Sentence)", 3),
  rep("Contextualized (Behavioral Sentence)", 3)
)

df_D <- data.frame(
  label = rownames(E_combined),
  role  = factor(role_D, levels = unique(role_D)),
  PC1   = pc_D$x[, 1],
  PC2   = pc_D$x[, 2],
  row.names = NULL
)

## Build trajectory table from static polysemous word -> its two contextual senses
traj_df <- rbind(
  data.frame(
    word    = poly_words,
    context = "Molecular Sentence",
    x_start = pc_D$x[poly_words, 1],
    y_start = pc_D$x[poly_words, 2],
    x_end   = pc_D$x[paste0(poly_words, " (in Mol context)"), 1],
    y_end   = pc_D$x[paste0(poly_words, " (in Mol context)"), 2]
  ),
  data.frame(
    word    = poly_words,
    context = "Behavioral Sentence",
    x_start = pc_D$x[poly_words, 1],
    y_start = pc_D$x[poly_words, 2],
    x_end   = pc_D$x[paste0(poly_words, " (in Behav context)"), 1],
    y_end   = pc_D$x[paste0(poly_words, " (in Behav context)"), 2]
  )
)

attn_pal <- c(
  "Molecular Anchor"                     = "#1f3b73",
  "Behavioral Anchor"                    = "#c88a2a",
  "Polysemous (Static Before Attention)" = "#7f7f7f",
  "Contextualized (Molecular Sentence)"  = "#2b8cbe",
  "Contextualized (Behavioral Sentence)" = "#e34a33"
)

par(mfrow = c(1, 1))
p4b <- ggplot(df_D, aes(x = PC1, y = PC2, color = role)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey85") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey85") +
  geom_curve(
    data = traj_df,
    aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
    inherit.aes = FALSE,
    curvature = 0.18, color = "grey30", linewidth = 0.75,
    arrow = grid::arrow(length = grid::unit(0.20, "cm"), type = "closed")
  ) +
  geom_point(aes(shape = role), size = 3.6, alpha = 0.95) +
  ggrepel::geom_text_repel(
    aes(label = label),
    size = 3.4, fontface = "bold", box.padding = 0.4,
    point.padding = 0.25, max.overlaps = Inf, seed = 42,
    show.legend = FALSE
  ) +
  scale_color_manual(values = attn_pal, name = "Token State") +
  scale_shape_manual(values = c(16, 16, 15, 17, 17), name = "Token State") +
  labs(
    title    = "How Self-Attention Resolves Polysemy: Static vs. Contextual Embeddings",
    subtitle = "Ambiguous tokens ('vector', 'expression', 'culture') split from a static midpoint toward their sentence context",
    x        = sprintf("PC1 (%.1f%% variance -- Molecular vs. Behavioral axis)", ve_D[1]),
    y        = sprintf("PC2 (%.1f%% variance)", ve_D[2])
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )
print(p4b)

cat("\n[18_LLM_Embeddings] done. Embeddings = multivariate vectors; reuse the",
    "whole toolbox.\n")
