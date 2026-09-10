make_all_coexpression_heatmaps <- function(THYGNS) {
       
    thyset.RFK <- DG[ DG[[THYGNS]] == 1, c("geneid","symbol","name","set") ];
    S <- length(thyset.RFK$geneid);
    SETS <- unique(thyset.RFK$set);

    cat("###\n### WORKING on ",THYGNS,"\n");
         
    # Heatmaps by REGION
    DO <- DT.Red[[THYGNS]] %>%
          select(c(thyset.RFK$geneid, "region_label"));
    TC <- nrow(DO);
                      
    CO.Red <- list();
    CN.Red <- list();
    for (i in HTMP.order) {
      # cat("# ", i, "\n");
      CO.Red[[i]] <- matrix(0L, nrow=S, ncol=S, dimnames=list(thyset.RFK$geneid, thyset.RFK$geneid));
      CN.Red[[i]] <- 0;
    };
    #
    for (d in 1:nrow(DO)) {
      V <- DO[d, ];
      RG <- V$region_label;
      CN.Red[[RG]] <- CN.Red[[RG]] + 1;
      for (i in 1:(S - 1)) {
        G1 <- thyset.RFK$geneid[[i]];
        V1 <- V[[G1]];
        if (G1 > 0) {
          CO.Red[[RG]][G1,G1] <- CO.Red[[RG]][G1,G1] + 1;
        };
        for (j in (i+1):S) {
          G2 <- thyset.RFK$geneid[[j]];
          V2 <- V[[G2]];
          if (V1 > 0 & V2 > 0) {
            CO.Red[[RG]][G1,G2] <- CO.Red[[RG]][G1,G2] + 1;
            CO.Red[[RG]][G2,G1] <- CO.Red[[RG]][G2,G1] + V1 + V2;
          };
        }; # j
      }; # i
      if (G2 > 0) {
        CO.Red[[RG]][G2,G2] <- CO.Red[[RG]][G2,G2] + 1;
      };
    }; # d
    CA.Red <- list();
    for (d in HTMP.order) {
      CA.Red[[d]] <- matrix(0L, nrow=S, ncol=S, dimnames=list(thyset.RFK$geneid, thyset.RFK$geneid));
      for (i in 1:S) {
        G1 <- thyset.RFK$geneid[[i]];
        for (j in 1:S) {
          G2 <- thyset.RFK$geneid[[j]];
          if (i <= j ) {
            CA.Red[[d]][G1,G2] <- CO.Red[[d]][G1,G2];
          } else { # i > j
            CA.Red[[d]][G1,G2] <- CO.Red[[d]][G1,G2]/(CO.Red[[d]][G2,G1]*2);
          };
        };
      };
    };

    
    for (d in HTMP.order) {
        cat("### Reddien REGION co-expression heatmaps: ",d,"\n");
  
        m <- CA.Red[[d]];
        m[ is.nan(m) ] <- NA;
        m[ is.infinite(m) ] <- NA;
        mu.min <- min(m[upper.tri(m)],na.rm=T);
        mu.max <- max(m[upper.tri(m)],na.rm=T);
        ml.min <- min(m[lower.tri(m)],na.rm=T);
        ml.max <- max(m[lower.tri(m)],na.rm=T);
        if ((is.infinite(mu.max) | is.infinite(mu.min)) | mu.max - mu.min == 0) {
            colrange1 <- as.numeric(map2(0, 1, seq, length.out=9)[[1]]);
        } else {
            colrange1 <- as.numeric(map2(mu.min, mu.max, seq, length.out=9)[[1]]);
        };
        if ((is.infinite(ml.max) | is.infinite(ml.min)) | ml.max - ml.min == 0) {
            colrange2 <- as.numeric(map2(0, 1, seq, length.out=256)[[1]]);
        } else {
            colrange2 <- as.numeric(map2(ml.min, ml.max, seq, length.out=256)[[1]]);
        };
        col1 <- colorRamp2(colrange1, brewer.pal(9, "Blues"));
        col2 <- colorRamp2(colrange2, paletteer::paletteer_c("viridis::inferno",256));
                                      # grDevices::Purple-Yellow # brewer.pal(9, "PuRd"));
        mmx <- ifelse(is.infinite(mu.max) | mu.max == 0, 1, mu.max);
        m.rows <- nrow(m);
        m.cols <- ncol(m);
        if (sum(m, na.rm=T) == 0) {
            cat("#---> Skipping empty matrix\n");
            next;
        };
        ht <- Heatmap(m, name = "coexpression",
                      rect_gp = gpar(type = "none"),
                      show_heatmap_legend = FALSE,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      row_labels = sprintf("%-25s",   thyset.RFK[ rownames(m), "symbol" ]),
                      column_labels = sprintf("%20s", thyset.RFK[ colnames(m), "symbol" ]),
                      cell_fun = function(j, i, x, y, w, h, fill) {
                          if (i < j) { # l = i < j
                              grid.rect(x, y, w, h, 
                                        gp = gpar(fill = ifelse(is.na(m[j,i]) | is.na(m[i,j]),
                                                                "#CCCCCC", col1(m[i,j])),
                                                  col = NA
                                                  ))
                              grid.text(sprintf("%s", ifelse(is.na(m[i,j]), 0, m[i,j])),
                                        x, y,
                                        gp = gpar(fontsize = 10,
                                                  col = ifelse(i==j | (!is.na(m[i,j]) & m[i,j] < 0.65*mmx),
                                                               "black", "white")
                                                  ))
                          } else if (i > j) {
                              grid.rect(x, y, w, h, 
                                        gp = gpar(fill = ifelse(is.na(m[i,j]),
                                                                "#CCCCCC", col2(m[i,j])),
                                                  col = NA))
                          }
                      },
                      na_col = "#CCCCCC"
                      );
        pdf(paste0(WD,"/coexpression_heatmaps/Reddien_coexpression_heatmap_",THYGNS,"_by_region_",d,".pdf"),
            paper = "special",
            width  = unit((0.55 * m.cols + 14) / 2.54, "cm"), # 1in~2.54cm
            height = unit((0.55 * m.rows +  8) / 2.54, "cm"),
            bg = "white");
        draw(ht,
             heatmap_legend_list = list(
                 Legend(title = "#Cells\nCo-expressed", col_fun = col1),
                 Legend(title = "Average\nExpression", col_fun = col2)
             ));
        decorate_heatmap_body("coexpression", {
                grid.text(d,
                          unit(1.015, "npc"), unit(-0.1, "npc"), hjust = 0, vjust = -0.25,
                          gp = gpar(fontsize = ifelse(S < 20, 12, 20), col = "black", fontface = "bold")
                          )
                grid.text(paste0(sprintf(  "Total #cells %6d", CN.Red[[d]]),
                                 sprintf("\n          of %6d", TC)),
                          unit(1.015, "npc"), unit(-0.1, "npc"), hjust = 0, vjust = 1.25,
                          gp = gpar(fontsize = ifelse(S < 20,  9, 16), fontfamily = "mono", col = "black")
                          )
            });
        dev.off();

    }; # for d (REGION)

         
    # Heatmaps by Cell type (Reddien) 
    DO <- DT.Red[[THYGNS]] %>%
          select(c(thyset.RFK$geneid, "label_short"));
    TC <- nrow(DO);
                  
    CO.Red <- list();
    CN.Red <- list();
    for (i in Red.CELL.ord) {
      if (i == "Other") next;
      #cat("# ", i, "\n");
      CO.Red[[i]] <- matrix(0L, nrow=S, ncol=S, dimnames=list(thyset.RFK$geneid, thyset.RFK$geneid));
      CN.Red[[i]] <- 0;
    };
    #
    for (d in 1:nrow(DO)) {
      V <- DO[d, ];
      RG <- V$label_short;
      CN.Red[[RG]] <- CN.Red[[RG]] + 1;
      for (i in 1:(S - 1)) {
        G1 <- thyset.RFK$geneid[[i]];
        V1 <- V[[G1]];
        if (G1 > 0) {
          CO.Red[[RG]][G1,G1] <- CO.Red[[RG]][G1,G1] + 1;
        };
        for (j in (i+1):S) {
          G2 <- thyset.RFK$geneid[[j]];
          V2 <- V[[G2]];
          if (V1 > 0 & V2 > 0) {
            CO.Red[[RG]][G1,G2] <- CO.Red[[RG]][G1,G2] + 1;
            CO.Red[[RG]][G2,G1] <- CO.Red[[RG]][G2,G1] + V1 + V2;
          };
        }; # j
      }; # i
      if (G2 > 0) {
        CO.Red[[RG]][G2,G2] <- CO.Red[[RG]][G2,G2] + 1;
      };
    }; # d
    CA.Red <- list();
    for (d in Red.CELL.ord) {
      if (d == "Other") next;
      CA.Red[[d]] <- matrix(0L, nrow=S, ncol=S, dimnames=list(thyset.RFK$geneid, thyset.RFK$geneid));
      for (i in 1:S) {
        G1 <- thyset.RFK$geneid[[i]];
        for (j in 1:S) {
          G2 <- thyset.RFK$geneid[[j]];
          if (i <= j ) {
            CA.Red[[d]][G1,G2] <- CO.Red[[d]][G1,G2];
          } else { # i > j
            CA.Red[[d]][G1,G2] <- CO.Red[[d]][G1,G2]/(CO.Red[[d]][G2,G1]*2);
          };
        };
      };
    };

    
    for (d in Red.CELL.ord) {
        if (d == "Other") next;
        cat("### Reddien CELL TYPE co-expression heatmaps: ",d,"\n");
  
        m <- CA.Red[[d]];
        m[ is.nan(m) ] <- NA;
        m[ is.infinite(m) ] <- NA;
        mu.min <- min(m[upper.tri(m)],na.rm=T);
        mu.max <- max(m[upper.tri(m)],na.rm=T);
        ml.min <- min(m[lower.tri(m)],na.rm=T);
        ml.max <- max(m[lower.tri(m)],na.rm=T);
        if ((is.infinite(mu.max) | is.infinite(mu.min)) | mu.max - mu.min == 0) {
            colrange1 <- as.numeric(map2(0, 1, seq, length.out=9)[[1]]);
        } else {
            colrange1 <- as.numeric(map2(mu.min, mu.max, seq, length.out=9)[[1]]);
        };
        if ((is.infinite(ml.max) | is.infinite(ml.min)) | ml.max - ml.min == 0) {
            colrange2 <- as.numeric(map2(0, 1, seq, length.out=256)[[1]]);
        } else {
            colrange2 <- as.numeric(map2(ml.min, ml.max, seq, length.out=256)[[1]]);
        };
        col1 <- colorRamp2(colrange1, brewer.pal(9, "Blues"));
        col2 <- colorRamp2(colrange2, paletteer::paletteer_c("viridis::inferno",256));
                                      # grDevices::Purple-Yellow # brewer.pal(9, "PuRd"));
        mmx <- ifelse(is.infinite(mu.max) | mu.max == 0, 1, mu.max);
        m.rows <- nrow(m);
        m.cols <- ncol(m);
        if (sum(m, na.rm=T) == 0) {
            cat("#---> Skipping empty matrix\n");
            next;
        };
        ht <- Heatmap(m, name = "coexpression",
                      rect_gp = gpar(type = "none"),
                      show_heatmap_legend = FALSE,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      row_labels = sprintf("%-25s",   thyset.RFK[ rownames(m), "symbol" ]),
                      column_labels = sprintf("%20s", thyset.RFK[ colnames(m), "symbol" ]),
                      cell_fun = function(j, i, x, y, w, h, fill) {
                          if (i < j) { # l = i < j
                              grid.rect(x, y, w, h, 
                                        gp = gpar(fill = ifelse(is.na(m[j,i]) | is.na(m[i,j]),
                                                                "#CCCCCC", col1(m[i,j])),
                                                  col = NA
                                                  ))
                              grid.text(sprintf("%s", ifelse(is.na(m[i,j]), 0, m[i,j])),
                                        x, y,
                                        gp = gpar(fontsize = 10,
                                                  col = ifelse(i==j | (!is.na(m[i,j]) & m[i,j] < 0.65*mmx),
                                                               "black", "white")
                                                  ))
                          } else if (i > j) {
                              grid.rect(x, y, w, h, 
                                        gp = gpar(fill = ifelse(is.na(m[i,j]),
                                                                "#CCCCCC", col2(m[i,j])),
                                                  col = NA))
                          }
                      },
                      na_col = "#CCCCCC"
                      );
        pdf(paste0(WD,"/coexpression_heatmaps/Reddien_coexpression_heatmap_",THYGNS,"_by_celltype_",d,".pdf"),
            paper = "special",
            width  = unit((0.55 * m.cols + 14) / 2.54, "cm"), # 1in~2.54cm
            height = unit((0.55 * m.rows +  8) / 2.54, "cm"),
            bg = "white");
        draw(ht,
             heatmap_legend_list = list(
                 Legend(title = "#Cells\nCo-expressed", col_fun = col1),
                 Legend(title = "Average\nExpression", col_fun = col2)
             ));
        decorate_heatmap_body("coexpression", {
                grid.text(d,
                          unit(1.015, "npc"), unit(-0.1, "npc"), hjust = 0, vjust = -0.25,
                          gp = gpar(fontsize = ifelse(S < 20, 12, 20), col = "black", fontface = "bold")
                          )
                grid.text(paste0(sprintf(  "Total #cells %6d", CN.Red[[d]]),
                                 sprintf("\n          of %6d", TC)),
                          unit(1.015, "npc"), unit(-0.1, "npc"), hjust = 0, vjust = 1.25,
                          gp = gpar(fontsize = ifelse(S < 20,  9, 16), fontfamily = "mono", col = "black")
                          )
            });
        dev.off();

    }; # for d (cell type)

         
    # Heatmaps by Cell type (Rajewsky) 
    DO <- DT.Raj[[THYGNS]] %>%
          select(c(thyset.RFK$geneid, "label_short"));
    TC <- nrow(DO);
                  
    CO.Raj <- list();
    CN.Raj <- list();
    for (i in Raj.CELL.ord2) {
      if (i == "Other") next;
      #cat("#a ", i, "\n");
      CO.Raj[[i]] <- matrix(0L, nrow=S, ncol=S, dimnames=list(thyset.RFK$geneid, thyset.RFK$geneid));
      CN.Raj[[i]] <- 0;
    };
    #
    for (d in 1:nrow(DO)) {
      #cat("#b ", d, "\n");
      V <- DO[d, ];
      RG <- V$label_short;
      CN.Raj[[RG]] <- CN.Raj[[RG]] + 1;
      for (i in 1:(S - 1)) {
        G1 <- thyset.RFK$geneid[[i]];
        V1 <- V[[G1]];
        if (G1 > 0) {
          CO.Raj[[RG]][G1,G1] <- CO.Raj[[RG]][G1,G1] + 1;
        };
        for (j in (i+1):S) {
          G2 <- thyset.RFK$geneid[[j]];
          V2 <- V[[G2]];
          if (V1 > 0 & V2 > 0) {
            CO.Raj[[RG]][G1,G2] <- CO.Raj[[RG]][G1,G2] + 1;
            CO.Raj[[RG]][G2,G1] <- CO.Raj[[RG]][G2,G1] + V1 + V2;
          };
        }; # j
      }; # i
      if (G2 > 0) {
        CO.Raj[[RG]][G2,G2] <- CO.Raj[[RG]][G2,G2] + 1;
      };
    }; # d
    CA.Raj <- list();
    for (d in Raj.CELL.ord2) {
      if (d == "Other") next;
      #cat("#c ", d, "\n");
      CA.Raj[[d]] <- matrix(0L, nrow=S, ncol=S, dimnames=list(thyset.RFK$geneid, thyset.RFK$geneid));
      for (i in 1:S) {
        G1 <- thyset.RFK$geneid[[i]];
        for (j in 1:S) {
          G2 <- thyset.RFK$geneid[[j]];
          if (i <= j ) {
            CA.Raj[[d]][G1,G2] <- CO.Raj[[d]][G1,G2];
          } else { # i > j
            CA.Raj[[d]][G1,G2] <- CO.Raj[[d]][G1,G2]/(CO.Raj[[d]][G2,G1]*2);
          };
        };
      };
    };

    
    for (d in Raj.CELL.ord2) {
        if (d == "Other") next;
        cat("### Rajewsky CELL TYPE co-expression heatmaps: ",d,"\n");
  
        m <- CA.Raj[[d]];
        m[ is.nan(m) ] <- NA;
        m[ is.infinite(m) ] <- NA;
        mu.min <- min(m[upper.tri(m)],na.rm=T);
        mu.max <- max(m[upper.tri(m)],na.rm=T);
        ml.min <- min(m[lower.tri(m)],na.rm=T);
        ml.max <- max(m[lower.tri(m)],na.rm=T);
        if ((is.infinite(mu.max) | is.infinite(mu.min)) | mu.max - mu.min == 0) {
            colrange1 <- as.numeric(map2(0, 1, seq, length.out=9)[[1]]);
        } else {
            colrange1 <- as.numeric(map2(mu.min, mu.max, seq, length.out=9)[[1]]);
        };
        if ((is.infinite(ml.max) | is.infinite(ml.min)) | ml.max - ml.min == 0) {
            colrange2 <- as.numeric(map2(0, 1, seq, length.out=256)[[1]]);
        } else {
            colrange2 <- as.numeric(map2(ml.min, ml.max, seq, length.out=256)[[1]]);
        };
        col1 <- colorRamp2(colrange1, brewer.pal(9, "Blues"));
        col2 <- colorRamp2(colrange2, paletteer::paletteer_c("viridis::inferno",256));
                                      # grDevices::Purple-Yellow # brewer.pal(9, "PuRd"));
        mmx <- ifelse(is.infinite(mu.max) | mu.max == 0, 1, mu.max);
        m.rows <- nrow(m);
        m.cols <- ncol(m);
        if (sum(m, na.rm=T) == 0) {
            cat("#---> Skipping empty matrix\n");
            next;
        };
        ht <- Heatmap(m, name = "coexpression",
                      rect_gp = gpar(type = "none"),
                      show_heatmap_legend = FALSE,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      row_labels = sprintf("%-25s",   thyset.RFK[ rownames(m), "symbol" ]),
                      column_labels = sprintf("%20s", thyset.RFK[ colnames(m), "symbol" ]),
                      cell_fun = function(j, i, x, y, w, h, fill) {
                          if (i < j) { # l = i < j
                              grid.rect(x, y, w, h, 
                                        gp = gpar(fill = ifelse(is.na(m[j,i]) | is.na(m[i,j]),
                                                                "#CCCCCC", col1(m[i,j])),
                                                  col = NA
                                                  ))
                              grid.text(sprintf("%s", ifelse(is.na(m[i,j]), 0, m[i,j])),
                                        x, y,
                                        gp = gpar(fontsize = 10,
                                                  col = ifelse(i==j | (!is.na(m[i,j]) & m[i,j] < 0.65*mmx),
                                                               "black", "white")
                                                  ))
                          } else if (i > j) {
                              grid.rect(x, y, w, h, 
                                        gp = gpar(fill = ifelse(is.na(m[i,j]),
                                                                "#CCCCCC", col2(m[i,j])),
                                                  col = NA))
                          }
                      },
                      na_col = "#CCCCCC"
                      );
        pdf(paste0(WD,"/coexpression_heatmaps/Rajewsky_coexpression_heatmap_",THYGNS,"_by_celltype_",d,".pdf"),
            paper = "special",
            width  = unit((0.55 * m.cols + 14) / 2.54, "cm"), # 1in~2.54cm
            height = unit((0.55 * m.rows +  8) / 2.54, "cm"),
            bg = "white");
        draw(ht,
             heatmap_legend_list = list(
                 Legend(title = "#Cells\nCo-expressed", col_fun = col1),
                 Legend(title = "Average\nExpression", col_fun = col2)
             ));
        decorate_heatmap_body("coexpression", {
                grid.text(d,
                          unit(1.015, "npc"), unit(-0.1, "npc"), hjust = 0, vjust = -0.25,
                          gp = gpar(fontsize = ifelse(S < 20, 12, 20), col = "black", fontface = "bold")
                          )
                grid.text(paste0(sprintf(  "Total #cells %6d", CN.Raj[[d]]),
                                 sprintf("\n          of %6d", TC)),
                          unit(1.015, "npc"), unit(-0.1, "npc"), hjust = 0, vjust = 1.25,
                          gp = gpar(fontsize = ifelse(S < 20,  9, 16), fontfamily = "mono", col = "black")
                          )
            });
        dev.off();

    }; # for d (cell type)

} # make_all_coexpression_heatmaps
