roundUp <- function(x) as.integer(ceiling(x/100)*100);

draw_key_line <- function(data, params, size) {
  grid::segmentsGrob(0.1, 0.5, 0.9, 0.5, 
                     gp = grid::gpar(
                              col = alpha(data$edge_colour, data$edge_alpha), 
                              lwd = data$edge_width * .pt, 
                              lty = data$edge_linetype, 
                              lineend = "butt"),
                     arrow = NULL);
} # draw_key_line
    
GeomEdgePath$draw_key <- draw_key_line;

thypg <- function(GR,THYTITLE,THYSET.RFK,NCELLS.max,AVGEXP.max) {
   thyedgecols <- unique(
       round(
           quantile(get.data.frame(GR, "E")$AVGEXP,
                    probs = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 1)),
           2));
   ggraph(GR, layout="fr") +
       geom_edge_fan(
           aes(width = NCELLS, color = AVGEXP,
               start_cap = circle(1, unit = "cm"),
               end_cap   = circle(1, unit = "cm")
               ),
           arrow = arrow(
               angle = 10,
               length = unit(0.15, "inches"),
               ends = "both",
               type = "closed"
               ),
           strength = 0.75
           ) +
       geom_node_point(
           aes(fill = GRP),
           shape = 21, color = "black", size = 20
           ) +
       geom_node_label(
           aes(label=THYSET.RFK[ name, "symbol" ]),
           size=10, show.legend = FALSE, label.padding = unit(0.1, "lines"),
           color = "red", fill = "#FFFFFFBB", label.size = NA,
           family = "Arial", fontface = "bold", repel = TRUE
           ) +
       theme_void() +
       scale_fill_manual(
           # labels = names(NODE.colors),
           breaks = names(NODE.colors),
           values = NODE.colors,
           name = "\nCell Type"
           ) +
       scale_edge_width_continuous(
           breaks = seq(0, NCELLS.max, 2.5*(10^floor(log10(NCELLS.max)-1))),
           limits = c(0, NCELLS.max),
           name = "\n#Cells Co-expressed"
           ) +
       scale_edge_color_viridis(
           option = "mako", direction = -1,
           guide = "legend",
           labels = sprintf("%04.2f", thyedgecols),
                  # seq(0, AVGEXP.max, 2*(10^floor(log10(AVGEXP.max)-1))),
           breaks = thyedgecols,
                  # seq(0, AVGEXP.max, 2*(10^floor(log10(AVGEXP.max)-1))),
           # limits = asinh(c(0, AVGEXP.max)),
           name = "\nAverage Expression"
           ) +
       coord_cartesian(clip="off") +
       ggtitle(THYTITLE) +
       guides(
           fill = guide_legend(override.aes = list(size = 15)),
           edge_color = guide_legend(override.aes = list(edge_width = 20))
           ) +
       theme(strip.text.x = element_text(size = 40, colour = "blue", 
                                         family = "Verdana", face = "bold"),
             legend.title = element_text(size = 32, face = "bold"), 
             legend.text  = element_text(size = 30),
             plot.title   = element_text(hjust = 1, size = 30),
             panel.spacing = unit(4, "lines"))
} # thypg


make_all_coexpression_graphs <- function(THYGNS) {

    thyset.RFK <- DG[ DG[[THYGNS]] == 1, c("geneid","symbol","name","set") ];
    S <- length(thyset.RFK$geneid);
    S.scale <- ifelse(S>30, 1, ifelse(S>20, 0.8, ifelse(S>10, 0.6, 0.4)));
    SETS <- unique(thyset.RFK$set);

    titletheme <- theme(plot.title = element_text(hjust = 0.5, 
                                                  size = 40 * S.scale, colour = "blue", 
                                                  family = "Verdana", face = "bold"));

    gcolnames <- c("GIDA","GIDB","NCELLS","AVGEXP","SET");

    cat("###\n### WORKING on ",THYGNS,"\n");
     
    # Graphs by REGION (Reddien) 
    DO <- DT.Red[[THYGNS]] %>%
          select(c(thyset.RFK$geneid, "region_label"));
    TC <- nrow(DO);

    CO.Red <- list();
    CN.Red <- list();
    for (i in HTMP.order) {
      # cat("# ", i, "\n");
      CO.Red[[i]] <- matrix(0L, nrow=S, ncol=S,
                            dimnames=list(thyset.RFK$geneid, thyset.RFK$geneid));
      CN.Red[[i]] <- 0;
    };
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
          if (j == S & G2 > 0) {
            CO.Red[[RG]][G2,G2] <- CO.Red[[RG]][G2,G2] + 1;
          };
        }; # j
      }; # i
    }; # d
    CA.Red <- list();
    for (d in HTMP.order) {
      CA.Red[[d]] <- matrix(0L, nrow=S, ncol=S, 
                            dimnames=list(thyset.RFK$geneid, thyset.RFK$geneid));
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

    write.table(CA.Red,
                file=paste0(WD,"/coexpression_graphs/Reddien_coexpression_matrix_",
                            THYGNS,"_byregion.tbl"),
                quote=TRUE, na="NA", col.names=TRUE, row.names=TRUE, sep="\t");

    G <- data.frame();
    for (d in HTMP.order) {
      for (i in 1:S) {
        G1 <- thyset.RFK$geneid[[i]];
        for (j in 1:S) {
          G2 <- thyset.RFK$geneid[[j]];
          if (i<j) {
            Ncell <- CA.Red[[d]][G1,G2];
            Nexpr <- CA.Red[[d]][G2,G1];
            # cat("Ncell: ",Ncell,"  Nexpr: ",Nexpr,"\n");
            if (!is.null(Ncell) & Ncell > 0) {
              if (nrow(G) > 0) {
                G <- rbind(G, c(G1, G2, Ncell, Nexpr, d));
              } else {
                G <- rbind(G, c(G1, G2, Ncell, Nexpr, d));
                colnames(G) <- gcolnames;
              };
            };
          };
        };
      };
    };

    G$NCELLS <- as.numeric(G$NCELLS);
    G$AVGEXP <- as.numeric(G$AVGEXP);

    write.table(G,
                file=paste0(WD,"/coexpression_graphs/Reddien_coexpression_matrix_",
                            THYGNS,"_byregion_alledges.tbl"),
                quote=TRUE, na="NA", col.names=TRUE, row.names=FALSE, sep="\t");

    cat("#--> REGION\n");

    ncells.avg <- mean(G$NCELLS, na.rm=TRUE);
    avgexp.avg <- mean(G$AVGEXP, na.rm=TRUE);
    ncells.max <- roundUp(max(G$NCELLS, na.rm=TRUE));
    avgexp.max <- roundUp(max(G$AVGEXP, na.rm=TRUE));

    gr.all_pre <- graph_from_data_frame(G, directed=TRUE);
    df_edges   <- as_data_frame(gr.all_pre, what = "edges");
    df_edges   <- df_edges[ order(df_edges$AVGEXP, df_edges$NCELLS), ];
    gr.all     <- graph_from_data_frame(
                      d = df_edges,
                      vertices = as_data_frame(gr.all_pre, what = "vertices"));
    V(gr.all)$GRP  <- DG[ as.character(attr(V(gr.all), "names")) , "set"];

    pga <- thypg(gr.all, "RAW",
                 thyset.RFK,
                 ncells.max, avgexp.max);

    gr.ncel_pre <- graph_from_data_frame(G[ G$NCELLS > ncells.avg, ],
                                         directed=TRUE);
    df_edges    <- as_data_frame(gr.ncel_pre, what = "edges");
    df_edges    <- df_edges[ order(df_edges$AVGEXP, df_edges$NCELLS), ];
    gr.ncel     <- graph_from_data_frame(
                       d = df_edges,
                       vertices = as_data_frame(gr.ncel_pre, what = "vertices"));
    V(gr.ncel)$GRP <- DG[ as.character(attr(V(gr.ncel), "names")) , "set"];
       
    pgb <- thypg(gr.ncel, paste0("NUM CELS > avg [avg_ncells=",
                                 round(ncells.avg,4),
                                 " / max_ncells=",
                                 round(ncells.max,4),"]"),
                 thyset.RFK,
                 ncells.max, avgexp.max);

    gr.aexp_pre <- graph_from_data_frame(G[ G$AVGEXP > avgexp.avg, ],
                                         directed=TRUE);
    df_edges    <- as_data_frame(gr.aexp_pre, what = "edges");
    df_edges    <- df_edges[ order(df_edges$AVGEXP, df_edges$NCELLS), ];
    gr.aexp     <- graph_from_data_frame(
                       d = df_edges,
                       vertices = as_data_frame(gr.aexp_pre, what = "vertices"));
    V(gr.aexp)$GRP <- DG[ as.character(attr(V(gr.aexp), "names")) , "set"];
       
    pgc <- thypg(gr.aexp, paste0("EXPRESSION > avg [avg_expr=",
                                 round(avgexp.avg,4),
                                 " / max_expr=",
                                 round(avgexp.max,4),"]"),
                 thyset.RFK,
                 ncells.max, avgexp.max);

    gr.cpx_pre <- graph_from_data_frame(G[ G$NCELLS > ncells.avg & G$AVGEXP > avgexp.avg, ],
                                        directed=TRUE);
    df_edges   <- as_data_frame(gr.cpx_pre, what = "edges");
    df_edges   <- df_edges[ order(df_edges$AVGEXP, df_edges$NCELLS), ];
    gr.cpx     <- graph_from_data_frame(
                      d = df_edges,
                      vertices = as_data_frame(gr.cpx_pre, what = "vertices"));
    V(gr.cpx)$GRP <- DG[ as.character(attr(V(gr.cpx), "names")) , "set"];
       
    pgt <- thypg(gr.cpx, paste0("NUM CELS > avg [avg_ncells=",
                                round(ncells.avg,4),
                                " / max_ncells=",
                                round(ncells.max,4),"]\n",
                                "EXPRESSION > avg [avg_expr=",
                                round(avgexp.avg,4),
                                " / max_expr=",
                                round(avgexp.max,4),"]"),
                 thyset.RFK,
                 ncells.max, avgexp.max);

    pgm <- ggarrange(pga + titletheme,
                     pgb + titletheme,
                     pgc + titletheme,
                     pgt + titletheme,
                     nrow = 1, ncol = 4,
                     common.legend = TRUE,
                     legend = "right");

    ggsave(plot=pgm,
           filename=paste0(WD,"/coexpression_graphs/Reddien_coexpression_graph_",
                           THYGNS,"_by_region_filters.pdf"),
           device = "pdf",
           width  = 120 * S.scale, height = 50 * S.scale,
           units= "in", dpi = 300,
           bg = "white", limitsize = FALSE);

    pgg <- pgt + facet_wrap(. ~ SET);

    ggsave(plot=pgg,
           filename=paste0(WD,"/coexpression_graphs/Reddien_coexpression_graph_",
                           THYGNS,"_by_region_NCels+AvgExpr.pdf"),
           device = "pdf",
           width  = 80 * S.scale, height = 50 * S.scale,
           units= "in", dpi = 300,
           bg = "white", limitsize = FALSE);


    # Graphs by CELL TYPE (Reddien) 
    DO <- DT.Red[[THYGNS]] %>%
          select(c(thyset.RFK$geneid, "label_short"));
    TC <- nrow(DO);

    CO.Red <- list();
    CN.Red <- list();
    for (i in Red.CELL.ord) {
      if (i == "Other") next;
      # cat("# ", i, "\n");
      CO.Red[[i]] <- matrix(0L, nrow=S, ncol=S,
                            dimnames=list(thyset.RFK$geneid, thyset.RFK$geneid));
      CN.Red[[i]] <- 0;
    };
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
          if (j == S & G2 > 0) {
            CO.Red[[RG]][G2,G2] <- CO.Red[[RG]][G2,G2] + 1;
          };
        }; # j
      }; # i
    }; # d
    CA.Red <- list();
    for (d in Red.CELL.ord) {
      if (d == "Other") next;
      CA.Red[[d]] <- matrix(0L, nrow=S, ncol=S, 
                            dimnames=list(thyset.RFK$geneid, thyset.RFK$geneid));
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

    write.table(CA.Red,
                file=paste0(WD,"/coexpression_graphs/Reddien_coexpression_matrix_",
                            THYGNS,"_bycelltype.tbl"),
                quote=TRUE, na="NA", col.names=TRUE, row.names=TRUE, sep="\t");

    G <- data.frame();
    for (d in Red.CELL.ord) {
      if (d == "Other") next;
      for (i in 1:S) {
        G1 <- thyset.RFK$geneid[[i]];
        for (j in 1:S) {
          G2 <- thyset.RFK$geneid[[j]];
          if (i < j) {
            Ncell <- CA.Red[[d]][G1,G2];
            Nexpr <- CA.Red[[d]][G2,G1];
            # cat("Ncell: ",Ncell,"  Nexpr: ",Nexpr,"\n");
            if (!is.null(Ncell) & Ncell > 0) {
              if (nrow(G) > 0) {
                G <- rbind(G, c(G1, G2, Ncell, Nexpr, d));
              } else {
                G <- rbind(G, c(G1, G2, Ncell, Nexpr, d));
                colnames(G) <- gcolnames;
              };
            };
          };
        };
      };
    };

    G$NCELLS <- as.numeric(G$NCELLS);
    G$AVGEXP <- as.numeric(G$AVGEXP);

    write.table(G,
                file=paste0(WD,"/coexpression_graphs/Reddien_coexpression_matrix_",
                            THYGNS,"_bycelltype_alledges.tbl"),
                quote=TRUE, na="NA", col.names=TRUE, row.names=FALSE, sep="\t");

    cat("#--> CELL TYPE REDDIEN\n");

    ncells.avg <- mean(G$NCELLS, na.rm=TRUE);
    avgexp.avg <- mean(G$AVGEXP, na.rm=TRUE);
    ncells.max <- roundUp(max(G$NCELLS, na.rm=TRUE));
    avgexp.max <- roundUp(max(G$AVGEXP, na.rm=TRUE));

    gr.all_pre <- graph_from_data_frame(G, directed=TRUE);
    df_edges   <- as_data_frame(gr.all_pre, what = "edges");
    df_edges   <- df_edges[ order(df_edges$AVGEXP, df_edges$NCELLS), ];
    gr.all     <- graph_from_data_frame(
                      d = df_edges,
                      vertices = as_data_frame(gr.all_pre, what = "vertices"));
    V(gr.all)$GRP  <- DG[ as.character(attr(V(gr.all), "names")) , "set"];

    pga <- thypg(gr.all, "RAW",
                 thyset.RFK,
                 ncells.max, avgexp.max);

    gr.ncel_pre <- graph_from_data_frame(G[ G$NCELLS > ncells.avg, ],
                                     directed=TRUE);
    df_edges    <- as_data_frame(gr.ncel_pre, what = "edges");
    df_edges    <- df_edges[ order(df_edges$AVGEXP, df_edges$NCELLS), ];
    gr.ncel     <- graph_from_data_frame(
                       d = df_edges,
                       vertices = as_data_frame(gr.ncel_pre, what = "vertices"));
    V(gr.ncel)$GRP <- DG[ as.character(attr(V(gr.ncel), "names")) , "set"];
       
    pgb <- thypg(gr.ncel, paste0("NUM CELS > avg [avg_ncells=",
                                 round(ncells.avg,4),
                                 " / max_ncells=",
                                 round(ncells.max,4),"]"),
                 thyset.RFK,
                 ncells.max, avgexp.max);

    gr.aexp_pre <- graph_from_data_frame(G[ G$AVGEXP > avgexp.avg, ],
                                         directed=TRUE);
    df_edges    <- as_data_frame(gr.aexp_pre, what = "edges");
    df_edges    <- df_edges[ order(df_edges$AVGEXP, df_edges$NCELLS), ];
    gr.aexp     <- graph_from_data_frame(
                       d = df_edges,
                       vertices = as_data_frame(gr.aexp_pre, what = "vertices"));
    V(gr.aexp)$GRP <- DG[ as.character(attr(V(gr.aexp), "names")) , "set"];
       
    pgc <- thypg(gr.aexp, paste0("EXPRESSION > avg [avg_expr=",
                                 round(avgexp.avg,4),
                                 " / max_expr=",
                                 round(avgexp.max,4),"]"),
                 thyset.RFK,
                 ncells.max, avgexp.max);
       
    gr.cpx_pre <- graph_from_data_frame(G[ G$NCELLS > ncells.avg &
                                        G$AVGEXP > avgexp.avg, ],
                                     directed=TRUE);
    df_edges   <- as_data_frame(gr.cpx_pre, what = "edges");
    df_edges   <- df_edges[ order(df_edges$AVGEXP, df_edges$NCELLS), ];
    gr.cpx     <- graph_from_data_frame(
                      d = df_edges,
                      vertices = as_data_frame(gr.cpx_pre, what = "vertices"));
    V(gr.cpx)$GRP <- DG[ as.character(attr(V(gr.cpx), "names")) , "set"];
       
    pgt <- thypg(gr.cpx, paste0("NUM CELS > avg [avg_ncells=",
                                round(ncells.avg,4),
                                " / max_ncells=",
                                round(ncells.max,4),"]\n",
                                "EXPRESSION > avg [avg_expr=",
                                round(avgexp.avg,4),
                                " / max_expr=",
                                round(avgexp.max,4),"]"),
                 thyset.RFK,
                 ncells.max, avgexp.max);

    # pgm <- ggarrange(pga + theme(plot.title = element_text(hjust = 0.5, size = 40, colour = "blue", 
    #                                                        family = "Verdana", face = "bold")),
    #                  pgb + theme(plot.title = element_text(hjust = 0.5, size = 40, colour = "blue", 
    #                                                        family = "Verdana", face = "bold")),
    #                  pgc + theme(plot.title = element_text(hjust = 0.5, size = 40, colour = "blue", 
    #                                                        family = "Verdana", face = "bold")),
    #                  pgt + theme(plot.title = element_text(hjust = 0.5, size = 40, colour = "blue", 
    #                                                        family = "Verdana", face = "bold")),
    pgm <- ggarrange(pga + titletheme,
                     pgb + titletheme,
                     pgc + titletheme,
                     pgt + titletheme,
                     nrow = 1, ncol = 4,
                     common.legend = TRUE,
                     legend = "right");

    ggsave(plot=pgm,
           filename=paste0(WD,"/coexpression_graphs/Reddien_coexpression_graph_",
                           THYGNS,"_by_celltype_filters.pdf"),
           device = "pdf",
           width  = 120 * S.scale, height = 50 * S.scale,
           units= "in", dpi = 300,
           bg = "white", limitsize = FALSE);

    pgg <- pgt + facet_wrap(. ~ SET);

    ggsave(plot=pgg,
           filename=paste0(WD,"/coexpression_graphs/Reddien_coexpression_graph_",
                           THYGNS,"_by_celltype_NCels+AvgExpr.pdf"),
           device = "pdf",
           width  = 100 * S.scale, height = 90 * S.scale,
           units= "in", dpi = 300,
           bg = "white", limitsize = FALSE);


    # Graphs by CELL TYPE (Rajewsky) 
    DO <- DT.Raj[[THYGNS]] %>%
          select(c(thyset.RFK$geneid, "label_short"));
    TC <- nrow(DO);
                  
    CO.Raj <- list();
    CN.Raj <- list();
    for (i in Raj.CELL.ord2) {
      if (i == "Other") next;
      #cat("#a ", i, "\n");
      CO.Raj[[i]] <- matrix(0L, nrow=S, ncol=S,
                            dimnames=list(thyset.RFK$geneid, thyset.RFK$geneid));
      CN.Raj[[i]] <- 0;
    };
    for (d in 1:nrow(DO)) {
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
          if (j == S & G2 > 0) {
            CO.Raj[[RG]][G2,G2] <- CO.Raj[[RG]][G2,G2] + 1;
          };
        }; # j
      }; # i
    }; # d
    CA.Raj <- list();
    for (d in Raj.CELL.ord2) {
      if (d == "Other") next;
      #cat("#c ", d, "\n");
      CA.Raj[[d]] <- matrix(0L, nrow=S, ncol=S, 
                            dimnames=list(thyset.RFK$geneid, thyset.RFK$geneid));
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

    write.table(CA.Raj,
                file=paste0(WD,"/coexpression_graphs/Rajewsky_coexpression_matrix_",
                            THYGNS,"_bycelltype.tbl"),
                quote=TRUE, na="NA", col.names=TRUE, row.names=TRUE, sep="\t");

    G <- data.frame();
    for (d in Raj.CELL.ord2) {
      if (d == "Other") next;
      for (i in 1:S) {
        G1 <- thyset.RFK$geneid[[i]];
        for (j in 1:S) {
          G2 <- thyset.RFK$geneid[[j]];
          if (i<j) {
            Ncell <- CA.Raj[[d]][G1,G2];
            Nexpr <- CA.Raj[[d]][G2,G1];
            # cat("Ncell: ",Ncell,"  Nexpr: ",Nexpr,"\n");
            if (!is.null(Ncell) & Ncell > 0) {
              if (nrow(G) > 0) {
                G <- rbind(G, c(G1, G2, Ncell, Nexpr, d));
              } else {
                G <- rbind(G, c(G1, G2, Ncell, Nexpr, d));
                colnames(G) <- gcolnames;
              };
            };
          };
        };
      };
    };

    G$NCELLS <- as.numeric(G$NCELLS);
    G$AVGEXP <- as.numeric(G$AVGEXP);

    write.table(G,
                file=paste0(WD,"/coexpression_graphs/Rajewsky_coexpression_matrix_",
                            THYGNS,"_bycelltype_alledges.tbl"),
                quote=TRUE, na="NA", col.names=TRUE, row.names=FALSE, sep="\t");

    cat("#--> CELL TYPE RAJEWSKY\n");

    ncells.avg <- mean(G$NCELLS, na.rm=TRUE);
    avgexp.avg <- mean(G$AVGEXP, na.rm=TRUE);
    ncells.max <- roundUp(max(G$NCELLS, na.rm=TRUE));
    avgexp.max <- roundUp(max(G$AVGEXP, na.rm=TRUE));

    gr.all_pre <- graph_from_data_frame(G, directed=TRUE);
    df_edges   <- as_data_frame(gr.all_pre, what = "edges");
    df_edges   <- df_edges[ order(df_edges$AVGEXP, df_edges$NCELLS), ];
    gr.all     <- graph_from_data_frame(
                      d = df_edges,
                      vertices = as_data_frame(gr.all_pre, what = "vertices"));
    V(gr.all)$GRP  <- DG[ as.character(attr(V(gr.all), "names")) , "set"];

    pga <- thypg(gr.all, "RAW",
                 thyset.RFK,
                 ncells.max, avgexp.max);

    gr.ncel_pre <- graph_from_data_frame(G[ G$NCELLS > ncells.avg, ],
                                         directed=TRUE);
    df_edges    <- as_data_frame(gr.ncel_pre, what = "edges");
    df_edges    <- df_edges[ order(df_edges$AVGEXP, df_edges$NCELLS), ];
    gr.ncel     <- graph_from_data_frame(
                       d = df_edges,
                       vertices = as_data_frame(gr.ncel_pre, what = "vertices"));
    V(gr.ncel)$GRP <- DG[ as.character(attr(V(gr.ncel), "names")) , "set"];
       
    pgb <- thypg(gr.ncel, paste0("NUM CELS > avg [avg_ncells=",
                                 round(ncells.avg,4),
                                 " / max_ncells=",
                                 round(ncells.max,4),"]"),
                 thyset.RFK,
                 ncells.max, avgexp.max);

    gr.aexp_pre <- graph_from_data_frame(G[ G$AVGEXP > avgexp.avg, ],
                                         directed=TRUE);
    df_edges    <- as_data_frame(gr.aexp_pre, what = "edges");
    df_edges    <- df_edges[ order(df_edges$AVGEXP, df_edges$NCELLS), ];
    gr.aexp     <- graph_from_data_frame(
                       d = df_edges,
                       vertices = as_data_frame(gr.aexp_pre, what = "vertices"));
    V(gr.aexp)$GRP <- DG[ as.character(attr(V(gr.aexp), "names")) , "set"];
       
    pgc <- thypg(gr.aexp, paste0("EXPRESSION > avg [avg_expr=",
                                 round(avgexp.avg,4),
                                 " / max_expr=",
                                 round(avgexp.max,4),"]"),
                 thyset.RFK,
                 ncells.max, avgexp.max);
       
    gr.cpx_pre <- graph_from_data_frame(G[ G$NCELLS > ncells.avg &
                                           G$AVGEXP > avgexp.avg, ],
                                        directed=TRUE);
    df_edges   <- as_data_frame(gr.cpx_pre, what = "edges");
    df_edges   <- df_edges[ order(df_edges$AVGEXP, df_edges$NCELLS), ];
    gr.cpx     <- graph_from_data_frame(
                      d = df_edges,
                      vertices = as_data_frame(gr.cpx_pre, what = "vertices"));
    V(gr.cpx)$GRP <- DG[ as.character(attr(V(gr.cpx), "names")) , "set"];
        
    pgt <- thypg(gr.cpx, paste0("NUM CELS > avg [avg_ncells=",
                                round(ncells.avg,4),
                                " / max_ncells=",round(ncells.max,4),"]\n",
                                "EXPRESSION > avg [avg_expr=",
                                round(avgexp.avg,4),
                                " / max_expr=",
                                round(avgexp.max,4),"]"),
                 thyset.RFK,
                 ncells.max, avgexp.max);

    # pgm <- ggarrange(pga + theme(plot.title = element_text(hjust = 0.5, size = 40, colour = "blue", 
    #                                                        family = "Verdana", face = "bold")),
    #                  pgb + theme(plot.title = element_text(hjust = 0.5, size = 40, colour = "blue", 
    #                                                        family = "Verdana", face = "bold")),
    #                  pgc + theme(plot.title = element_text(hjust = 0.5, size = 40, colour = "blue", 
    #                                                        family = "Verdana", face = "bold")),
    #                  pgt + theme(plot.title = element_text(hjust = 0.5, size = 40, colour = "blue", 
    #                                                        family = "Verdana", face = "bold")),
    pgm <- ggarrange(pga + titletheme,
                     pgb + titletheme,
                     pgc + titletheme,
                     pgt + titletheme,
                     nrow = 1, ncol = 4,
                     common.legend = TRUE,
                     legend = "right");

    ggsave(plot=pgm,
           filename=paste0(WD,"/coexpression_graphs/Rajewsky_coexpression_graph_",
                           THYGNS,"_by_celltype_filters.pdf"),
           device = "pdf",
           width  = 120 * S.scale, height = 50 * S.scale,
           units= "in", dpi = 300,
           bg = "white", limitsize = FALSE);

    pgg <- pgt + facet_wrap(. ~ SET);

    ggsave(plot=pgg,
           filename=paste0(WD,"/coexpression_graphs/Rajewsky_coexpression_graph_",
                           THYGNS,"_by_celltype_NCels+AvgExpr.pdf"),
           device = "pdf",
           width  = 100 * S.scale, height = 90 * S.scale,
           units= "in", dpi = 300,
           bg = "white", limitsize = FALSE);
   
} # make_all_coexpression_graphs
