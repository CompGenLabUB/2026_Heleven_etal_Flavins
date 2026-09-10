make_all_heatmaps <- function(THYGNS) {
        
    thyset.RFK <- DG[ DG[[THYGNS]] == 1, c("geneid","symbol","name","set") ];
    SETS <- unique(thyset.RFK$set);

    cat("###\n### WORKING on ",THYGNS,"\n");
         
    # Heatmaps by REGION
    MT.Red <- DT.Red[[THYGNS]] %>%
              select(c(thyset.RFK$geneid,"region_label")) %>%  # select(c(9:49,8)) %>%
              pivot_longer(1:length(thyset.RFK$geneid),names_to="gene",values_to="expression") %>%
              group_by(region_label, gene) %>% 
              summarize(count=sum(expression>0),
                        mean=mean(expression,na.rm=TRUE),
                        meangt=mean(expression[expression>0], na.rm=TRUE));
    MT.Red.count  <- MT.Red %>%
                     select(region_label,gene,count) %>%
                     pivot_wider(names_from = region_label, values_from = count ) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    MT.Red.mean   <- MT.Red %>%
                     select(region_label,gene,mean) %>%
                     pivot_wider(names_from = region_label, values_from = mean  ) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    MT.Red.maxmean <- max(MT.Red$mean,     na.rm=TRUE);
    MT.Red.meangt <- MT.Red %>%
                     select(region_label,gene,meangt) %>%
                     pivot_wider(names_from = region_label, values_from = meangt) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    MT.Red.maxmeangt <- max(MT.Red$meangt, na.rm=TRUE);
    DT.Red.count  <- as.data.frame(MT.Red.count);
    rownames(DT.Red.count)  <- DT.Red.count$gene;
    DT.Red.count  <- DT.Red.count[ thyset.RFK$geneid, ];
    DT.Red.mean   <- as.data.frame(MT.Red.mean);
    rownames(DT.Red.mean)   <- DT.Red.mean$gene;
    DT.Red.mean   <- DT.Red.mean[ thyset.RFK$geneid, ];
    DT.Red.meangt <- as.data.frame(MT.Red.meangt);
    rownames(DT.Red.meangt) <- DT.Red.meangt$gene;
    DT.Red.meangt <- DT.Red.meangt[ thyset.RFK$geneid, ];
    DT.Red.count  <- as.matrix(subset(DT.Red.count,  select = -gene));
    DT.Red.mean   <- as.matrix(subset(DT.Red.mean,   select = -gene));
    DT.Red.meangt <- as.matrix(subset(DT.Red.meangt, select = -gene));
    DT.Red.meangtlg2 <- log2(DT.Red.meangt);
    DT.Red.meangtlg2[is.infinite(DT.Red.meangtlg2)] <- NA

    DT.Red.rows <- nrow(DT.Red.count);
    DT.Red.cols <- ncol(DT.Red.count);
    DT.col_fun <- colorRamp2(seq(0,MT.Red.maxmean,MT.Red.maxmean/24), HMPcols);
    htmp <- Heatmap(
            DT.Red.mean,
            na_col = "#CCCCCC",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            cluster_row_slices = FALSE,
            column_order = HTMP.order,
            row_order = rownames(DT.Red.mean), # thyset.RFK$geneid, # 
            row_labels = sprintf("%25s", thyset.RFK[ rownames(DT.Red.mean), "symbol" ]),
            row_names_side = "left",
            row_split = factor(thyset.RFK[ rownames(DT.Red.mean), "set" ], levels=SETS),
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            col = DT.col_fun,
            column_names_rot = 45,
            column_names_gp = gpar(fontsize = 15, fontface = "bold"),
            heatmap_width  = unit(1.25 * DT.Red.cols + 12,  "cm"),
            heatmap_height = unit(0.5 * DT.Red.rows + 8, "cm"),
            cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%d", as.numeric(as.matrix(DT.Red.count)[i, j])),
                                    x, y,
                                    gp = gpar(fontsize = 10,
                                              col = ifelse(
                                                      is.numeric(as.matrix(DT.Red.mean)[i, j]) &
                                                      as.matrix(DT.Red.mean)[i, j] < 0.5*MT.Red.maxmean,
                                                      "white", "black")
                                              ))
                       },
            row_title = NULL,
            heatmap_legend_param = list(
                col_fun = DT.col_fun, title = "Expression\nLevel",
                legend_width = unit(6, "cm"), 
                legend_height = unit(1, "cm"), 
                direction = "horizontal",
                title_position = "leftcenter"),
            right_annotation = rowAnnotation(
              foo = anno_block(
                  labels = thyset.RFK$set, # SETS,
                  panel_fun = function(index, levels) {
                      # txt = gsub(" ", "\n", thyset.RFK[ rownames(DT.Red.mean)[index], "set" ]);
                      txt = thyset.RFK[ rownames(DT.Red.mean)[index], "set" ];
                      grid.text(txt, 0, 0.5, rot = 0, just = "left",
                                gp = gpar(col = "black", 
                                          fontsize = 14, fontface = "bold"))
                  },
                  width = unit(6, "cm")
              ))
        );
    pdf(paste0(WD,"/heatmaps/Reddien_heatmap_",THYGNS,"_by_region.pdf"),
        paper = "special",
        width  = unit((1.25 * DT.Red.cols + 12) / 2.54, "cm"), # 1in~2.54cm
        height = unit((0.55 * DT.Red.rows + 10) / 2.54, "cm"),
        bg = "white");
    draw(htmp, heatmap_legend_side="top");
    dev.off();
    cat("#-> Reddien heatmap by region\n");

    DT.col_fun2 <- colorRamp2(seq(1,MT.Red.maxmeangt,(MT.Red.maxmeangt-1)/24), HMPcols);
    htmp <- Heatmap(
            DT.Red.meangt,
            na_col = "#CCCCCC",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            cluster_row_slices = FALSE,
            column_order = HTMP.order,
            row_order = rownames(DT.Red.meangt), # thyset.RFK$geneid, # 
            row_labels = sprintf("%25s", thyset.RFK[ rownames(DT.Red.mean), "symbol" ]),
            row_names_side = "left",
            row_split = factor(thyset.RFK[rownames(DT.Red.meangt), "set" ], levels=SETS),
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            col = DT.col_fun2,
            column_names_rot = 45,
            column_names_gp = gpar(fontsize = 15, fontface = "bold"),
            heatmap_width  = unit(1.25 * DT.Red.cols + 12,  "cm"),
            heatmap_height = unit(0.5 * DT.Red.rows + 8, "cm"),
            cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%d", as.numeric(as.matrix(DT.Red.count)[i, j])),
                                    x, y,
                                    gp = gpar(fontsize = 10,
                                              col = ifelse(
                                                      is.numeric(as.matrix(DT.Red.meangt)[i, j]) &
                                                      as.matrix(DT.Red.meangt)[i, j] < 0.5*MT.Red.maxmeangt,
                                                      "white", "black")
                                              ))
                       },
            row_title = NULL,
            heatmap_legend_param = list(
                col_fun = DT.col_fun2, title = "Expression\nLevel",
                legend_width = unit(6, "cm"), 
                legend_height = unit(1, "cm"), 
                direction = "horizontal",
                title_position = "leftcenter"),
            right_annotation = rowAnnotation(
              foo = anno_block(
                  labels = thyset.RFK$set, # SETS,
                  panel_fun = function(index, levels) {
                      # txt = gsub(" ", "\n", thyset.RFK[rownames(DT.Red.meangt)[index], "set"]);
                      txt = thyset.RFK[rownames(DT.Red.meangt)[index], "set"];
                      grid.text(txt, 0, 0.5, rot = 0, just = "left",
                                gp = gpar(col = "black", 
                                          fontsize = 14, fontface = "bold"))
                  },
                  width = unit(6, "cm")
              ))
        );
    pdf(paste0(WD,"/heatmaps/Reddien_heatmap_",THYGNS,"_by_region_meangt0.pdf"),
        paper = "special",
        width  = unit((1.25 * DT.Red.cols + 12) / 2.54, "cm"), # 1in~2.54cm
        height = unit((0.55 * DT.Red.rows + 10) / 2.54, "cm"),
        bg = "white");
    draw(htmp, heatmap_legend_side="top");
    dev.off();
    cat("#-> Reddien heatmap by region (mean gt0)\n");

            
    # Heatmaps by LABEL_SHORT
    ML.Red <- DT.Red[[THYGNS]] %>%
              select(c(thyset.RFK$geneid,"label_short")) %>%  # select(c(9:49,52)) %>%
              pivot_longer(1:length(thyset.RFK$geneid),names_to="gene",values_to="expression") %>%
              group_by(label_short, gene) %>% 
              summarize(count=sum(expression>0),
                        mean=mean(expression,na.rm=TRUE),
                        meangt=mean(expression[expression>0], na.rm=TRUE));
    ML.Red.count  <- ML.Red %>%
                     select(label_short,gene,count) %>%
                     pivot_wider(names_from = label_short, values_from = count ) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    ML.Red.mean   <- ML.Red %>%
                     select(label_short,gene,mean) %>%
                     pivot_wider(names_from = label_short, values_from = mean  ) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    ML.Red.maxmean <- max(ML.Red$mean,     na.rm=TRUE);
    ML.Red.meangt <- ML.Red %>%
                     select(label_short,gene,meangt) %>%
                     pivot_wider(names_from = label_short, values_from = meangt) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    ML.Red.maxmeangt <- max(ML.Red$meangt, na.rm=TRUE);
    DL.Red.count  <- as.data.frame(ML.Red.count);
    rownames(DL.Red.count)  <- DL.Red.count$gene;
    DL.Red.count  <- DL.Red.count[ thyset.RFK$geneid, ];
    DL.Red.mean   <- as.data.frame(ML.Red.mean);
    rownames(DL.Red.mean)   <- DL.Red.mean$gene;
    DL.Red.mean   <- DL.Red.mean[ thyset.RFK$geneid, ];
    DL.Red.meangt <- as.data.frame(ML.Red.meangt);
    rownames(DL.Red.meangt) <- DL.Red.meangt$gene;
    DL.Red.meangt <- DL.Red.meangt[ thyset.RFK$geneid, ];
    DL.Red.count  <- as.matrix(subset(DL.Red.count,  select = -gene));
    DL.Red.mean   <- as.matrix(subset(DL.Red.mean,   select = -gene));
    DL.Red.meangt <- as.matrix(subset(DL.Red.meangt, select = -gene));
    DL.Red.meangtlg2 <- log2(DL.Red.meangt);
    DL.Red.meangtlg2[is.infinite(DL.Red.meangtlg2)] <- NA

    DL.Red.rows <- nrow(DL.Red.count);
    DL.Red.cols <- ncol(DL.Red.count);
    DL.Red.col_fun <- colorRamp2(seq(0,ML.Red.maxmean,ML.Red.maxmean/24), HMPcols);
    htmp <- Heatmap(
            DL.Red.mean,
            na_col = "#CCCCCC",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            cluster_row_slices = FALSE,
            column_order = colnames(DL.Red.mean),
            row_order = rownames(DL.Red.mean), # thyset.RFK$geneid, # 
            row_labels = sprintf("%25s", thyset.RFK[ rownames(DL.Red.mean), "symbol" ]),
            row_names_side = "left",
            row_split = factor(thyset.RFK[ rownames(DL.Red.mean), "set" ], levels=SETS),
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            col = DL.col_fun,
            column_names_rot = 45,
            column_names_gp = gpar(fontsize = 15, fontface = "bold"),
            heatmap_width  = unit(1.25 * DL.Red.cols + 12,  "cm"),
            heatmap_height = unit(0.5 * DL.Red.rows + 8, "cm"),
            cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%d", as.numeric(as.matrix(DL.Red.count)[i, j])),
                                    x, y,
                                    gp = gpar(fontsize = 10,
                                              col = ifelse(
                                                      is.numeric(as.matrix(DL.Red.mean)[i, j]) &
                                                      as.matrix(DL.Red.mean)[i, j] < 0.5*ML.Red.maxmean,
                                                      "white", "black")
                                              ))
                       },
            row_title = NULL,
            heatmap_legend_param = list(
                col_fun = DL.col_fun, title = "Expression\nLevel",
                legend_width = unit(6, "cm"), 
                legend_height = unit(1, "cm"), 
                direction = "horizontal",
                title_position = "leftcenter"),
            right_annotation = rowAnnotation(
              foo = anno_block(
                  labels = thyset.RFK$set, # SETS,
                  panel_fun = function(index, levels) {
                      # txt = gsub(" ", "\n", thyset.RFK[ rownames(DL.Red.mean)[index], "set" ]);
                      txt = thyset.RFK[ rownames(DL.Red.mean)[index], "set" ];
                      grid.text(txt, 0, 0.5, rot = 0, just = "left",
                                gp = gpar(col = "black", 
                                          fontsize = 14, fontface = "bold"))
                  },
                  width = unit(6, "cm")
              ))
        );
    pdf(paste0(WD,"/heatmaps/Reddien_heatmap_",THYGNS,"_by_celltype.pdf"),
        paper = "special",
        width  = unit((1.25 * DL.Red.cols + 12) / 2.54, "cm"), # 1in~2.54cm
        height = unit((0.55 * DL.Red.rows + 10) / 2.54, "cm"),
        bg = "white");
    draw(htmp, heatmap_legend_side="top");
    dev.off();
    cat("#-> Reddien heatmap by label_short\n");

    DL.Red.col_fun2 <- colorRamp2(seq(1,ML.Red.maxmeangt,(ML.Red.maxmeangt-1)/24), HMPcols);
    htmp <- Heatmap(
            DL.Red.meangt,
            na_col = "#CCCCCC",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            cluster_row_slices = FALSE,
            column_order = colnames(DL.Red.meangt),
            row_order = rownames(DL.Red.meangt), # thyset.RFK$geneid, # 
            row_labels = sprintf("%25s", thyset.RFK[ rownames(DL.Red.meangt), "symbol" ]),
            row_names_side = "left",
            row_split = factor(thyset.RFK[ rownames(DL.Red.meangt), "set" ], levels=SETS),
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            col = DL.Red.col_fun2,
            column_names_rot = 45,
            column_names_gp = gpar(fontsize = 15, fontface = "bold"),
            heatmap_width  = unit(1.25 * DL.Red.cols + 12,  "cm"),
            heatmap_height = unit(0.5 * DL.Red.rows + 8, "cm"),
            cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%d", as.numeric(as.matrix(DL.Red.count)[i, j])),
                                    x, y,
                                    gp = gpar(fontsize = 10,
                                              col = ifelse(
                                                      is.numeric(as.matrix(DL.Red.meangt)[i, j]) &
                                                      as.matrix(DL.Red.meangt)[i, j] < 0.5*ML.Red.maxmeangt,
                                                      "white", "black")
                                              ))
                       },
            row_title = NULL,
            heatmap_legend_param = list(
                col_fun = DL.Red.col_fun2, title = "Expression\nLevel",
                legend_width = unit(6, "cm"), 
                legend_height = unit(1, "cm"), 
                direction = "horizontal",
                title_position = "leftcenter"),
            right_annotation = rowAnnotation(
              foo = anno_block(
                  labels = thyset.RFK$set, # SETS,
                  panel_fun = function(index, levels) {
                      # txt = gsub(" ", "\n", thyset.RFK[ rownames(DL.Red.meangt)[index], "set" ]);
                      txt = thyset.RFK[ rownames(DL.Red.meangt)[index], "set" ];
                      grid.text(txt, 0, 0.5, rot = 0, just = "left",
                                gp = gpar(col = "black", 
                                          fontsize = 14, fontface = "bold"))
                  },
                  width = unit(6, "cm")
              ))
        );
    pdf(paste0(WD,"/heatmaps/Reddien_heatmap_",THYGNS,"_by_celltype_meangt0.pdf"),
        paper = "special",
        width  = unit((1.25 * DL.Red.cols + 12) / 2.54, "cm"), # 1in~2.54cm
        height = unit((0.55 * DL.Red.rows + 10) / 2.54, "cm"),
        bg = "white");
    draw(htmp, heatmap_legend_side="top");
    dev.off();
    cat("#-> Reddien heatmap by label_short (mean gt0)\n");


    ML.Raj <- DT.Raj[[THYGNS]] %>%
              select(c(thyset.RFK$geneid,"label_short")) %>%  # select(c(9:49,51)) %>%
              pivot_longer(1:length(thyset.RFK$geneid),names_to="gene",values_to="expression") %>%
              group_by(label_short, gene) %>% 
              summarize(count=sum(expression>0),
                        mean=mean(expression,na.rm=TRUE),
                        meangt=mean(expression[expression>0], na.rm=TRUE));
    ML.Raj.count  <- ML.Raj %>%
                     select(label_short,gene,count) %>%
                     pivot_wider(names_from = label_short, values_from = count ) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    ML.Raj.mean   <- ML.Raj %>%
                     select(label_short,gene,mean) %>%
                     pivot_wider(names_from = label_short, values_from = mean  ) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    ML.Raj.maxmean <- max(ML.Raj$mean,     na.rm=TRUE);
    ML.Raj.meangt <- ML.Raj %>%
                     select(label_short,gene,meangt) %>%
                     pivot_wider(names_from = label_short, values_from = meangt) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    ML.Raj.maxmeangt <- max(ML.Raj$meangt, na.rm=TRUE);
    DL.Raj.count  <- as.data.frame(ML.Raj.count);
    rownames(DL.Raj.count)  <- DL.Raj.count$gene;
    DL.Raj.count  <- DL.Raj.count[ thyset.RFK$geneid, ];
    DL.Raj.mean   <- as.data.frame(ML.Raj.mean);
    rownames(DL.Raj.mean)   <- DL.Raj.mean$gene;
    DL.Raj.mean   <- DL.Raj.mean[ thyset.RFK$geneid, ];
    DL.Raj.meangt <- as.data.frame(ML.Raj.meangt);
    rownames(DL.Raj.meangt) <- DL.Raj.meangt$gene;
    DL.Raj.meangt <- DL.Raj.meangt[ thyset.RFK$geneid, ];
    DL.Raj.count  <- as.matrix(subset(DL.Raj.count,  select = -gene));
    DL.Raj.mean   <- as.matrix(subset(DL.Raj.mean,   select = -gene));
    DL.Raj.meangt <- as.matrix(subset(DL.Raj.meangt, select = -gene));
    DL.Raj.meangtlg2 <- log2(DL.Raj.meangt);
    DL.Raj.meangtlg2[is.infinite(DL.Raj.meangtlg2)] <- NA

    DL.Raj.rows <- nrow(DL.Raj.count);
    DL.Raj.cols <- ncol(DL.Raj.count);
    DL.Raj.col_fun <- colorRamp2(seq(0,ML.Raj.maxmean,ML.Raj.maxmean/24), HMPcols);
    htmp <- Heatmap(
            DL.Raj.mean,
            na_col = "#CCCCCC",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            cluster_row_slices = FALSE,
            column_order = colnames(DL.Raj.mean),
            row_order = rownames(DL.Raj.mean), # thyset.RFK$geneid, # 
            row_labels = sprintf("%25s", thyset.RFK[ rownames(DL.Raj.mean), "symbol" ]),
            row_names_side = "left",
            row_split = factor(thyset.RFK[ rownames(DL.Raj.mean), "set" ], levels=SETS),
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            col = DL.Raj.col_fun,
            column_names_rot = 70,
            column_names_gp = gpar(fontsize = 15, fontface = "bold"),
            heatmap_width  = unit(1.25 * DL.Raj.cols + 12,  "cm"),
            heatmap_height = unit(0.5 * DL.Raj.rows + 8, "cm"),
            cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%d", as.numeric(as.matrix(DL.Raj.count)[i, j])),
                                    x, y,
                                    gp = gpar(fontsize = 10,
                                              col = ifelse(
                                                      is.numeric(as.matrix(DL.Raj.mean)[i, j]) &
                                                      as.matrix(DL.Raj.mean)[i, j] < 0.5*ML.Raj.maxmean,
                                                      "white", "black")
                                              ))
                       },
            row_title = NULL,
            heatmap_legend_param = list(
                col_fun = DL.Raj.col_fun, title = "Expression\nLevel",
                legend_width = unit(6, "cm"), 
                legend_height = unit(1, "cm"), 
                direction = "horizontal",
                title_position = "leftcenter"),
            right_annotation = rowAnnotation(
              foo = anno_block(
                  labels = thyset.RFK$set, # SETS,
                  panel_fun = function(index, levels) {
                      # txt = gsub(" ", "\n", thyset.RFK[ rownames(DL.Raj.mean)[index], "set" ]);
                      txt = thyset.RFK[ rownames(DL.Raj.mean)[index], "set" ];
                      grid.text(txt, 0, 0.5, rot = 0, just = "left",
                                gp = gpar(col = "black", 
                                          fontsize = 14, fontface = "bold"))
                  },
                  width = unit(6, "cm")
              ))
        );
    pdf(paste0(WD,"/heatmaps/Rajewsky_heatmap_",THYGNS,"_by_celltype.pdf"),
        paper = "special",
        width  = unit((1.25 * DL.Raj.cols + 12) / 2.54, "cm"), # 1in~2.54cm
        height = unit((0.55 * DL.Raj.rows + 15) / 2.54, "cm"),
        bg = "white");
    draw(htmp, heatmap_legend_side="top");
    dev.off();
    cat("#-> Rajewsky heatmap by label_short\n");

    DL.Raj.col_fun2 <- colorRamp2(seq(1,ML.Raj.maxmeangt,(ML.Raj.maxmeangt - 1)/24), HMPcols);
    htmp <- Heatmap(
            DL.Raj.meangt,
            na_col = "#CCCCCC",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            cluster_row_slices = FALSE,
            column_order = colnames(DL.Raj.meangt),
            row_order = rownames(DL.Raj.meangt), # thyset.RFK$geneid, # 
            row_labels = sprintf("%25s", thyset.RFK[ rownames(DL.Raj.meangt), "symbol" ]),
            row_names_side = "left",
            row_split = factor(thyset.RFK[ rownames(DL.Raj.meangt), "set" ], levels=SETS),
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            col = DL.Raj.col_fun2,
            column_names_rot = 70,
            column_names_gp = gpar(fontsize = 15, fontface = "bold"),
            heatmap_width  = unit(1.25 * DL.Raj.cols + 12,  "cm"),
            heatmap_height = unit(0.5 * DL.Raj.rows + 8, "cm"),
            cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%d", as.numeric(as.matrix(DL.Raj.count)[i, j])),
                                    x, y,
                                    gp = gpar(fontsize = 10,
                                              col = ifelse(
                                                      is.numeric(as.matrix(DL.Raj.meangt)[i, j]) &
                                                      as.matrix(DL.Raj.meangt)[i, j] < 0.5*ML.Raj.maxmeangt,
                                                      "white", "black")
                                              ))
                       },
            row_title = NULL,
            heatmap_legend_param = list(
                col_fun = DL.Raj.col_fun2, title = "Expression\nLevel",
                legend_width = unit(6, "cm"), 
                legend_height = unit(1, "cm"), 
                direction = "horizontal",
                title_position = "leftcenter"),
            right_annotation = rowAnnotation(
              foo = anno_block(
                  labels = thyset.RFK$set, # SETS,
                  panel_fun = function(index, levels) {
                      # txt = gsub(" ", "\n", thyset.RFK[ rownames(DL.Raj.meangt)[index], "set" ]);
                      txt = thyset.RFK[ rownames(DL.Raj.meangt)[index], "set" ];
                      grid.text(txt, 0, 0.5, rot = 0, just = "left",
                                gp = gpar(col = "black", 
                                          fontsize = 14, fontface = "bold"))
                  },
                  width = unit(6, "cm")
              ))
        );
    pdf(paste0(WD,"/heatmaps/Rajewsky_heatmap_",THYGNS,"_by_celltype_meangt0.pdf"),
        paper = "special",
        width  = unit((1.25 * DL.Raj.cols + 12) / 2.54, "cm"), # 1in~2.54cm
        height = unit((0.55 * DL.Raj.rows + 15) / 2.54, "cm"),
        bg = "white");
    draw(htmp, heatmap_legend_side="top");
    dev.off();
    cat("#-> Rajewsky heatmap by label_short (mean gt0)\n");


    # Heatmaps by LABEL_NUM
    MN.Red <- DT.Red[[THYGNS]] %>%
              select(c(thyset.RFK$geneid,"label_num")) %>%  # select(c(9:49,51)) %>%
              pivot_longer(1:length(thyset.RFK$geneid),names_to="gene",values_to="expression") %>%
              group_by(label_num, gene) %>% 
              summarize(count=sum(expression>0),
                        mean=mean(expression,na.rm=TRUE),
                        meangt=mean(expression[expression>0], na.rm=TRUE));
    MN.Red.count  <- MN.Red %>%
                     select(label_num,gene,count) %>%
                     pivot_wider(names_from = label_num, values_from = count ) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    MN.Red.mean   <- MN.Red %>%
                     select(label_num,gene,mean) %>%
                     pivot_wider(names_from = label_num, values_from = mean  ) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    MN.Red.maxmean <- max(MN.Red$mean,     na.rm=TRUE);
    MN.Red.meangt <- MN.Red %>%
                     select(label_num,gene,meangt) %>%
                     pivot_wider(names_from = label_num, values_from = meangt) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    MN.Red.maxmeangt <- max(MN.Red$meangt, na.rm=TRUE);
    DN.Red.count  <- as.data.frame(MN.Red.count);
    rownames(DN.Red.count)  <- DN.Red.count$gene;
    DN.Red.count  <- DN.Red.count[ thyset.RFK$geneid, ];
    DN.Red.mean   <- as.data.frame(MN.Red.mean);
    rownames(DN.Red.mean)   <- DN.Red.mean$gene;
    DN.Red.mean   <- DN.Red.mean[ thyset.RFK$geneid, ];
    DN.Red.meangt <- as.data.frame(MN.Red.meangt);
    rownames(DN.Red.meangt) <- DN.Red.meangt$gene;
    DN.Red.meangt <- DN.Red.meangt[ thyset.RFK$geneid, ];
    DN.Red.count  <- as.matrix(subset(DN.Red.count,  select = -gene));
    DN.Red.mean   <- as.matrix(subset(DN.Red.mean,   select = -gene));
    DN.Red.meangt <- as.matrix(subset(DN.Red.meangt, select = -gene));
    DN.Red.meangtlg2 <- log2(DN.Red.meangt);
    DN.Red.meangtlg2[is.infinite(DN.Red.meangtlg2)] <- NA

    DN.Red.rows <- nrow(DN.Red.count);
    DN.Red.cols <- ncol(DN.Red.count);
    DN.Red.col_fun <- colorRamp2(seq(0,MN.Red.maxmean,MN.Red.maxmean/24), HMPcols);
    htmp <- Heatmap(
            DN.Red.mean,
            na_col = "#CCCCCC",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            cluster_row_slices = FALSE,
            column_order = Red.CELL.ordnum[ Red.CELL.ordnum %in% colnames(DN.Red.mean) ], # Red.CELL.ordnum[1:length(colnames(DN.Red.mean))], # colnames(DN.Red.mean),
            row_order = rownames(DN.Red.mean), # thyset.RFK$geneid, # 
            row_labels = sprintf("%25s", thyset.RFK[ rownames(DN.Red.mean), "symbol" ]),
            row_names_side = "left",
            row_split = factor(thyset.RFK[ rownames(DN.Red.mean), "set" ], levels=SETS),
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            col = DN.Red.col_fun,
            column_names_rot = 45,
            column_names_gp = gpar(fontsize = 15, fontface = "bold"),
            heatmap_width  = unit(1.25 * DN.Red.cols + 12,  "cm"),
            heatmap_height = unit(0.5 * DN.Red.rows + 8, "cm"),
            cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%d", as.numeric(as.matrix(DN.Red.count)[i, j])),
                                    x, y,
                                    gp = gpar(fontsize = 10,
                                              col = ifelse(
                                                      is.numeric(as.matrix(DN.Red.mean)[i, j]) &
                                                      as.matrix(DN.Red.mean)[i, j] < 0.5*MN.Red.maxmean,
                                                      "white", "black")
                                              ))
                       },
            row_title = NULL,
            heatmap_legend_param = list(
                col_fun = DN.Red.col_fun, title = "Expression\nLevel",
                legend_width = unit(6, "cm"), 
                legend_height = unit(1, "cm"), 
                direction = "horizontal",
                title_position = "leftcenter"),
            right_annotation = rowAnnotation(
              foo = anno_block(
                  labels = thyset.RFK$set, # SETS,
                  panel_fun = function(index, levels) {
                      # txt = gsub(" ", "\n", thyset.RFK[ rownames(DN.Red.mean)[index], "set" ]);
                      txt = thyset.RFK[ rownames(DN.Red.mean)[index], "set" ];
                      grid.text(txt, 0, 0.5, rot = 0, just = "left",
                                gp = gpar(col = "black", 
                                          fontsize = 14, fontface = "bold"))
                  },
                  width = unit(6, "cm")
              ))
        );
    pdf(paste0(WD,"/heatmaps/Reddien_heatmap_",THYGNS,"_by_cellnum.pdf"),
        paper = "special",
        width  = unit((1.25 * DN.Red.cols + 12) / 2.54, "cm"), # 1in~2.54cm
        height = unit((0.55 * DN.Red.rows + 10) / 2.54, "cm"),
        bg = "white");
    draw(htmp, heatmap_legend_side="top");
    dev.off();
    cat("#-> Reddien heatmap by label_num\n");

    DN.Red.col_fun2 <- colorRamp2(seq(1,MN.Red.maxmeangt,(MN.Red.maxmeangt-1)/24), HMPcols);
    htmp <- Heatmap(
            DN.Red.meangt,
            na_col = "#CCCCCC",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            cluster_row_slices = FALSE,
            column_order = Red.CELL.ordnum[ Red.CELL.ordnum %in% colnames(DN.Red.meangt) ], # Red.CELL.ordnum[1:length(colnames(DN.Red.meangt))], # colnames(DN.Red.meangt),
            row_order = rownames(DN.Red.meangt), # thyset.RFK$geneid, # 
            row_labels = sprintf("%25s", thyset.RFK[ rownames(DN.Red.meangt), "symbol" ]),
            row_names_side = "left",
            row_split = factor(thyset.RFK[ rownames(DN.Red.meangt), "set" ], levels=SETS),
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            col = DN.Red.col_fun2,
            column_names_rot = 45,
            column_names_gp = gpar(fontsize = 15, fontface = "bold"),
            heatmap_width  = unit(1.25 * DN.Red.cols + 12,  "cm"),
            heatmap_height = unit(0.5 * DN.Red.rows + 8, "cm"),
            cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%d", as.numeric(as.matrix(DN.Red.count)[i, j])),
                                    x, y,
                                    gp = gpar(fontsize = 10,
                                              col = ifelse(
                                                      is.numeric(as.matrix(DN.Red.meangt)[i, j]) &
                                                      as.matrix(DN.Red.meangt)[i, j] < 0.5*MN.Red.maxmeangt,
                                                      "white", "black")
                                              ))
                       },
            row_title = NULL,
            heatmap_legend_param = list(
                col_fun = DN.Red.col_fun2, title = "Expression\nLevel",
                legend_width = unit(6, "cm"), 
                legend_height = unit(1, "cm"), 
                direction = "horizontal",
                title_position = "leftcenter"),
            right_annotation = rowAnnotation(
              foo = anno_block(
                  labels = thyset.RFK$set, # SETS,
                  panel_fun = function(index, levels) {
                      # txt = gsub(" ", "\n", thyset.RFK[ rownames(DN.Red.meangt)[index], "set" ]);
                      txt = thyset.RFK[ rownames(DN.Red.meangt)[index], "set" ];
                      grid.text(txt, 0, 0.5, rot = 0, just = "left",
                                gp = gpar(col = "black", 
                                          fontsize = 14, fontface = "bold"))
                  },
                  width = unit(6, "cm")
              ))
        );
    pdf(paste0(WD,"/heatmaps/Reddien_heatmap_",THYGNS,"_by_cellnum_meangt0.pdf"),
        paper = "special",
        width  = unit((1.25 * DN.Red.cols + 12) / 2.54, "cm"), # 1in~2.54cm
        height = unit((0.55 * DN.Red.rows + 10) / 2.54, "cm"),
        bg = "white");
    draw(htmp, heatmap_legend_side="top");
    dev.off();
    cat("#-> Reddien heatmap by label_num (mean gt0)\n");


    MN.Raj <- DT.Raj[[THYGNS]] %>%
              select(c(thyset.RFK$geneid,"label_num")) %>%  # select(c(9:49,52)) %>%
              pivot_longer(1:length(thyset.RFK$geneid),names_to="gene",values_to="expression") %>%
              group_by(label_num, gene) %>% 
              summarize(count=sum(expression>0),
                        mean=mean(expression,na.rm=TRUE),
                        meangt=mean(expression[expression>0], na.rm=TRUE));
    MN.Raj.count  <- MN.Raj %>%
                     select(label_num,gene,count) %>%
                     pivot_wider(names_from = label_num, values_from = count ) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    MN.Raj.mean   <- MN.Raj %>%
                     select(label_num,gene,mean) %>%
                     pivot_wider(names_from = label_num, values_from = mean  ) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    MN.Raj.maxmean <- max(MN.Raj$mean,     na.rm=TRUE);
    MN.Raj.meangt <- MN.Raj %>%
                     select(label_num,gene,meangt) %>%
                     pivot_wider(names_from = label_num, values_from = meangt) %>%
                     mutate(across(where(is.numeric), ~replace_na(., 0)));
    MN.Raj.maxmeangt <- max(MN.Raj$meangt, na.rm=TRUE);
    DN.Raj.count  <- as.data.frame(MN.Raj.count);
    rownames(DN.Raj.count)  <- DN.Raj.count$gene;
    DN.Raj.count  <- DN.Raj.count[ thyset.RFK$geneid, ];
    DN.Raj.mean   <- as.data.frame(MN.Raj.mean);
    rownames(DN.Raj.mean)   <- DN.Raj.mean$gene;
    DN.Raj.mean   <- DN.Raj.mean[ thyset.RFK$geneid, ];
    DN.Raj.meangt <- as.data.frame(MN.Raj.meangt);
    rownames(DN.Raj.meangt) <- DN.Raj.meangt$gene;
    DN.Raj.meangt <- DN.Raj.meangt[ thyset.RFK$geneid, ];
    DN.Raj.count  <- as.matrix(subset(DN.Raj.count,  select = -gene));
    DN.Raj.mean   <- as.matrix(subset(DN.Raj.mean,   select = -gene));
    DN.Raj.meangt <- as.matrix(subset(DN.Raj.meangt, select = -gene));
    DN.Raj.meangtlg2 <- log2(DN.Raj.meangt);
    DN.Raj.meangtlg2[is.infinite(DN.Raj.meangtlg2)] <- NA

    DN.Raj.rows <- nrow(DN.Raj.count);
    DN.Raj.cols <- ncol(DN.Raj.count);
    DN.Raj.col_fun <- colorRamp2(seq(0,MN.Raj.maxmean,MN.Raj.maxmean/24), HMPcols);
    htmp <- Heatmap(
            DN.Raj.mean,
            na_col = "#CCCCCC",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            cluster_row_slices = FALSE,
            column_order = Raj.CELL.ordnum[ Raj.CELL.ordnum %in% colnames(DN.Raj.mean) ], # Raj.CELL.ordnum[1:length(colnames(DN.Raj.mean))], # colnames(DN.Raj.mean),
            row_order = rownames(DN.Raj.mean), # thyset.RFK$geneid, # 
            row_labels = sprintf("%25s", thyset.RFK[ rownames(DN.Raj.mean), "symbol" ]),
            row_names_side = "left",
            row_split = factor(thyset.RFK[ rownames(DN.Raj.mean), "set" ], levels=SETS),
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            col = DN.Raj.col_fun,
            column_names_rot = 70,
            column_names_gp = gpar(fontsize = 15, fontface = "bold"),
            heatmap_width  = unit(1.25 * DN.Raj.cols + 12,  "cm"),
            heatmap_height = unit(0.5 * DN.Raj.rows + 8, "cm"),
            cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%d", as.numeric(as.matrix(DN.Raj.count)[i, j])),
                                    x, y,
                                    gp = gpar(fontsize = 10,
                                              col = ifelse(
                                                      is.numeric(as.matrix(DN.Raj.mean)[i, j]) &
                                                      as.matrix(DN.Raj.mean)[i, j] < 0.5*MN.Raj.maxmean,
                                                      "white", "black")
                                              ))
                       },
            row_title = NULL,
            heatmap_legend_param = list(
                col_fun = DN.Raj.col_fun, title = "Expression\nLevel",
                legend_width = unit(6, "cm"), 
                legend_height = unit(1, "cm"), 
                direction = "horizontal",
                title_position = "leftcenter"),
            right_annotation = rowAnnotation(
              foo = anno_block(
                  labels = thyset.RFK$set, # SETS,
                  panel_fun = function(index, levels) {
                      # txt = gsub(" ", "\n", thyset.RFK[ rownames(DN.Raj.mean)[index], "set" ]);
                      txt = thyset.RFK[ rownames(DN.Raj.mean)[index], "set" ];
                      grid.text(txt, 0, 0.5, rot = 0, just = "left",
                                gp = gpar(col = "black", 
                                          fontsize = 14, fontface = "bold"))
                  },
                  width = unit(6, "cm")
              ))
        );
    pdf(paste0(WD,"/heatmaps/Rajewsky_heatmap_",THYGNS,"_by_cellnum.pdf"),
        paper = "special",
        width  = unit((1.25 * DN.Raj.cols + 12) / 2.54, "cm"), # 1in~2.54cm
        height = unit((0.55 * DN.Raj.rows + 16) / 2.54, "cm"),
        bg = "white");
    draw(htmp, heatmap_legend_side="top");
    dev.off();
    cat("#-> Rajewsky heatmap by label_num\n");

    DN.Raj.col_fun2 <- colorRamp2(seq(1,MN.Raj.maxmeangt,(MN.Raj.maxmeangt-1)/24), HMPcols);
    htmp <- Heatmap(
            DN.Raj.meangt,
            na_col = "#CCCCCC",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            cluster_row_slices = FALSE,
            column_order = Raj.CELL.ordnum[ Raj.CELL.ordnum %in% colnames(DN.Raj.meangt) ], # Raj.CELL.ordnum[1:length(colnames(DN.Raj.meangt))], # colnames(DN.Raj.meangt),
            row_order = rownames(DN.Raj.meangt), # thyset.RFK$geneid, # 
            row_labels = sprintf("%25s", thyset.RFK[ rownames(DN.Raj.meangt), "symbol" ]),
            row_names_side = "left",
            row_split = factor(thyset.RFK[ rownames(DN.Raj.meangt), "set" ], levels=SETS),
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            col = DN.Raj.col_fun2,
            column_names_rot = 70,
            column_names_gp = gpar(fontsize = 15, fontface = "bold"),
            heatmap_width  = unit(1.25 * DN.Raj.cols + 12,  "cm"),
            heatmap_height = unit(0.5 * DN.Raj.rows + 8, "cm"),
            cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%d", as.numeric(as.matrix(DN.Raj.count)[i, j])),
                                    x, y,
                                    gp = gpar(fontsize = 10,
                                              col = ifelse(
                                                      is.numeric(as.matrix(DN.Raj.meangt)[i, j]) &
                                                      as.matrix(DN.Raj.meangt)[i, j] < 0.5*MN.Raj.maxmeangt,
                                                      "white", "black")
                                              ))
                       },
            row_title = NULL,
            heatmap_legend_param = list(
                col_fun = DN.Raj.col_fun2, title = "Expression\nLevel",
                legend_width = unit(6, "cm"), 
                legend_height = unit(1, "cm"), 
                direction = "horizontal",
                title_position = "leftcenter"),
            right_annotation = rowAnnotation(
              foo = anno_block(
                  labels = thyset.RFK$set, # SETS,
                  panel_fun = function(index, levels) {
                      # txt = gsub(" ", "\n", thyset.RFK[ rownames(DN.Raj.meangt)[index], "set" ]);
                      txt = thyset.RFK[ rownames(DN.Raj.meangt)[index], "set" ];
                      grid.text(txt, 0, 0.5, rot = 0, just = "left",
                                gp = gpar(col = "black", 
                                          fontsize = 14, fontface = "bold"))
                  },
                  width = unit(6, "cm")
              ))
        );
    pdf(paste0(WD,"/heatmaps/Rajewsky_heatmap_",THYGNS,"_by_cellnum_meangt0.pdf"),
        paper = "special",
        width  = unit((1.25 * DN.Raj.cols + 12) / 2.54, "cm"), # 1in~2.54cm
        height = unit((0.55 * DN.Raj.rows + 16) / 2.54, "cm"),
        bg = "white");
    draw(htmp, heatmap_legend_side="top");
    dev.off();
    cat("#-> Rajewsky heatmap by label_num (mean gt0)\n");


} # make_all_heatmaps