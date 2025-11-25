library(shiny)
library(bslib)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)

# TEMA Y COLORES
###########################################################################

# Paleta de colores profesional para un diseño moderno y oscuro
bg_main_color = "#1A1A2E"      # Azul oscuro imperial (fondo principal de la app)
bg_sidebar_color = "#16213E"   # Azul más oscuro para la barra lateral
fg_text_color = "#EAEAEA"       # Blanco humo para texto general
primary_accent_color = "#00BFFF" # Azul cielo profundo para acentos y títulos
color_s = "#4A9EFF"            # Azul para susceptibles
color_i = "#FF6B6B"            # Rojo para infectados
color_r = "#39F7A1"            # Verde menta para recuperados
color_v = "#C77DFF"            # Morado para vacunados

theme = bs_theme(
  bg = bg_main_color,
  fg = fg_text_color,
  primary = primary_accent_color,
  base_font = font_google("Space Mono"),
  code_font = font_google("Space Mono")
)

# Mapeo a nombres de variables existentes para compatibilidad
primary_color = primary_accent_color
fg_color = fg_text_color
bg_color = bg_sidebar_color


###########################################################################
# FUNCIONES OPTIMIZADAS Y PRECOMPUTACION (Trabajan sobre vectores)
###########################################################################

# calculateNodePropensities: recibe vector estado y lista de vecinos
calculateNodePropensities <- function(state, neighbors_list, beta, miu, nu) {
  N <- length(state)
  inf_rate <- numeric(N)
  vac_rate <- numeric(N)
  
  # Para susceptibles: beta * (# vecinos infectados) + nu
  sus_idx <- which(state == 0)
  if (length(sus_idx) > 0) {
    # contar infectados en vecinos
    infected_neighbors_count <- vapply(neighbors_list[sus_idx],
                                       function(nb) sum(state[nb] == 1),
                                       integer(1))
    inf_rate[sus_idx] <- beta * infected_neighbors_count
    vac_rate[sus_idx] <- nu
  }
  
  # Infectados: tasa de recuperación = miu
  lambda <- numeric(N)
  lambda[state == 0] <- inf_rate[state == 0] + vac_rate[state == 0]
  lambda[state == 1] <- miu
  # lambda for recovered/vaccinated stays 0
  
  list(lambda = lambda, inf_rate = inf_rate, vac_rate = vac_rate,
       Lambda = sum(lambda))
}

# draw_next_event_direct: retorna índice y tau; usa cumsum sobre lambda
draw_next_event_direct <- function(lambda, Lambda) {
  # Si Lambda 0, no hay evento
  if (Lambda <= 0) return(list(i = NA_integer_, tau = Inf))
  u1 <- runif(1); u2 <- runif(1)
  tau <- -log(1 - u1) / Lambda
  target <- u2 * Lambda
  cum <- cumsum(lambda)
  # Encontrar primer índice donde cum >= target
  i <- which(cum >= target)[1]
  list(i = i, tau = tau)
}

# update_states: trabaja con vectores para eficiencia; devuelve vectores actualizados y Lambda
update_states <- function(state, lambda, inf_rate, vac_rate, X, i, neighbors_i, beta, miu, nu) {
  # state: integer vector; i: índice seleccionado
  s_before <- state[i]
  N <- length(state)
  
  if (is.na(s_before)) return(list(state = state, lambda = lambda,
                                   inf_rate = inf_rate, vac_rate = vac_rate, X = X))
  # S -> I o S -> V
  if (s_before == 0) {
    infection_rate <- inf_rate[i]
    vacc_rate <- vac_rate[i]
    total <- infection_rate + vacc_rate
    if (total <= 0) {
      return(list(state = state, lambda = lambda, inf_rate = inf_rate,
                  vac_rate = vac_rate, X = X))
    }
    if (runif(1) * total < infection_rate) {
      # S -> I
      state[i] <- 1
      X["S"] <- X["S"] - 1
      X["I"] <- X["I"] + 1
      
      # actualizar lambda del nodo i
      lambda[i] <- miu
      inf_rate[i] <- 0
      vac_rate[i] <- 0
      
      # vecinos susceptibles ganan beta en su inf_rate & lambda
      if (length(neighbors_i) > 0) {
        sus_nb <- neighbors_i[state[neighbors_i] == 0]
        if (length(sus_nb) > 0) {
          inf_rate[sus_nb] <- inf_rate[sus_nb] + beta
          lambda[sus_nb] <- lambda[sus_nb] + beta
        }
      }
      
    } else {
      # S -> V
      state[i] <- 3
      X["S"] <- X["S"] - 1
      X["V"] <- X["V"] + 1
      
      # quitar tasas
      lambda[i] <- 0
      inf_rate[i] <- 0
      vac_rate[i] <- 0
    }
    
  } else if (s_before == 1) {
    # I -> R
    state[i] <- 2
    X["I"] <- X["I"] - 1
    X["R"] <- X["R"] + 1
    
    # quitar lambda del nodo i
    lambda[i] <- 0
    inf_rate[i] <- 0
    vac_rate[i] <- 0
    
    # vecinos susceptibles pierden beta
    if (length(neighbors_i) > 0) {
      sus_nb <- neighbors_i[state[neighbors_i] == 0]
      if (length(sus_nb) > 0) {
        inf_rate[sus_nb] <- inf_rate[sus_nb] - beta
        lambda[sus_nb] <- lambda[sus_nb] - beta
        # evitar errores numéricos negativos
        inf_rate[sus_nb][inf_rate[sus_nb] < 0] <- 0
        lambda[sus_nb][lambda[sus_nb] < 0] <- 0
      }
    }
  }
  
  # recalcular Lambda (sum lambda)
  Lambda_new <- sum(lambda)
  list(state = state, lambda = lambda, inf_rate = inf_rate,
       vac_rate = vac_rate, X = X, Lambda = Lambda_new)
}

# direct_method_SIRV_graph: ahora recibe neighbors_list y trabaja con vectores
direct_method_SIRV_graph <- function(neighbors_list, state_init, beta, miu, nu, T) {
  N <- length(state_init)
  
  state <- state_init
  X <- c(S = sum(state == 0), I = sum(state == 1),
         R = sum(state == 2), V = sum(state == 3))
  
  # inicializar tasas usando la función vectorizada
  props <- calculateNodePropensities(state, neighbors_list, beta, miu, nu)
  lambda <- props$lambda
  inf_rate <- props$inf_rate
  vac_rate <- props$vac_rate
  Lambda <- props$Lambda
  
  t <- 0
  # guardamos resultados en lista para evitar rbind en cada evento
  rows <- list()
  rows[[1]] <- data.frame(t = 0, S = X["S"], I = X["I"], R = X["R"], V = X["V"])
  idx_row <- 2
  
  while (t < T && Lambda > 1e-10) {
    ev <- draw_next_event_direct(lambda, Lambda)
    i <- ev$i
    tau <- ev$tau
    if (is.na(i)) break
    t <- t + tau
    # vecinos del nodo i (precomputados)
    neighbors_i <- neighbors_list[[i]]
    res <- update_states(state, lambda, inf_rate, vac_rate, X, i, neighbors_i,
                         beta, miu, nu)
    state <- res$state
    lambda <- res$lambda
    inf_rate <- res$inf_rate
    vac_rate <- res$vac_rate
    X <- res$X
    Lambda <- res$Lambda
    
    rows[[idx_row]] <- data.frame(t = t, S = X["S"], I = X["I"], R = X["R"], V = X["V"])
    idx_row <- idx_row + 1
  }
  
  # si se acabaron eventos pero tiempo < T, añadir punto final
  if (t < T) {
    rows[[idx_row]] <- data.frame(t = T, S = X["S"], I = X["I"], R = X["R"], V = X["V"])
  }
  
  bind_rows(rows)
}

###########################################################################
# TEMA VISUAL (igual que tu original; usa variables globales de color)
###########################################################################

theme_custom = function() {
  theme_minimal(base_size = 14, base_family = "Space Mono") +
    theme(
      plot.background = element_rect(fill = bg_main_color, color = NA),
      panel.background = element_rect(fill = bg_color, color = NA),
      panel.grid.major = element_line(color = "#3A3A3A", linewidth = 0.3),
      panel.grid.minor = element_line(color = "#2F3030", linewidth = 0.2),
      text = element_text(color = fg_color),
      axis.text = element_text(color = fg_color),
      axis.title = element_text(color = fg_color, face = "bold"),
      plot.title = element_text(color = primary_color, face = "bold", 
                                hjust = 0.5, size = 16),
      plot.subtitle = element_text(color = fg_color, hjust = 0.5, size = 12),
      legend.background = element_rect(fill = bg_color, color = NA),
      legend.text = element_text(color = fg_color),
      legend.title = element_text(color = fg_color, face = "bold"),
      legend.key = element_rect(fill = bg_color, color = NA),
      panel.border = element_rect(color = "#3A3A3A", fill = NA, linewidth = 1)
    )
}

###########################################################################
# INTERFAZ (UI) - sin cambios funcionales
###########################################################################

ui = page_sidebar(
  theme = theme,
  title = div(
    style = "display: flex; align-items: center; gap: 10px;",
    span("Modelo SIRV en Redes", style = "font-size: 24px;"),
  ),
  
  sidebar = sidebar(
    width = 320,
    bg = bg_color,
    
    card(
      card_header(
        "Parámetros Epidemiológicos",
        style = paste0("background-color: ", bg_color, "; color: ", primary_color, 
                       "; font-weight: bold;")
      ),
      
      sliderInput("beta", 
                  "Tasa de infección (β):",
                  min = 0.01, max = 0.5, value = 0.1, step = 0.01),
      
      sliderInput("miu", 
                  "Tasa de recuperación (μ):",
                  min = 0.01, max = 0.2, value = 0.03, step = 0.01),
      
      sliderInput("nu",
                  "Tasa de vacunación (ν):",
                  min = 0, max = 0.5, value = 0.01, step = 0.01)
    ),
    
    card(
      card_header(
        "Parámetros de Red",
        style = paste0("background-color: ", bg_color, "; color: ", primary_color, 
                       "; font-weight: bold;")
      ),
      
      sliderInput("N", 
                  "Número de nodos (N):",
                  min = 50, max = 500, value = 100, step = 50),
      
      sliderInput("k", 
                  "Grado promedio (k):",
                  min = 2, max = 20, value = 5, step = 1),
      
      selectInput("network_type", 
                  "Tipo de red:",
                  choices = c(
                    "Regular" = "regular",
                    "Erdős-Rényi" = "erdos_renyi",
                    "Barabási-Albert" = "barabasi_albert"
                  ),
                  selected = "regular")
    ),
    
    card(
      card_header(
        "Parámetros de Simulación",
        style = paste0("background-color: ", bg_color, "; color: ", primary_color, 
                       "; font-weight: bold;")
      ),
      
      sliderInput("T", 
                  "Tiempo de simulación (T):",
                  min = 50, max = 500, value = 100, step = 50),
      
      sliderInput("n_sims", 
                  "Número de simulaciones:",
                  min = 1, max = 100, value = 50, step = 1),
      
      sliderInput("initial_infected",
                  "Nodos inicialmente infectados:",
                  min = 1, max = 10, value = 1, step = 1),
      
      numericInput("seed",
                   "Semilla aleatoria:",
                   value = 42, min = 1, max = 10000)
    ),
    
    actionButton("simulate", 
                 "Ejecutar Simulación", 
                 class = "btn-primary",
                 style = paste0("width: 100%; font-weight: bold; ",
                                "background-color: ", primary_color, 
                                "; color: ", bg_sidebar_color, "; border: none;",
                                " font-size: 16px; padding: 12px;")),
    
    hr(style = paste0("border-color: ", primary_color)),
    
    card(
      card_header(
        "Información del Sistema",
        style = paste0("background-color: ", bg_color, "; color: ", primary_color, 
                       "; font-weight: bold;")
      ),
      div(
        style = paste0("color: ", fg_color, "; font-family: 'Space Mono';"),
        uiOutput("r0_display"),
        hr(style = "border-color: #3A3A3A; margin: 10px 0;"),
        uiOutput("network_info_display")
      )
    )
  ),
  
  navset_card_pill(
    nav_panel(
      "Simulaciones",
      plotOutput("sim_plot", height = "650px")
    ),
    nav_panel(
      "Red",
      plotOutput("network_plot", height = "650px")
    ),
    nav_panel(
      "Estadísticas",
      plotOutput("stats_plot", height = "650px")
    )
  )
)

###########################################################################
# SERVIDOR (usa precomputación)
###########################################################################

server = function(input, output, session) {
  
  simulation_data = reactiveVal(NULL)
  network_graph = reactiveVal(NULL)
  
  observeEvent(input$simulate, {
    withProgress(message = 'Ejecutando simulaciones...', value = 0, {
      set.seed(input$seed)
      
      # Construcción de la red base G0
      if (input$network_type == "regular") {
        G0 <- sample_k_regular(input$N, input$k)
      } else if (input$network_type == "erdos_renyi") {
        m <- round(input$k * input$N / 2)
        G0 <- sample_gnm(input$N, m)
      } else {
        G0 <- sample_pa(input$N, m = floor(input$k / 2))
      }
      
      network_graph(G0) # para graficar la red
      Nnodes <- vcount(G0)
      
      # ---------- PRECOMPUTACION ----------
      # Lista de vecinos (vector de enteros por nodo), precomputada una vez
      # usamos as_adj(sparse = FALSE) para acceso rápido; para N grande podrías usar sparse
      adj_mat <- as_adj(G0, sparse = FALSE)
      neighbors_list <- lapply(1:Nnodes, function(i) which(adj_mat[i, ] != 0))
      # -----------------------------------
      
      X_array <- vector("list", input$n_sims)
      
      for (q in seq_len(input$n_sims)) {
        incProgress(1/input$n_sims, detail = paste("Simulación", q))
        
        # Iniciar vectores de estado y tasas (más rápido que manipular el grafo)
        state <- integer(Nnodes) # 0 = S por defecto
        lambda <- numeric(Nnodes)
        inf_rate <- numeric(Nnodes)
        vac_rate <- numeric(Nnodes)
        
        # infectados iniciales
        initial_infected <- sample.int(Nnodes, input$initial_infected)
        state[initial_infected] <- 1
        
        # inicializar propiedades (vectorizado)
        X_t <- direct_method_SIRV_graph(neighbors_list, state,
                                        input$beta, input$miu, input$nu, input$T)
        
        X_t$simulation <- q
        X_array[[q]] <- X_t
      }
      
      all_sims <- bind_rows(X_array)
      simulation_data(all_sims)
    })
  })
  
  ###########################################################################
  # GRAFICOS (igual que antes)
  ###########################################################################
  
  output$network_plot = renderPlot({
    req(network_graph())
    G <- network_graph()
    
    if (is.null(V(G)$state)) {
      V(G)$state <- 0
      initial_infected <- sample(1:vcount(G), input$initial_infected)
      V(G)$state[initial_infected] <- 1
    }
    
    colors_map <- c("0" = color_s, "1" = color_i, "2" = color_r, "3" = color_v)
    node_colors <- colors_map[as.character(V(G)$state)]
    
    par(bg = bg_main_color, col = fg_color, col.axis = fg_color, 
        col.lab = fg_color, col.main = primary_color, family = "sans")
    
    plot(G,
         vertex.color = node_colors,
         vertex.size = 8,
         vertex.label = NA,
         vertex.frame.color = "#3A3A3A",
         edge.color = "#4A4A4A",
         edge.width = 0.8,
         main = paste("Red inicial -", 
                      switch(input$network_type,
                             "regular" = "Regular",
                             "erdos_renyi" = "Erdős-Rényi",
                             "barabasi_albert" = "Barabási-Albert")),
         layout = layout_with_fr(G))
    
    legend("topright",
           legend = c("Susceptible (S)", "Infectado (I)", "Recuperado (R)", "Vacunado (V)"),
           col = c(color_s, color_i, color_r, color_v),
           pch = 19,
           pt.cex = 2,
           text.col = fg_color,
           bg = bg_color,
           box.col = "#3A3A3A",
           cex = 1.1)
  }, bg = bg_main_color)
  
  output$sim_plot = renderPlot({
    req(simulation_data())
    
    sims_long <- simulation_data() %>%
      pivot_longer(cols = c(S, I, R, V),
                   names_to = "Compartimento",
                   values_to = "N")
    
    ggplot(sims_long, aes(x = t, y = N, group = interaction(simulation, Compartimento), color = Compartimento)) +
      geom_line(alpha = 0.3) +
      scale_color_manual(
        values = c("S" = color_s, "I" = color_i, "R" = color_r, "V" = color_v),
        labels = c("S" = "Susceptibles", "I" = "Infectados", "R" = "Recuperados", "V" = "Vacunados")
      ) +
      labs(
        title = "Trayectorias de las Simulaciones Estocásticas",
        subtitle = paste(input$n_sims, "simulaciones individuales"),
        x = "Tiempo",
        y = "Número de individuos",
        color = "Estado"
      ) +
      theme_custom() +
      guides(color = guide_legend(override.aes = list(alpha = 1)))
  }, bg = bg_main_color)
  
  output$stats_plot = renderPlot({
    req(simulation_data())
    
    stats <- simulation_data() %>%
      group_by(t) %>%
      summarise(
        S_mean = mean(S), I_mean = mean(I), R_mean = mean(R), V_mean = mean(V),
        S_sd = sd(S), I_sd = sd(I), R_sd = sd(R), V_sd = sd(V)
      ) %>%
      pivot_longer(cols = c(S_mean, I_mean, R_mean, V_mean),
                   names_to = "Compartimento",
                   values_to = "Media") %>%
      mutate(
        SD = case_when(
          Compartimento == "S_mean" ~ S_sd,
          Compartimento == "I_mean" ~ I_sd,
          Compartimento == "R_mean" ~ R_sd,
          Compartimento == "V_mean" ~ V_sd
        ),
        Compartimento = case_when(
          Compartimento == "S_mean" ~ "S",
          Compartimento == "I_mean" ~ "I",
          Compartimento == "R_mean" ~ "R",
          Compartimento == "V_mean" ~ "V"
        )
      )
    
    ggplot(stats, aes(x = t, y = Media, color = Compartimento, fill = Compartimento)) +
      geom_line(linewidth = 1.8) +
      geom_ribbon(aes(ymin = Media - SD, ymax = Media + SD),
                  alpha = 0.25, color = NA) +
      scale_color_manual(
        values = c("S" = color_s, "I" = color_i, "R" = color_r, "V" = color_v),
        labels = c("S" = "Susceptibles", "I" = "Infectados", "R" = "Recuperados", "V" = "Vacunados")
      ) +
      scale_fill_manual(
        values = c("S" = color_s, "I" = color_i, "R" = color_r, "V" = color_v),
        labels = c("S" = "Susceptibles", "I" = "Infectados", "R" = "Recuperados", "V" = "Vacunados")
      ) +
      labs(
        title = "Media ± Desviación Estándar",
        subtitle = "Agregado de todas las simulaciones",
        x = "Tiempo",
        y = "Número de individuos",
        color = "Estado",
        fill = "Estado"
      ) +
      theme_custom()
  }, bg = bg_main_color)
  
  output$r0_display = renderUI({
    R0 = input$beta * input$k / input$miu
    color = if(R0 > 1) color_i else primary_color
    
    div(
      style = "text-align: center; padding: 10px;",
      div(
        style = "font-size: 14px; color: #999999;",
        "Número Reproductivo Básico"
      ),
      div(
        style = paste0("font-size: 32px; font-weight: bold; color: ", color, ";"),
        paste0("R₀ = ", round(R0, 3))
      ),
      div(
        style = paste0("font-size: 12px; margin-top: 5px; color: ", 
                       if(R0 > 1) color_i else primary_color, ";"),
        if(R0 > 1) "Epidemia esperada" else "✓ Extinción esperada"
      )
    )
  })
  
  output$network_info_display = renderUI({
    req(network_graph())
    G <- network_graph()
    
    div(
      style = "padding: 10px;",
      div(
        style = "display: flex; justify-content: space-between; margin-bottom: 8px;",
        span("Nodos:", style = "color: #999999;"),
        span(vcount(G), style = paste0("color: ", primary_color, "; font-weight: bold;"))
      ),
      div(
        style = "display: flex; justify-content: space-between; margin-bottom: 8px;",
        span("Aristas:", style = "color: #999999;"),
        span(ecount(G), style = paste0("color: ", primary_color, "; font-weight: bold;"))
      ),
      div(
        style = "display: flex; justify-content: space-between;",
        span("Grado promedio:", style = "color: #999999;"),
        span(round(mean(degree(G)), 2), 
             style = paste0("color: ", primary_color, "; font-weight: bold;"))
      )
    )
  })
}

# CORRER APP
shinyApp(ui = ui, server = server)
