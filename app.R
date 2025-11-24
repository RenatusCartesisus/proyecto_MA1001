library(shiny)
library(bslib)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)

#TEMA
###########################################################################

theme = bs_theme(
  bg = "#222222",
  fg = "#F2F5E9",
  primary = "#C8F739",
  baseFont = font_google("Space Mono"),
  codeFont = font_google("Space Mono")
)
primary_color = "#C8F739"
fg_color = "#F2F5E9"
bg_color = "#2F3030"

#NODOS DE LA RED
###########################################################################

calculateNodePropensities = function(G, beta, miu, nu) {
  Lambda = 0
  N = vcount(G)
  
  for (i in 1:N) {
    state_i = V(G)$state[i]
    
    if (state_i == 0) {
      # Susceptible: puede infectarse por vecinos infectados, o vacunarse a tasa nu
      neighbors = neighbors(G, i)
      neighbor_states = V(G)$state[neighbors]
      infection_rate = beta * sum(neighbor_states == 1)
      vacc_rate = nu
      lambda_i = infection_rate + vacc_rate
      V(G)$lambda[i] = lambda_i
      V(G)$inf_rate[i] = infection_rate
      V(G)$vac_rate[i] = vacc_rate
      Lambda = Lambda + lambda_i
      
    } else if (state_i == 1) {
      # Infectados: recuperacion a tasa miu
      V(G)$lambda[i] = miu
      V(G)$inf_rate[i] = 0
      V(G)$vac_rate[i] = 0
      Lambda = Lambda + miu
      
    } else {
      # Recuperados (2) o Vacunados (3): no generan eventos
      V(G)$lambda[i] = 0
      V(G)$inf_rate[i] = 0
      V(G)$vac_rate[i] = 0
    }
  }
  
  return(list(G = G, Lambda = Lambda))
}

###########################################################################

#PROCESO ESTOCASTICO
###########################################################################

draw_next_event_direct = function(Lambda, G) {
  u1 = runif(1)
  u2 = runif(1)
  tau = -log(1 - u1) / Lambda
  target_sum = u2 * Lambda
  sum_i = 0
  i_selected = NA
  N = vcount(G)
  
  for (i in 1:N) {
    sum_i = sum_i + V(G)$lambda[i]
    if (sum_i >= target_sum) {
      i_selected = i
      break
    }
  }
  
  return(list(i = i_selected, tau = tau))
}

###########################################################################

#ACTUALIZAR LOS ESTADOS (ahora con vacunacion)
###########################################################################

update_states = function(G, X, Lambda, i, beta, miu, nu) {
  # X: vector de conteos [S, I, R, V]
  state_before = V(G)$state[i]
  
  if (state_before == 0) {
    # Si es susceptible, puede pasar a infectado o vacunado.
    infection_rate = V(G)$inf_rate[i]
    vacc_rate = V(G)$vac_rate[i]
    total = infection_rate + vacc_rate
    
    if (total <= 0) {
      # sin propensiones, no deberia ocurrir, pero lo manejamos
      return(list(X = X, Lambda = Lambda, G = G))
    }
    
    # decidir tipo de evento proporcional a las tasas
    if (runif(1) * total < infection_rate) {
      # S -> I
      X[1] = X[1] - 1
      X[2] = X[2] + 1
      V(G)$state[i] = 1
      
      old_lambda = V(G)$lambda[i]
      Lambda = Lambda - old_lambda + miu
      V(G)$lambda[i] = miu
      
      # cuando un nodo se infecta, aumenta la tasa de sus vecinos susceptibles
      neighbors_i = neighbors(G, i)
      susceptible_neighbors = neighbors_i[V(G)$state[neighbors_i] == 0]
      
      if (length(susceptible_neighbors) > 0) {
        for (j in susceptible_neighbors) {
          V(G)$lambda[j] = V(G)$lambda[j] + beta
          V(G)$inf_rate[j] = V(G)$inf_rate[j] + beta
          Lambda = Lambda + beta
        }
      }
      
    } else {
      # S -> V (vacunacion perfecta)
      X[1] = X[1] - 1
      X[4] = X[4] + 1
      V(G)$state[i] = 3
      
      old_lambda = V(G)$lambda[i]
      Lambda = Lambda - old_lambda
      V(G)$lambda[i] = 0
      V(G)$inf_rate[i] = 0
      V(G)$vac_rate[i] = 0
      
      # No se ajustan vecinos porque nodo era susceptible previamente (no aportaba beta)
    }
    
  } else if (state_before == 1) {
    # I -> R
    X[2] = X[2] - 1
    X[3] = X[3] + 1
    V(G)$state[i] = 2
    Lambda = Lambda - miu
    V(G)$lambda[i] = 0
    
    # Al recuperarse, sus vecinos susceptibles pierden la contribucion de beta
    neighbors_i = neighbors(G, i)
    susceptible_neighbors = neighbors_i[V(G)$state[neighbors_i] == 0]
    
    if (length(susceptible_neighbors) > 0) {
      for (j in susceptible_neighbors) {
        V(G)$lambda[j] = V(G)$lambda[j] - beta
        V(G)$inf_rate[j] = V(G)$inf_rate[j] - beta
        Lambda = Lambda - beta
      }
    }
  }
  
  return(list(X = X, Lambda = Lambda, G = G))
}

###########################################################################

#METODO DIRECTO SIRV
###########################################################################

direct_method_SIRV_graph = function(G, beta, miu, nu, T) {
  N = vcount(G)
  S = sum(V(G)$state == 0)
  I = sum(V(G)$state == 1)
  R = sum(V(G)$state == 2)
  Vv = sum(V(G)$state == 3)
  
  result = calculateNodePropensities(G, beta, miu, nu)
  G = result$G
  Lambda = result$Lambda
  t = 0
  X_t = data.frame(t = t, S = S, I = I, R = R, V = Vv)
  
  while (t < T) {
    if (Lambda < 1e-10) {
      X_t = rbind(X_t, c(T, S, I, R, Vv))
      break
    }
    
    event = draw_next_event_direct(Lambda, G)
    i = event$i
    tau = event$tau
    t = t + tau
    
    update_result = update_states(G, c(S, I, R, Vv), Lambda, i, beta, miu, nu)
    S = update_result$X[1]
    I = update_result$X[2]
    R = update_result$X[3]
    Vv = update_result$X[4]
    Lambda = update_result$Lambda
    G = update_result$G
    
    X_t = rbind(X_t, c(t, S, I, R, Vv))
  }
  
  return(X_t)
}

###########################################################################

#TEMA PERSONALIZADO (igual que antes)
###########################################################################

theme_custom = function() {
  theme_minimal(base_size = 14, base_family = "Space Mono") +
    theme(
      plot.background = element_rect(fill = "#222222", color = NA),
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


#INTERFAZ (se agrega nu)
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
                                "; color: #222222; border: none;",
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

#SERVIDOR
###########################################################################

server = function(input, output, session) {
  
  simulation_data = reactiveVal(NULL)
  network_graph = reactiveVal(NULL)
  
  observeEvent(input$simulate, {
    withProgress(message = 'Ejecutando simulaciones...', value = 0, {
      set.seed(input$seed)
      
      if (input$network_type == "regular") {
        G0 = sample_k_regular(input$N, input$k)
      } else if (input$network_type == "erdos_renyi") {
        m = round(input$k * input$N / 2)
        G0 = sample_gnm(input$N, m)
      } else {
        G0 = sample_pa(input$N, m = floor(input$k / 2))
      }
      
      network_graph(G0)
      X_array = list()
      
      for (q in 1:input$n_sims) {
        incProgress(1/input$n_sims, detail = paste("Simulación", q))
        G = G0
        V(G)$state = 0
        V(G)$lambda = 0
        V(G)$inf_rate = 0
        V(G)$vac_rate = 0
        
        initial_infected = sample(1:input$N, input$initial_infected)
        V(G)$state[initial_infected] = 1
        
        X_t = direct_method_SIRV_graph(G, input$beta, input$miu, input$nu, input$T)
        X_t$simulation = q
        X_array[[q]] = X_t
      }
      
      all_sims = bind_rows(X_array)
      simulation_data(all_sims)
    })
  })
  
  output$network_plot = renderPlot({
    req(network_graph())
    
    G = network_graph()
    
    # Si algún nodo no tiene estado asignado (solo aparece al inicio)
    if (is.null(V(G)$state)) {
      V(G)$state = 0
      initial_infected = sample(1:vcount(G), input$initial_infected)
      V(G)$state[initial_infected] = 1
    }
    
    # Colores por estado:
    # 0 = S, 1 = I, 2 = R, 3 = V
    colors_map = c(
      "0" = "#4A9EFF",         # S  
      "1" = "#FF6B6B",         # I  
      "2" = primary_color,     # R  
      "3" = "#C77DFF"          # V  
    )
    
    node_colors = colors_map[as.character(V(G)$state)]
    
    par(bg = "#222222", col = fg_color, col.axis = fg_color, 
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
           col = c("#4A9EFF", "#FF6B6B", primary_color, "#C77DFF"),
           pch = 19,
           pt.cex = 2,
           text.col = fg_color,
           bg = bg_color,
           box.col = "#3A3A3A",
           cex = 1.1)
  }, bg = "#222222")
  
  
  output$network_plot = renderPlot({
    req(network_graph())
    
    G = network_graph()
    V(G)$state = 0
    initial_infected = sample(1:vcount(G), input$initial_infected)
    V(G)$state[initial_infected] = 1
    
    colors = c("#4A9EFF", "#FF6B6B", primary_color, "#7BD389")[V(G)$state + 1]
    
    par(bg = "#222222", col = fg_color, col.axis = fg_color, 
        col.lab = fg_color, col.main = primary_color, family = "sans")
    
    plot(G,
         vertex.color = colors,
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
           legend = c("Susceptible", "Infectado"),
           col = c("#4A9EFF", "#FF6B6B"),
           pch = 19,
           pt.cex = 2,
           text.col = fg_color,
           bg = bg_color,
           box.col = "#3A3A3A",
           cex = 1.1)
  }, bg = "#222222")
  
  output$stats_plot = renderPlot({
    req(simulation_data())
    
    stats = simulation_data() %>%
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
        values = c("S" = "#4A9EFF", "I" = "#FF6B6B", "R" = primary_color, "V" = "#7BD389"),
        labels = c("Susceptibles", "Infectados", "Recuperados", "Vacunados")
      ) +
      scale_fill_manual(
        values = c("S" = "#4A9EFF", "I" = "#FF6B6B", "R" = primary_color, "V" = "#7BD389"),
        labels = c("Susceptibles", "Infectados", "Recuperados", "Vacunados")
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
  }, bg = "#222222")
  
  output$r0_display = renderUI({
    R0 = input$beta * input$k / input$miu
    color = if(R0 > 1) "#FF6B6B" else primary_color
    
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
                       if(R0 > 1) "#FF6B6B" else primary_color, ";"),
        if(R0 > 1) "Epidemia esperada" else "✓ Extinción esperada"
      )
    )
  })
  
  output$network_info_display = renderUI({
    req(network_graph())
    G = network_graph()
    
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


#CORRER APP
###########################################################################

shinyApp(ui = ui, server = server)

