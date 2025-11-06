ui <- fluidPage(
  useShinyjs(),
  
  tags$head(tags$style(HTML("
    body {
      background: linear-gradient(to right, #f0f4f7, #d9e2ec);
      font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
      color: #2c3e50;
    }
    h2 {
      font-weight: 500;
      letter-spacing: 0.3px;
      color: #3e3f4a;
      margin-bottom: 20px;
    }
    .sidebarPanel, .mainPanel { margin-top: 10px; }
    .btn-run {
      background-color: #28a745;
      color: white;
      font-weight: 600;
      font-size: 1.05em;
      border: none;
      border-radius: 8px;
      width: 100%;
      padding: 10px 0;
      box-shadow: 0 2px 5px rgba(0,0,0,0.1);
      transition: background 0.25s ease;
    }
    .btn-run:hover { background-color: #218838; }
    .btn-sample, .btn-secondary {
      background-color: #6c757d;
      color: white;
      font-weight: 600;
      border: none;
      border-radius: 6px;
      width: 100%;
      padding: 8px 0;
      box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      transition: background 0.25s ease;
      margin-bottom: 20px;
    }
    .btn-sample:hover, .btn-secondary:hover { background-color: #5a6268; }
    .card {
      background: white;
      border-radius: 12px;
      box-shadow: 0 4px 10px rgba(0,0,0,0.08);
      margin-bottom: 25px;
      overflow: hidden;
    }
    .card-header {
      background: #4a90e2;
      color: white;
      font-weight: 700;
      padding: 10px 15px;
      font-size: 1.1em;
      letter-spacing: 0.3px;
    }
    .card-content {
      padding: 15px;
    }
    table.dataTable th {
      white-space: nowrap !important;
      background-color: #f8f9fa !important;
    }
    .dataTables_wrapper .dataTables_paginate .paginate_button {
      background: #4a90e2 !important;
      color: white !important;
      border-radius: 5px !important;
      margin: 2px;
    }
    .dataTables_wrapper .dataTables_paginate .paginate_button:hover {
      background: #357ABD !important;
      color: white !important;
    }
    .shiny-download-link.btn {
      margin-top: 10px;
      background-color: #007bff;
      color: white;
      font-weight: 600;
      border-radius: 6px;
      border: none;
      padding: 8px 12px;
    }
    .shiny-download-link.btn:hover {
      background-color: #0056b3;
    }
    .form-control, .selectize-input {
      border-radius: 6px !important;
    }
    footer {
      margin-top: 30px;
      border-top: 1px solid #ccc;
      color: gray;
      font-size: 0.9em;
      padding-top: 10px;
      line-height: 1.5;
    }
    a { color: #4a90e2; text-decoration: none; }
    a:hover { text-decoration: underline; }
  "))),
  
  h2("Conserved Human MHC-I Epitopes"),
  
  sidebarLayout(
    sidebarPanel(
      style = "background:rgba(255,255,255,0.85);border-radius:10px;padding:20px;box-shadow:0 3px 6px rgba(0,0,0,0.08);",
      fileInput("fasta", "Upload FASTA:", accept = c(".fasta", ".fa")),
      actionButton("use_sample", "Use Sample Data", class = "btn btn-secondary"),
      textInput("uniprot_id", "UniProt ID (for published epitopes):", ""),
      selectInput("alleles", "Select HLA Alleles:",
                  c("HLA-A*01:01", "HLA-A*02:01", "HLA-B*07:02", "HLA-B*08:01"),
                  "HLA-A*02:01", multiple = TRUE),
      sliderInput("peptide_length", "Peptide Length:", 8, 14, value = c(9, 10)),
      numericInput("min_conserved_len", "Min Conserved Length:", 10, 5, 100),
      checkboxInput("highlight_regions", "Highlight conserved regions (blue)", TRUE),
      checkboxInput("debug_skip_api", "Skip IEDB Predictions", FALSE),
      actionButton("submit", "Run Prediction", class = "btn-run", style = "margin-top:10px;")
    ),
    
    mainPanel(
      uiOutput("sample_info"),
      uiOutput("prediction_prompt"),
      
      div(class = "card",
          div(class = "card-header", "Interactive Overview Plot"),
          div(class = "card-content", withSpinner(plotlyOutput("sequence_plot")))),
      
      div(class = "card",
          div(class = "card-header", "Predicted Epitopes"),
          div(class = "card-content",
              withSpinner(DTOutput("top_table")),
              downloadButton("download_table", "Download Predicted", class = "btn btn-primary"))),
      
      div(class = "card",
          div(class = "card-header", "Published Epitopes (mapped to FASTA consensus)"),
          div(class = "card-content",
              withSpinner(DTOutput("published_table")),
              downloadButton("download_published", "Download Published", class = "btn btn-primary")))
    )
  ),
  
  tags$footer(
    HTML("
      <div>Version 0.7.4 © 2025 LumaVista Bio. See README.md for license details.</div>
      <div style='margin-top:5px;'>
        <strong>Predictions powered by 
          <a href='https://www.iedb.org' target='_blank' rel='noopener noreferrer'>
            IEDB MHC Binding API
          </a>
        </strong><br/>
        Research and educational use only. Please cite IEDB when generating data with this tool.
      </div>
    ")
  )
)
