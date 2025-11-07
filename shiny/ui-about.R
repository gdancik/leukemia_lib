aboutPage <- list(
  h1("AML-BET and its Developers:"),
  hr(color = "blue", style = "height: 2px"),
  
  h2("Developer and Maintainer"),
  h4("Nathaniel Gauvin, B.S. in Computer Science, Minor in Bioinformatics"), 
  h4("Eastern Connecticut State University Graduate"),
  hr(color = "blue", style = "height: 2px"),
  
  h2("Advisor"),
  h4("Dr. Garrett Dancik, Ph.D. in Bioinformatics, M.S. in Statistics"), 
  h4("Eastern Connecticut State University Professor and Department Chair"),
  hr(color = "blue", style = "height: 2px"),
  
  h2("Description"),
  p("AML-BET was developed as an honors thesis project for Eastern Connecticut
    State University's Honors Scholars program. AML-BET aims to make biomarker
    evaluation as accessible and accurate as possible. It allows researchers
    interested in AML-associated genes to investigate these genes at their 
    convenience, removing the potentially long wait times that comes with 
    collaboration. AML-BET contains nine microarray and RNA-Seq datasets, 
    with all RNA-Seq datasets normalized using TMM. For additional information, 
    please read the thesis paper associated with AML-BET below.", 
    style = "font-size: 18px"),
  hr(color = "blue", style = "height: 2px"),
  h5(a("Thesis Paper", href = "https://github.com/NateGauvin/AML-BET"))
)