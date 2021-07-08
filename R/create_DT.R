

#' Create an interactive data table with download buttons 
create_DT <- function (DF,
                       caption = "",
                       scrollY = 400) 
{
    data <- DT::datatable(DF, caption = caption, extensions = "Buttons", 
                          options = list(dom = "Bfrtip", buttons = c("copy", "csv", 
                                                                     "excel", "pdf", "print"), scrollY = scrollY, scrollX = T, 
                                         scrollCollapse = T, paging = F, columnDefs = list(list(className = "dt-center", 
                                                                                                targets = "_all"))))
    return(data)
}

