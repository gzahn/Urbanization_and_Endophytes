googlemap_json_to_string <- 
  function (style_list) 
{
  style_string <- ""
  for (i in 1:length(style_list)) {
    if ("featureType" %in% names(style_list[[i]])) {
      style_string <- paste0(style_string, "feature:", 
                             style_list[[i]]$featureType, "|")
    }
    if ("elementType" %in% names(style_list[[i]])) {
      style_string <- paste0(style_string, "element:", 
                             style_list[[i]]$elementType, "|")
    }
    elements <- style_list[[i]]$stylers
    a <- lapply(elements, function(x) paste0(names(x), ":", 
                                             x)) %>% unlist() %>% paste0(collapse = "|")
    style_string <- paste0(style_string, a)
    if (i < length(style_list)) {
      style_string <- paste0(style_string, "&style=")
    }
  }
  style_string <- gsub("#", "0x", style_string)
  style_string <- gsub("[|]", "%7C", style_string)
  return(style_string)
}